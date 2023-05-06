import numpy as np
import astropy
import astropy.units as u
from astropy.coordinates import EarthLocation
import warnings

from exo_transit_tracker.exoplanet import Exoplanet
from exo_transit_tracker.utils import query_nea
from exo_transit_tracker.visibility import calc_vis

def next_transit(st_name, pl_id, source='NEA', P=None, T0=None, loc=None):
	"""
		Determine the timing of the next transit of a given exoplanet.
	
		Inputs
		------
		
			st_name: str
				Name of the host star. Must be resolvable by the CDS
				Simbad service.
	
			pl_id: str
				Letter or identifier of the exoplanet.
			
			source: str
				Where to draw the exoplanet transit ephemeris from. 
					Options include:
						'NEA': exoplanet data are downloaded from the 
							   NASA Exoplanet Archive and queried.
						None: the transit ephemeris (P and T0) must 
							  be provided.
				
			P: float or NoneType
				The orbital period (in days) for the planet. Must be
				provided if source is None.
			
			T0: float or NoneType
				The transit (or conjunction) time for the planet in 
				Barycentric Julian Days. Must be provided if source 
				is None.
		
			loc: tuple / list or None
				A tuple containing the longitude (in degrees), 
				latitude (in degrees), and elevation (in meters) for 
				the location on Earth to be treated as the observer.
				If provided, the next_transit will be that which is 
				visible to this observer.
					
		Outputs
		-------
	
			nextTransit: astropy.time.Time
				A Time object of the next transit mid-time for this planet, for
				the given (or determined) ephemeris, visible from the given 
				location (if applicable).			
	"""

	# Parse args
	if source not in ['NEA', None]:
		raise ValueError('Source of ephemeris must be "NEA" or None.')
	elif source is None:
		if P is None or T0 is None:
			raise ValueError('Ephmeris must be provided if source is None.')
	
	if loc is not None:
		if len(loc) != 3:
			msg = 'loc must a be 3 item tuple or list containing'
			msg += 'longitude, latitude, and elevation.'
			raise ValueError(msg)
		try:
			loc = EarthLocation.from_geodetic(loc[0], loc[1], height=loc[2])
		except Exception as ee:
			print(ee)
			
	# If source is NEA and ephemeris is given, default to ephemeris
	if P is not None and T0 is not None and source == 'NEA':
		msg = 'Using given ephemeris as opposed to '
		msg += 'querying NEA.'
		warnings.warn(msg)
		source = None
	
	# Create an Exoplanet object and calculate the next transit from now
	xo = Exoplanet(st_name, pl_id, P=P, T0=T0)
	 	
	# CASE: determine transit timing with given ephemeris
	if source is None and loc is None:
		xo.ephemeris_to_next_transit()
		return xo.next_transit
		
	# CASE: determine transit timing with given ephemeris with a location
	elif source is None and loc is not None:
		# Query simbad for hostName
		xo.get_radec_simbad()
		xo.simbad_known = True
		
		# Calculate next transit time
		xo.ephemeris_to_next_transit()
		
		# Calculate visibility
		vis = calc_vis(xo, loc)
		
		if not vis:
			print('Sorry, this star is never visible from this location.' \
				  ' No transits visible.')
			return
		
		return xo.next_transit
		
	# CASE: get ephemeris from NEA
	else:
		# Query simbad for hostName
		xo.get_radec_simbad()
		xo.simbad_known = True
		
		# Query the NEA for all exoplanet RA/Dec
		qryCols = 'pl_name,ra,dec'
		qryTab = 'ps'
		qryConst = {'col': ['default_flag'],
					'oper': ['='],
					'val': ['1']}
		df_coords = query_nea(col=qryCols, tab=qryTab, const=qryConst)
							
		# Cross match with the target system's RA/Dec
		df_coords['dist'] = np.sqrt((df_coords.ra - xo.ra_deg) ** 2 +
								    (df_coords.dec - xo.dec_deg) ** 2)
		matchInd = df_coords.where(df_coords.dist ==
								   min(df_coords.dist))['dist'].dropna().index
		matchName = df_coords.loc[matchInd, 'pl_name'].values[0]
		if xo.st_name not in matchName:
			msg = 'Matched {} on'.format(xo.st_name)
			msg += ' the NEA to {}. '.format(matchName)
			msg += 'Assuming alias and proceeding.'
			warnings.warn(msg)
		
		# Gather all the ephemerides for this planet from the
		# NEA with another query. 
		qryCols = 'pl_name,default_flag,tran_flag,pl_orbper,pl_tranmid'
		qyrTab = 'ps'
		df_pl = query_nea(col=qryCols, tab=qryTab)
		
		# Trim other planets
		df_pl = df_pl.loc[df_pl.where(df_pl['pl_name'] == 
								      matchName)['pl_name'].dropna().index]
		
		
		# Evaluate the results and pick an ephemeris
		df_sub = df_pl.loc[df_pl.where(df_pl['tran_flag'] == 
								       1)['tran_flag'].dropna().index]
		if len(df_sub) == 0:
			print('According to the NEA, {} does not transit. '.format(xo.pl_name))
			print('If you have discovered a transit of this planet, '
				  'please consult your nearest observatory immediately.')
			return 
		xo.shows_transits = True
		
		# Check the default entry
		df_default = df_pl.loc[df_pl.where(df_pl['default_flag'] == 
								       1)['default_flag'].dropna().index]
		if df_default['pl_orbper'].values[0] is not None and \
			df_default['pl_tranmid'].values[0] is not None:
				xo.P = df_default['pl_orbper'].values[0]
				xo.T0 = df_default['pl_tranmid'].values[0]
		else:
			# Use an entry that gives both P and T0 or use multiple entries
			df_other = df_pl.loc[df_pl.where(df_pl['default_flag'] == 
								       0)['default_flag'].dropna().index]
			df_other['full_ephem'] = [False] * len(df_other)
			final_P, final_T0 = None, None
			for i in df_other.index:
				if df_other.loc[i, 'pl_orbper'] is not None and \
					df_other.loc[i, 'pl_tranmid'] is not None:
						xo.P = df_other.loc[i, 'pl_orbper']
						xo.T0 = df_other.loc[i, 'pl_tranmid']
						break
				else:
					if df_other.loc[i, 'pl_orbper'] is not None and final_P is None:
						final_P = df_other.loc[i, 'pl_orbper']
					if df_other.loc[i, 'pl_tranmid'] is not None and final_T0 is None:
						final_T0 = df_other.loc[i, 'pl_tranmid']
			
			if xo.P is None and final_P is not None:
				xo.P = final_P
			if xo.T0 is None and final_T0 is not None:
				xo.T0 = final_T0
			if np.any([xo.P is None, xo.T0 is None]):
				print('Could not extract a complete ephemeris from the NEA.'
					  ' No transits can be reported.')
				return
								
		# Calculate the next transit from now
		xo.ephemeris_to_next_transit()
									
		# If no location, simply return. Otherwise, check visibility
		if loc is None:
			return xo.next_transit
		else:
			# Calculate visibility
			vis = calc_vis(xo, loc)
		
			if not vis:
				print('Sorry, this star is never visible from this location.' \
				  	  ' No transits visible.')
				return
		
			return xo.next_transit
			
			






