import numpy as np
import os
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord, Longitude, Latitude
from astropy.time import Time
import astropy.units as u
import warnings

##################################################
#
# This script defines the Exoplanet class.
# 
# The Exoplanet class contains basic information
# about the planet in question independent of its
# transit and visibility to the Observer.
#
##################################################

class Exoplanet():
	
	shows_transits = False
	alias_list = []
	simbad_known = False
	ra_deg = None
	dec_deg = None
	ra = None
	dec = None
	st_coords = None
	P = None
	T0 = None
	
	def __init__(self, st_name, pl_id, **kwargs):
		"""
			Constructor of the Exoplanet class.
		
			Inputs
			------
				st_name: str
					Name of the host star.
			
				pl_id: str
					Identifier of the planet.
	
		"""
	
		self.st_name = st_name
		self.pl_id = pl_id
		self.pl_name = '{} {}'.format(st_name, pl_id)
	
		# Define attributes with kwargs
		for key, val in kwargs.items():
			setattr(self, key, val)
		
		
		
		
		
	def get_radec_simbad(self):
		"""
			Query the CDS Simbad database for coordinates. Star
			must be recognized by Simbad in order for RA/Dec to 
			be found. 
	
		"""
	
		res = Simbad.query_object(self.st_name)
		if res is None:
			raise Exception('The name {} cannot '.format(self.st_name) + \
							'be matched to a known ojbect in Simbad.')
	
		# Store RA and Dec in useful units
		self.ra = ':'.join(res['RA'][0].split())
		self.dec = ':'.join(res['DEC'][0].split())
		self.st_coords = SkyCoord(self.ra, self.dec, unit=[u.hourangle, u.deg])
		self.ra_deg = self.st_coords.ra.value
		self.dec_deg = self.st_coords.dec.value
			
# 	def ephemeris_to_next_transit(self, P, T0, t=None):
# 		"""
# 			Calculate the timing of the next transit according to a given ephemeris. 
# 	
# 			Inputs
# 			------
# 	
# 				P: float
# 					Orbital period in days.
# 			
# 				T0: float
# 					Transit (conjunction) time in Barycentric Julian Days.
# 			
# 				t: None or time
# 					The reference time for calculating the "next" transit. If None,
# 					the current time will be used.
# 			
# 			Outputs
# 			-------
# 	
# 				astropy.Time object containing the timing of the next transit relative
# 				to t (or now). 
# 		"""
# 		# Get current time or evaluate passed time
# 		if t is not None:
# 			try:
# 				t = Time(t)
# 			except:
# 				warnings.warn('Cannot understand {} as a time.'.format(t) + \
# 							  ' Using the present time instead.')
# 				t = Time.now()
# 		else:
# 			t = Time.now()
# 			
# 		# Determine next transit time
# 		t0 = Time(T0, format='jd', scale='tdb')
# 		epoch = (t - t0).jd // P + 1
# 	
# 		return Time(T0 + P * epoch, format='jd', scale='tdb')
# 		 
# 	