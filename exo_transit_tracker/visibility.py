from astropy.time import Time
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation
import astroplan as ap
import warnings

def calc_vis(xo, loc, night='nautical', altMin=15, moonSep=30):
    """
        Calculate the visibility of the host star from the given site and a
        list of constraints and determine when the next visible transit occurs.
        
        Inputs
        ------
        
            xo: Exoplanet object
                Exoplanet object containing a time in the next_transit attribute
                
            loc: EarthLocation 
                astropy.coordinates.EarthLocation object for the observer.
                
            night: str (Optional)
                Identifier of the kind of night to be used for checking visibility.
                Must be one of:
                    sunset: Sun 0 deg below horizon
                    civil: Sun 6 deg below horizon
                    nautical: Sun 12 deg below horizon
                    astronomical: Sun 18 deg below horizon
                    
            altMin: float (Optional)
                Minimum altitude (in degrees) to be used for checking visibility.
                Must be in the range [0, 90].
                
            moonSep: float (Optional)
                Minimum angular separation from the Moon to be used for checking
                visibility. Must be greater than 0.
        
        Outputs
        -------
            True, if the current transit in xo.next_transits is visible
            False, if the star is never visible from this source
    """
    
    # Check visibility args
    if night.lower() not in ['sunset', 'civil', 'nautical', 'astronomical']:
        raise ValueError('{} is not a recognized definition of night time.'.format(night))
    if altMin < 0:
        warnings.warn('Minimum altitude cannot be less than 0 deg. Setting to 0 deg.')
        altMin = 0
    elif altMin >= 90:
        raise ValueError('Minimum altitude cannot be >= 90 deg.')
    if moonSep < 0:
        msg = 'Minimum Moon separation cannot be less than 0 deg. '
        msg += 'Setting to 0 deg.'
        warnings.warn(msg)
        moonSep = 0
    
    # Make sure Exoplanet object has a next_transit attribute computed
    if xo.next_transit is None:
        msg = 'Visibility is being evaluated before a transit '
        msg += 'time has been computed. Computing now...'
        warnings.warn(msg)
        xo.ephemeris_to_next_transit()
        
    # Build astroplan constraints list
    Constraints = []
    if night == 'sunset':
        Constraints.append(ap.AtNightConstraint(max_solar_altitude = 0 * u.deg))
    elif night == 'civil':
        Constraints.append(ap.AtNightConstraint.twilight_civil())
    elif night == 'nautical':
        Constraints.append(ap.AtNightConstraint.twilight_nautical())
    else:
        Constraints.append(ap.AtNightConstraint.twilight_astronomical())
    
    Constraints.append(ap.AltitudeConstraint(min=altMin * u.deg))
    Constraints.append(ap.MoonSeparationConstraint(min=moonSep * u.deg))
    
    Observer = ap.Observer(location=loc, name='Observer')
    
    # First check if this source is ever observable from this location
    all_year = Time(['2023-01-01T00:00:00', '2023-12-31T23:59:59'])
    ever_observable = ap.is_observable(Constraints, Observer, xo.st_coords,
                                       time_range=all_year)[0]
    if not ever_observable:
        return False
    
    # Continue to determine the next observable transit. Start by checking the 
    # transit time already determined for the Exoplanet object
    vis = ap.is_observable(Constraints, Observer, xo.st_coords,
                           times=xo.next_transit)[0]

    if vis:
        return True
        
    while not vis:
        xo.ephemeris_to_next_transit(t=xo.next_transit + 1 * u.hour)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            vis = ap.is_observable(Constraints, Observer, xo.st_coords,
                                   times=xo.next_transit)[0]
        
        # Rip-cord for this while loop
        if xo.next_transit.byear > 10000:
            msg = 'No visibility found up to the year 10,000 (!!).'
            msg += ' Stopping analysis here.'
            raise RuntimeError(msg)
    
    return True