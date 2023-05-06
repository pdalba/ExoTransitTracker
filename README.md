# ExoTransitTracker
Discover ongoing and upcoming exoplanet transit events visible to you.

An exoplanet "transit" occurs when a distant exoplanet crosses in front of its star from your point of view. This is one of the key ways that exoplanets are discovered and studied. 

With this tool, you can input a planet name and the code will query the CDS Simbad Astronomical Database to determine its right ascension and declination. It will then crossmatch that against the NASA Exoplanet Archive in order to determine the known orbital ephemeris of the planet. Then, the next transit is determined automatically. 

You can also specify an observation location (e.g., longitude, latitude, and altitude) and the code will determine the next transit that is visible to that specific location (if any). 

This package relies on astropy and specifically its subpackage astroplan.

Example call:

>>>from exo_transit_tracker import next_transit
>>>next_transit('Kepler-51', 'b').isot
'2023-05-30T00:51:24.480'

Current version = 1.0
