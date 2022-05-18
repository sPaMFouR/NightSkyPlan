#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxx----------CALCULATION OF TWILIGHT TIMES AND NIGHT DURATION---------xxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #


# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import glob
import ephem
import numpy as np
import pandas as pd
import easygui as eg
from astropy.io import fits
from astropy import units as u
from datetime import datetime, timedelta
from astropy.coordinates import Angle, SkyCoord
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables
# ------------------------------------------------------------------------------------------------------------------- #
list_dates = 'DateList.dat'
list_telescopes = 'TelescopeList.dat'
datetime_format = '%Y-%m-%dT%H:%M:%S'

TIME_keyword = 'DATE-OBS'
RA_keyword = 'RA'
DEC_keyword = 'DEC'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions to help compute Moon Phase and Moon Angle
# ------------------------------------------------------------------------------------------------------------------- #

def get_timeradec(filename):
    """
    Extract time of observation 'time_obs', Right Ascension 'ra' and Declination 'dec'
    from the header of the file 'filename'.
    Args:
        filename : FITS file whose header has to be appended
    Returns:
        time_obs : Time of observation for which moon angle and moon separation have to be computed
        ra       : Right Ascension of the source for which moon angle has to be computed
        dec      : Declination of the source for which moon angle has to be computed
    """
    hdulist = fits.open(filename)
    header = hdulist[0].header

    box_msg = 'Verify Header Keywords for Time, RA, DEC of Observation'
    box_title = 'Details of Header Keywords'
    field_names = ['TIME', 'RA', 'DEC']
    field_values = [TIME_keyword, RA_keyword, DEC_keyword]
    header_keywords = eg.multenterbox(msg=box_msg, title=box_title, fields=field_names, values=field_values)

    time_obs = header[header_keywords[0]]
    ra = header[header_keywords[1]]
    dec = header[header_keywords[2]]
    
    return time_obs, ra, dec


def get_galacticcoord(ra, dec):
    """
    Computes Galactic Coordinates (l, b) of the source from Right Ascension 'ra' and Declination 'dec'.
    Args:
        ra   : Right Ascension of the source for which galactic coordinates have to be computed
        dec  : Declination of the source for which galactic coordinates have to be computed
    Returns:
        l    : Galactic latitude of the source
        b    : Declination of the source
    """
    ra = Angle(ra + ' hours').degree
    dec = Angle(dec + ' degrees').degree
    coord = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
    galcoord = coord.galactic
    l, b = round(galcoord.l.degree, 2), round(galcoord.b.degree, 2)
    return l, b


def get_moonphaseandangle(time, ra, dec):
    """
    The function computes moon separation for the object ('ra', 'dec') at the time of observation 'time' 
    and the moon phase.
    Args:
        time       : Time of observation for which moon angle and moon separation have to be computed
        ra         : Right Ascension of the source for which moon angle has to be computed
        dec        : Declination of the source for which moon angle has to be computed
    Returns:
        moon_sep   : Moon separation for the object ('ra', 'dec') at the time of observation
        moon_phase : Moon Phase at the time of observation
        altitude   : Altitude of the source at the time of observation
        airmass    : Airmass for the source at the time of observation
        l          : Galactic longitude of the source at the time of observation
        b          : Galactic latitude for the source at the time of observation
    """
    # Declare the Object of Observation
    object_obs = ephem.FixedBody()
    object_obs._epoch = ephem.J2000
    object_obs._ra = ra
    object_obs._dec = dec

    time_obs = datetime.strptime(time, datetime_format)
    if not utc:
        time_obs += timedelta(hours=OBS_TIMEZONE)

    # Calculation of Moon Phase and Moon Angle
    telescope.date = time_obs
    object_obs.compute(telescope)
    moon = ephem.Moon(time_obs)
    moon_phase = '{0:.2f}'.format(moon.phase)
    moon_seprad = ephem.separation(object_obs, moon)
    moon_sep = '{0:.2f}'.format(Angle(str(moon_seprad) + ' degrees').degree)

    # Calculation of Altitude, Airmass and Galactic Coordinates
    altitude = round(Angle(str(object_obs.alt) + ' degrees').degree, 2)
    airmass = round(1 / np.cos(np.radians(90 - altitude)), 2)
    l, b = get_galacticcoord(ra, dec)

    return moon_sep, moon_phase, altitude, airmass, l, b

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Manual Setup - GUI Code
# ------------------------------------------------------------------------------------------------------------------- #
# Choice of Input
bool_choice = eg.boolbox(msg='Enter The Type of Input', title='Input Data', choices=['DateList.dat', 'FITS Files'])

if bool_choice:
    datelist_df = pd.read_csv(list_dates, sep='\s+', comment='#')
else:
    ctext = eg.enterbox(msg='Enter The Common Text of FITS Files', title='Common Text', default='*.fits')

# Telescope Details
telescope = eg.enterbox(msg='Enter The Name of the Telescope', title='Name of the Telescope', default='DOT')
telescope_df = pd.read_csv(list_telescopes, sep='\s+', comment='#').set_index('ShortName')

if telescope in telescope_df.index.values:
    (OBS_NAME, OBS_LONG, OBS_LAT, OBS_ALT, OBS_TIMEZONE, _, _) = telescope_df.loc[telescope].values
else:
    telescope = eg.multenterbox('Enter the Necessary Details of the Telescope', title='Details of the Telescope',
                                fields=['Name', 'Longitude', 'Latitude', 'Altitude', 'TimeZone'],
                                values=telescope_df.loc['DOT'].values)
    (OBS_NAME, OBS_LONG, OBS_LAT, OBS_ALT, OBS_TIMEZONE) = telescope

# Are Input Times in UTC or Local Time?
utc = eg.boolbox(msg='Are Input Times in UTC or Local Time?', title='Input Times', choices=['UTC', 'Local Time'])
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Declaring Object 'telescope'
# ------------------------------------------------------------------------------------------------------------------- #
telescope = ephem.Observer()
telescope.pressure = 0
telescope.lon = OBS_LONG
telescope.lat = OBS_LAT
telescope.elevation = OBS_ALT
telescope.epoch = ephem.J2000
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculate the Moon Phase and Moon Angle for the Input (DateList.dat or FITS Files)
# ------------------------------------------------------------------------------------------------------------------- #
output_cols = ['FileName', 'Date', 'RA', 'DEC', 'MoonAngle', 'MoonPhase', 'Altitude', 'Airmass', 'GLON', 'GLAT']

if bool_choice:
    moon_df = datelist_df.copy()
    for index in datelist_df.index.values:
        time, ra, dec = tuple(datelist_df.loc[index])
        output_val = get_moonphaseandangle(time, ra, dec)
        modoutput_cols = output_cols[4:]
        for idx, column in enumerate(modoutput_cols):
            moon_df.loc[index, modoutput_cols[idx]] = output_val[idx]
else:
    moon_df = pd.DataFrame()
    for index, filename in enumerate(glob.glob(ctext)):
        time, ra, dec = get_timeradec(filename)
        output_val = get_moonphaseandangle(time, ra, dec)
        output_list = [filename, time, ra, dec] + list(output_val)
        print (output_list, output_cols)
        for idx, column in enumerate(output_cols):
            moon_df.loc[index, output_cols[idx]] = output_list[idx]

moon_df.to_csv('MoonPhaseAngle.asc', sep=' ', index=False)
# ------------------------------------------------------------------------------------------------------------------- #
