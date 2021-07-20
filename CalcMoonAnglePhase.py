#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxx----------CALCULATION OF TWILIGHT TIMES AND NIGHT DURATION---------xxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #


# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import glob
import ephem
import pandas as pd
import easygui as eg
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import Angle
from datetime import datetime, timedelta
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
# Manual Setup - GUI Code
# ------------------------------------------------------------------------------------------------------------------- #
# Choice of Input
bool_choice = eg.boolbox(msg='Enter The Type of Input', title='Input Data', choices=['Date List', 'FITS Files'])

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

    time_obs = header[TIME_keyword]
    ra = header[RA_keyword]
    dec = header[DEC_keyword]
    
    return time_obs, ra, dec


def get_moonphaseandangle(time_obs, ra, dec):
    """
    The function computes moon separation for the object ('ra', 'dec') at the time of observation 
    and the moon phase at the same time.
    Args:
        time_obs   : Time of observation for which moon angle and moon separation have to be computed
        ra         : Right Ascension of the source for which moon angle has to be computed
        dec        : Declination of the source for which moon angle has to be computed
    Returns:
        moon_sep   : Computed Moon separation for the object ('ra', 'dec') at the time of observation
        moon_phase : Computed Moon Phase at the time of observation
    """
    # Declare the Object of Observation
    object = ephem.FixedBody()
    object._epoch = ephem.J2000
    object._ra = ra
    object._dec = dec

    time_obs = datetime.strptime(time_obs, datetime_format)
    if not utc:
        time_obs += timedelta(hours=OBS_TIMEZONE)

    # Calculation of Moon Phase and Moon Angle
    telescope.date = time_obs
    object.compute(telescope)
    moon_phase = '{0:.2f}'.format(ephem.Moon(time_obs).phase)

    moon_pos = ephem.Moon(str(time_obs))
    moon_seprad = ephem.separation(object, moon_pos)
    moon_sep = '{0:.2f}'.format(Angle(str(moon_seprad) + ' degrees').degree)

    return moon_sep, moon_phase

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculate the Moon Phase and Moon Angle for the Input (DateList.dat or FITS Files)
# ------------------------------------------------------------------------------------------------------------------- #
if bool_choice:
    moon_df = datelist_df.copy()
    for index in datelist_df.index.values:
        time_obs, ra, dec = tuple(datelist_df.loc[index])
        moon_sep, moon_phase = get_moonphaseandangle(time_obs, ra, dec)
        moon_df.loc[index, 'MoonAngle'] = moon_sep
        moon_df.loc[index, 'MoonPhase'] = moon_phase
else:
    moon_df = pd.DataFrame()
    for index, filename in enumerate(glob.glob(ctext)):
        time_obs, ra, dec = get_timeradec(filename)
        moon_sep, moon_phase = get_moonphaseandangle(time_obs, ra, dec)
        moon_df.loc[index, 'FileName'] = filename
        moon_df.loc[index, 'Date'] = time_obs
        moon_df.loc[index, 'RA'] = ra
        moon_df.loc[index, 'DEC'] = dec
        moon_df.loc[index, 'MoonAngle'] = moon_sep
        moon_df.loc[index, 'MoonPhase'] = moon_phase
        print (moon_df)
        
moon_df.to_csv('MoonPhaseAngle.asc', sep=' ', index=False)
# ------------------------------------------------------------------------------------------------------------------- #