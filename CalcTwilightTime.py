#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxx-------------------------NIGHT SKY PLANNER-------------------------xxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #


# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import ephem
import numpy as np
import pandas as pd
import easygui as eg
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import Angle
from datetime import datetime, timedelta

from matplotlib import rc
import matplotlib.colors as mcolors
from matplotlib import pyplot as plt
from matplotlib.ticker import FixedLocator
from matplotlib.dates import DateFormatter, MinuteLocator, HourLocator

from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

plt.style.use('bmh')
plt.rc('font', family='sans-serif')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables
# ------------------------------------------------------------------------------------------------------------------- #
list_telescopes = 'TelescopeList.dat'

dict_twilights = {'Civil' : ['-6', True], 'Nautical': ['-12', True], 'Astronomical': ['-18', True],
                  'Sunset/Sunrise': ['-0.34', False], 'Moonset/Moonrise': ['-0.34', False]}
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Manual Setup - GUI Code
# ------------------------------------------------------------------------------------------------------------------- #
# Telescope Details
telescope = eg.enterbox(msg='Enter The Name of the Telescope', title='Name of the Telescope', default='DOT')
telescope_df = pd.read_csv(list_telescopes, sep='\s+', comment='#').set_index('ShortName')

if telescope in telescope_df.index.values:
    (OBS_NAME, OBS_LONG, OBS_LAT, OBS_ALT, OBS_TIMEZONE) = telescope_df.loc[telescope].values
else:
    telescope = eg.multenterbox('Enter the Site Details of the Telescope [Default = DOT]',
                                title='Details of the Telescope',
                                fields=['Name', 'Longitude', 'Latitude', 'Altitude', 'TimeZone'],
                                values=telescope_df.loc['DOT'].values)
    (OBS_NAME, OBS_LONG, OBS_LAT, OBS_ALT, OBS_TIMEZONE) = telescope

# Range of Dates for which Twilight Times have to be calculated
date_start = str(Time(Time.now(), format='iso', out_subfmt='date'))
date_end = str(Time(date_start, format='iso', out_subfmt='date') + 120 * u.d)
daterange = eg.multenterbox('Enter the Time Duration for which the Twilight Times need to be computed',
                            title='Time Duration', fields=['Date Start', 'Date End'], values=[date_start, date_end])

# Choose which columns to be printed to the output file
printcols = eg.multchoicebox('Choose the Appropriate Columns for the Output file', title='Output Columns',
                             choices=['Sunset', 'DuskCivil', 'DuskNautical', 'DuskAstronomical',
                                      'DawnAstronomical', 'DawnNautical', 'DawnCivil', 'Sunrise',
                                      'NightDuration', 'MoonIllumination', 'Moonrise', 'Moonset'],
                             preselect=[2, 5, 8, 9])

# Plot Times in UTC or Local Time?
utc = eg.boolbox(msg='Plot Times in UTC or Local Time?', title='UTC Or Local Time?', choices=['UTC', 'Local Time'])
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
# Calculate Sunset, Sunrise and Twilight Times
# ------------------------------------------------------------------------------------------------------------------- #

def calculate_twilighttime(category='Sunset/Sunrise'):
    """
    It will compute Rising and Setting times for the twilight time specified by 'twilight'.
    Args:
        category : Specifies which Twilight times have to be coomputed
    Returns:
        setting  : Setting time for the twilight
        rising   : Rising time for the twilight
    """
    telescope.horizon = dict_twilights[category][0]
    setting = telescope.previous_setting(ephem.Sun(), use_center=dict_twilights[category][1])
    rising = telescope.next_rising(ephem.Sun(), use_center=dict_twilights[category][1])
    
    return setting, rising


def calculate_moontime(time_midnight):
    """
    It will compute Moonrise and Moonset times.
    Args:
        time_midnight : Midnight time for the date on which Moonrise and Moonset is to be computed
    Returns:
        rising        : Time of Moonrise
        setting       : Time of Moonset
    """
    category='Moonset/Moonrise'
    telescope.horizon = dict_twilights[category][0]
    
    previousrise = telescope.previous_rising(ephem.Moon(), use_center=dict_twilights[category][1])
    previousset = telescope.previous_setting(ephem.Moon(), use_center=dict_twilights[category][1])
    nextset = telescope.next_setting(ephem.Moon(), use_center=dict_twilights[category][1])
    nextrise = telescope.next_rising(ephem.Moon(), use_center=dict_twilights[category][1])

    if previousrise.datetime() > previousset.datetime():
        rising = previousrise
        setting = nextsset
    if time_midnight - previousset.datetime() > nextrise.datetime() - time_midnight:
        rising = nextrise
        setting = nextset
    else:
        rising = previousrise
        setting = previousset

    return rising, setting


def get_twilighttimes(time_df, date_obs):
    """
    Obtain Sunset, Sunrise and Twilight times for the observatory on a given date 'date_obs'.
    Args:
        time_df  : Input Pandas DataFrame to which twilight times have to be appended
        date_obs : Date of observation for which twilight times have to be computed
    Returns:
        time_df  : Output Pandas DataFrame to which twilight times have been appended
    """
    time_prevnight = (Time(date_obs) - abs(OBS_TIMEZONE) * u.hour).datetime 
    time_midnight = (Time(date_obs) + 1 * u.day - abs(OBS_TIMEZONE) * u.hour).datetime
    telescope.date = time_midnight

    # Calculation Of Local MoonRise & MoonSet [Refraction Correction Of Moon = -0.34 Degrees]
    moon_rise, moon_set = calculate_moontime(time_midnight)

    # Calculation Of Local Sunset & Sunrise [Refractrion Correction for Sun = -0.34 Degrees]
    sun_set, sun_rise = calculate_twilighttime('Sunset/Sunrise')
    
    # Calculation Of Civil Twilight [Elevation Of Sun = -6 Degrees]
    dusk_civil, dawn_civil = calculate_twilighttime('Civil')

    # Calculation Of Nautical Twilight [Elevation Of Sun = -12 Degrees]
    dusk_nauti, dawn_nauti = calculate_twilighttime('Nautical')
        
    # Calculation Of Astronomical Twilight [Elevation Of Sun = -18 Degrees]
    dusk_astro, dawn_astro = calculate_twilighttime('Astronomical')
    
    if utc:
        twilighttimes = [datetime.strftime(time.datetime(), '%H:%M:%S') for time in
                         [sun_set, dusk_civil, dusk_nauti, dusk_astro, dawn_astro,
                          dawn_nauti, dawn_civil, sun_rise, moon_rise, moon_set]]
    else:
        twilighttimes = [datetime.strftime(time.datetime() + timedelta(hours=OBS_TIMEZONE), '%H:%M:%S')
                         for time in [sun_set, dusk_civil, dusk_nauti, dusk_astro, dawn_astro,
                                      dawn_nauti, dawn_civil, sun_rise, moon_rise, moon_set]]

    (sunset, duskcivil, dusknauti, duskastro, dawnastro,
        dawnnauti, dawncivil, sunrise, moonrise, moonset) = twilighttimes

    # Calculation of Maximum Illumination of Moon During the Night
    list_illumination = []
    for time in np.arange(time_prevnight, time_midnight, timedelta(minutes=15)).astype(datetime):
        time = datetime.strftime(time, '%Y/%m/%d %H:%M:%S')
        list_illumination.append(ephem.Moon(time).phase)

    nightduration = str(dawn_astro.datetime() - dusk_astro.datetime()).split('.')[0]
    illumination = "{0:.2f}".format(np.max(list_illumination))

    dict_columns = {'Sunset': sunset, 'DuskCivil': duskcivil, 'DuskNautical': dusknauti,
                    'DuskAstronomical': duskastro, 'DawnAstronomical': dawnastro, 'DawnNautical': dawnnauti,
                    'DawnCivil': dawncivil, 'Sunrise': sunrise, 'NightDuration': nightduration,
                    'MoonIllumination': illumination, 'Moonrise': moonrise, 'Moonset': moonset}

    for column, time in dict_columns.items():
        if column in printcols:
            date_obs = str(Time(date_obs, format='iso', out_subfmt='date'))
            time_df.loc[date_obs, column] = time

    return time_df

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculate the Twilight Times for the Range of Dates
# ------------------------------------------------------------------------------------------------------------------- #
[datestart, dateend] = daterange
localdate_duration = np.arange(Time(datestart), Time(dateend), 1 * u.d)

twilight_df = pd.DataFrame()
for localdate in localdate_duration:
    twilight_df = get_twilighttimes(twilight_df, localdate)
twilight_df

twilight_df.index.name = 'Date'
twilight_df.to_csv('TwilightTimes_{0}To{1}'.format(datestart, dateend), sep='\t')
# ------------------------------------------------------------------------------------------------------------------- #
