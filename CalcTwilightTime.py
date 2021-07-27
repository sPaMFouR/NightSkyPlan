#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxx----------CALCULATION OF TWILIGHT TIMES AND NIGHT DURATION---------xxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #


# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import sys
import ephem
import numpy as np
import pandas as pd
import easygui as eg
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import Angle
from datetime import datetime, timedelta

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.dates import DateFormatter, MonthLocator

from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

plt.style.use('bmh')
plt.rc('font', family='sans-serif')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables
# ------------------------------------------------------------------------------------------------------------------- #
list_telescopes = 'TelescopeList.dat'

dict_twilights = {'Civil': ['-6', True], 'Nautical': ['-12', True], 'Astronomical': ['-18', True],
                  'Sunset/Sunrise': ['-0.34', False], 'Moonset/Moonrise': ['-0.34', False]}
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Manual Setup - GUI Code
# ------------------------------------------------------------------------------------------------------------------- #
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

# Range of Dates for which Twilight Times have to be calculated
date_start = str(Time(Time.now(), format='iso', out_subfmt='date'))
date_end = str(Time(date_start, format='iso', out_subfmt='date') + 365 * u.d)
daterange = eg.multenterbox('Enter the Time Duration for which the Twilight Times need to be computed',
                            title='Time Duration', fields=['Date Start', 'Date End'], values=[date_start, date_end])

# Choose which columns to be printed to the Output File
columns = ['Sunset', 'DuskCivil', 'DuskNautical', 'DuskAstronomical', 'DawnAstronomical', 'DawnNautical', 'DawnCivil',
           'Sunrise', 'Moonrise', 'Moonset', 'NightDuration', 'NauticalNight', 'AstronomicalNight', 'MoonPhase']
printcols = eg.multchoicebox('Choose the Appropriate Columns for the Output file', title='Output Columns',
                             choices=columns, preselect=[2, 5, 11, 12, 13])

# Calculate Times in UTC or Local Time?
utc = eg.boolbox(msg='Calculate Times in UTC or Local Time?', title='Time Zone', choices=['UTC', 'Local Time'])
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
# Utility Functions
# ------------------------------------------------------------------------------------------------------------------- #

def sign(value):
    """
    Returns the sign of the input 'value'
    """
    return (float(value) > 0) - (float(value) < 0)


def display_text(text_to_display):
    """
    Displays text mentioned in the string 'text_to_display'
    Args:
        text_to_display : Text to be displayed
    Returns:
        None
    """
    print("\n" + "# " + "-" * (12 + len(text_to_display)) + " #")
    print("# " + "-" * 5 + " " + str(text_to_display) + " " + "-" * 5 + " #")
    print("# " + "-" * (12 + len(text_to_display)) + " #" + "\n")

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculate Sunset, Sunrise and Twilight Times
# ------------------------------------------------------------------------------------------------------------------- #

def calculate_twilighttime(category='Sunset/Sunrise'):
    """
    Computes Rising and Setting times for the twilight time specified by 'category'.
    Args:
        category : Specifies which Twilight times have to be coomputed
    Returns:
        setting  : Setting time for the twilight
        rising   : Rising time for the twilight
    """
    if category in dict_twilights:
        telescope.horizon = dict_twilights[category][0]
        setting = telescope.previous_setting(ephem.Sun(), use_center=dict_twilights[category][1])
        rising = telescope.next_rising(ephem.Sun(), use_center=dict_twilights[category][1])
        return setting, rising
    else:
        display_text("ERROR: Invalid Category Chosen '{0}'".format(category))
        sys.exit(1)


def calculate_moontime(time_midnight):
    """
    Computes Moonrise and Moonset at the time specified by 'time_midnight'.
    Args:
        time_midnight : Midnight time for the date on which Moonrise and Moonset is to be computed
    Returns:
        rising        : Time of Moonrise
        setting       : Time of Moonset
    """
    category = 'Moonset/Moonrise'
    telescope.horizon = dict_twilights[category][0]

    previousrise = telescope.previous_rising(ephem.Moon(), use_center=dict_twilights[category][1])
    previousset = telescope.previous_setting(ephem.Moon(), use_center=dict_twilights[category][1])
    nextset = telescope.next_setting(ephem.Moon(), use_center=dict_twilights[category][1])
    nextrise = telescope.next_rising(ephem.Moon(), use_center=dict_twilights[category][1])

    if previousrise.datetime() > previousset.datetime():
        rising = previousrise
        setting = nextset
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

    # Calculation of Maximum Phase of Moon During the Night
    list_moonphases = []
    for time in np.arange(time_prevnight, time_midnight, timedelta(minutes=15)).astype(datetime):
        time = datetime.strftime(time, '%Y/%m/%d %H:%M:%S')
        list_moonphases.append(ephem.Moon(time).phase)

    nightduration = str(sun_rise.datetime() - sun_set.datetime()).split('.')[0]
    nauticalnight = str(dawn_nauti.datetime() - dusk_nauti.datetime()).split('.')[0]
    astronomicalnight = str(dawn_astro.datetime() - dusk_astro.datetime()).split('.')[0]
    moonphase = '{0:.2f}'.format(np.max(list_moonphases))

    if utc:
        twilighttimes = [datetime.strftime(time.datetime(), '%H:%M:%S') for time in
                         [sun_set, dusk_civil, dusk_nauti, dusk_astro, dawn_astro,
                          dawn_nauti, dawn_civil, sun_rise, moon_rise, moon_set]]
    else:
        twilighttimes = [datetime.strftime(time.datetime() + timedelta(hours=OBS_TIMEZONE), '%H:%M:%S')
                         for time in [sun_set, dusk_civil, dusk_nauti, dusk_astro, dawn_astro,
                                      dawn_nauti, dawn_civil, sun_rise, moon_rise, moon_set]]

    dict_columns = {'Sunset': 0, 'DuskCivil': 1, 'DuskNautical': 2, 'DuskAstronomical': 3,
                    'DawnAstronomical': 4, 'DawnNautical': 5, 'DawnCivil': 6, 'Sunrise': 7, 'Moonrise': 8,
                    'Moonset': 9, 'NightDuration': nightduration, 'NauticalNight': nauticalnight,
                    'AstronomicalNight': astronomicalnight, 'MoonPhase': moonphase}

    for column, index in dict_columns.items():
        date_obs = str(Time(date_obs, format='iso', out_subfmt='date'))
        if type(index) == int:
            time_df.loc[date_obs, column] = twilighttimes[index]
        else:
            time_df.loc[date_obs, column] = index

    return time_df

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculate the Twilight Times for the Range of Dates
# ------------------------------------------------------------------------------------------------------------------- #
[datestart, dateend] = daterange
localdate_duration = np.arange(Time(datestart), Time(dateend) + 1 * u.d, 1 * u.d)

twilight_df = pd.DataFrame()
for localdate in localdate_duration:
    twilight_df = get_twilighttimes(twilight_df, localdate)

twilight_df.index.name = 'Date'
twilight_df[printcols].to_csv('TwilightTimes_{0}To{1}.asc'.format(datestart, dateend), sep='\t')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Night Duration as a Function of Date of Observation
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(18, 13))
ax = fig.add_subplot(111)

dates = [datetime.strptime(x, '%Y-%m-%d') for x in twilight_df.index.values]

nightduration = [(datetime.strptime(x, '%H:%M:%S') - datetime(1900, 1, 1)).total_seconds() / 3600
                 for x in twilight_df['NightDuration']]
nauticalnight = [(datetime.strptime(x, '%H:%M:%S') - datetime(1900, 1, 1)).total_seconds() / 3600
                 for x in twilight_df['NauticalNight']]
astronomicalnight = [(datetime.strptime(x, '%H:%M:%S') - datetime(1900, 1, 1)).total_seconds() / 3600
                     for x in twilight_df['AstronomicalNight']]

ax.plot(dates, nightduration, marker='o', ls='', mfc='dimgrey', mew=2, ms=7, c='orangered',
        alpha=0.5, label='Sunset-To-Sunrise')
ax.plot(dates, nauticalnight, marker='^', ls='', mfc='dimgrey', mew=2, ms=7, c='dodgerblue',
        alpha=0.5, label='Nautical Night')
ax.plot(dates, astronomicalnight, marker='s', ls='', mfc='dimgrey', mew=2, ms=7, c='navy',
        alpha=0.5, label='Astronomical Night')
# ax.fill_between(dates, nightduration, color='orangered', alpha=0.2)
# ax.fill_between(dates, nauticalnight, color='dodgerblue', alpha=0.2)
# ax.fill_between(dates, astronomicalnight, color='navy', alpha=0.2)

ax.legend(markerscale=2, frameon=True, fancybox=True, shadow=True, fontsize=14)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.yaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(0.1))
ax.xaxis.set_major_locator(MonthLocator())
ax.xaxis.set_major_formatter(DateFormatter('%Y-%m'))
ax.tick_params(axis='both', which='major', direction='in', width=1.6, length=9, labelsize=12)
ax.tick_params(axis='both', which='minor', direction='in', width=0.9, length=5, labelsize=12)

# Print Observatory Details
lat_deg = '%7.4f' % Angle(OBS_LAT + ' degrees').degree
long_deg = '%7.4f' % Angle(OBS_LONG + ' degrees').degree

text_ns = 'N'
text_ew = 'E'

if not sign(lat_deg):
    text_ns = 'S'
if not sign(long_deg):
    text_ew = 'W'

degree_sign = '$^\circ$'
text_name = OBS_NAME + ' [+' + str(OBS_TIMEZONE) + 'h]\n'
text_lat = 'Latitude : ' + lat_deg + degree_sign + text_ns
text_long = ', Longitude : ' + long_deg + degree_sign + text_ew
text_alt = ', Altitude : ' + str(OBS_ALT) + ' m'
display_text = text_name + text_lat + text_long + text_alt + '\n'

ax.set_xlabel('Date [YYYY-MM]', fontsize=16)
ax.set_ylabel('Nautical Night Duration [In Hours]', fontsize=16)
ax.set_title(display_text, fontsize=16)

fig.savefig('NightDuration_{0}To{1}.pdf'.format(datestart, dateend), format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #
