#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxx-------------------------NIGHT SKY PLANNER-------------------------xxxxxxxxxxxxxxxxxxxxxxxxx #
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
from matplotlib.ticker import FixedLocator
from matplotlib.dates import DateFormatter, MinuteLocator, HourLocator, DayLocator, MonthLocator

from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

plt.style.use('bmh')
plt.rc('font', family='sans-serif')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables
# ------------------------------------------------------------------------------------------------------------------- #
# List of Targets and Telescopes
list_targets = 'TargetList.dat'
list_telescopes = 'TelescopeList.dat'

# Defaults for computing Twilights and Setting/Rising Times
dict_twilights = {'Civil': ['-6', True], 'Nautical': ['-12', True], 'Astronomical': ['-18', True],
                  'Sunset/Sunrise': ['-0.34', False], 'Moonset/Moonrise': ['-0.34', False]}
                  
# Defaults Used In Plotting Trajectories (Including Colors)
time_offset = 0
object_count = 0
colors = ['r', 'sandybrown', 'gold', 'darkorange', 'salmon', 'deeppink', 'limegreen', 'teal', 'y', 'brown', 'c']
markers = ['o', '^', 'v', 'd', 'P', 'X', 'p', 'h', 'D', '4', '+', 's']
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Utility Functions
# ------------------------------------------------------------------------------------------------------------------- #

def sign(value):
    """
    Returns the sign of the input 'value'
    """
    return (float(value) > 0) - (float(value) < 0)


def remove_empty_values(list_values):
    """
    Args:
        list_values : Python list from which empty entries are to be removed
    Returns:
        list_values : Python list with empty entries removed
    """
    while True:
        try:
            list_values.remove('')
        except ValueError:
            break
    return list_values


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
# Class 'ObjectToObs' For Declaring Objects To Be Observed
# ------------------------------------------------------------------------------------------------------------------- #

class ObjectToObs:
    def __init__(self, name, ra, dec, plot_ax):
        self.name = name
        self.object = ephem.FixedBody()
        self.object._epoch = ephem.J2000
        self.object._ra = ra
        self.object._dec = dec
        self.ax = plot_ax
        self.list_alt = []

    def get_altitude(self, time_obs):
        global telescope
        telescope.date = str(time_obs)
        self.object.compute(telescope)
        altitude = Angle(str(self.object.alt) + ' degrees').degree
        return altitude

    def get_moonsep(self, time_obs):
        global telescope
        telescope.date = str(time_obs)
        self.object.compute(telescope)
        moon_pos = ephem.Moon(str(time_obs))
        angle_ephem = ephem.separation(self.object, moon_pos)
        angle_sep = int(Angle(str(angle_ephem) + ' degrees').degree)
        return angle_sep

    def plot_objtrack(self, utc=True):
        for utcdate in list(utcdate_duration):
            self.list_alt.append(self.get_altitude(str(utcdate.value)))
        if utc:
            self.plot_in_utc()
        else:
            self.plot_in_local()

    def plot_in_utc(self):
        global object_count
        self.ax.plot(utcdate_intervals, self.list_alt, c=colors[object_count],
                     marker=markers[object_count], ls='-', lw=1, ms=7, alpha=0.5, label=self.name)
        object_count += 1

    def plot_in_local(self):
        global object_count
        self.ax.plot(localdate_intervals, self.list_alt, color=colors[object_count],
                     ls='-', lw=1, ms=7, marker=markers[object_count], alpha=0.5, label=self.name)
        object_count += 1

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Trajectories of Target List as viewed from a Specified Observatory
# ------------------------------------------------------------------------------------------------------------------- #

def plot_obsplan(ax_obj, utc=True):
    """
    Sets plot parameters for plotting the trajectory of targets in sky.
    Args:
        ax_obj  : Axes object over which the observatory planning plot is displayed
        utc     : Boolean value to be determine whether UTC or Local Time is to be used for plotting
    Returns:
        None
    """
    global sunset, sunrise, duskcivil, dawncivil, dusknauti, dawnnauti, duskastro, dawnastro

    # Print Observatory Details
    # ------------------------------------------------------------------------------------------------------------- #
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
    # ------------------------------------------------------------------------------------------------------------- #

#     # Compute Time
#     # ------------------------------------------------------------------------------------------------------------- #
#     currenttime = Time.now()

#     if not utc:
#         sunset += OBS_TIMEZONE * u.hour
#         sunrise += OBS_TIMEZONE * u.hour
#         duskcivil += OBS_TIMEZONE * u.hour
#         dawncivil += OBS_TIMEZONE * u.hour
#         dusknauti += OBS_TIMEZONE * u.hour
#         dawnnauti += OBS_TIMEZONE * u.hour
#         duskastro += OBS_TIMEZONE * u.hour
#         dawnastro += OBS_TIMEZONE * u.hour
#         currenttime += OBS_TIMEZONE * u.hour
#     # ------------------------------------------------------------------------------------------------------------- #

#     # Print Text In The Plot
#     # ------------------------------------------------------------------------------------------------------------- #
#     ax_obj.text(sunset.value, 92, 'Sunset', rotation=50, c='orangered', fontsize=10)
#     ax_obj.text(sunrise.value - timedelta(minutes=10), 92, 'Sunrise', rotation=+50, c='orangered', fontsize=10)
#     ax_obj.text(duskcivil.value, 7, 'Civil Twilight', rotation=-90, c='navy', alpha=1, fontsize=10)
#     ax_obj.text(dawncivil.value, 7, 'Civil Twilight', rotation=-90, c='navy', alpha=1, fontsize=10)
#     ax_obj.text(dusknauti.value, 5, 'Nautical Twilight', rotation=-90, c='navy', alpha=1, fontsize=10)
#     ax_obj.text(dawnnauti.value, 5, 'Nautical Twilight', rotation=-90, c='navy', alpha=1, fontsize=10)
#     ax_obj.text(duskastro.value, 3, 'Astronomical Twilight', rotation=-90, c='navy', alpha=1, fontsize=10)
#     ax_obj.text(dawnastro.value, 3, 'Astronomical Twilight', rotation=-90, c='navy', alpha=1, fontsize=10)

#     nightspan = dawnnauti.value - dusknauti.value
#     midnight = dusknauti.value + nightspan / 2
#     printtime = datetime.strptime(str(currenttime).split('.')[0], '%Y-%m-%d %H:%M:%S')
#     printnightspan = 'Night Span = {0}'.format(str(nightspan))
#     printmoonphase = 'Moon Illumination = {0:.1f}%'.format(ephem.Moon(midnight).phase)

    ax_obj.text(datestart, HORIZON - 3, 'Telescope Horizon', c='navy', fontsize=10)
    ax_obj.text(datestart, ZENITH + 1, 'Telescope Zenith', c='navy', fontsize=10)

    ax_obj.set_facecolor('lightgrey')
    ax_obj.fill_between(ax_obj.get_xbound(), HORIZON - 0.5, HORIZON + 0.5, fc='dimgrey')
    ax_obj.fill_between(ax_obj.get_xbound(), ZENITH - 0.5, ZENITH + 0.5, fc='dimgrey')
    ax_obj.fill_between(ax_obj.get_xbound(), HORIZON + 0.5, ZENITH - 0.5, fc='k')

#     if sunset.value < currenttime.utc.datetime < sunrise.value:
#         ax_obj.axvline(x=currenttime.value, ls='--', color='k')
#         ax_obj.text(currenttime.value, 7, 'Current Time', rotation=-90, c='k', fontsize=10)
#     # ------------------------------------------------------------------------------------------------------------- #

#     # Fill Color In Sectors of Observation/Non-Observation
#     # ------------------------------------------------------------------------------------------------------------- #
#     ax_obj.axvline(x=sunset.value, ls='-', lw=1, color='k')
#     ax_obj.axvline(x=sunrise.value, ls='-', lw=1, color='k')
#     ax_obj.axvline(x=duskcivil.value, ls='-', lw=1, color='k')
#     ax_obj.axvline(x=dawncivil.value, ls='-', lw=1, color='k')
#     ax_obj.axvline(x=dusknauti.value, ls='-', lw=1, color='k')
#     ax_obj.axvline(x=dawnnauti.value, ls='-', lw=1, color='k')
#     ax_obj.axvline(x=duskastro.value, ls='-', lw=1, color='k')
#     ax_obj.axvline(x=dawnastro.value, ls='-', lw=1, color='k')

#     ax_obj.fill_between([sunset.value, duskcivil.value], HORIZON + 0.5, ZENITH - 0.5, fc='royalblue', alpha=1)
#     ax_obj.fill_between([dawncivil.value, sunrise.value], HORIZON + 0.5, ZENITH - 0.5, fc='royalblue', alpha=1)
#     ax_obj.fill_between([duskcivil.value, dusknauti.value], HORIZON + 0.5, ZENITH - 0.5, fc='royalblue', alpha=0.7)
#     ax_obj.fill_between([dawnnauti.value, dawncivil.value], HORIZON + 0.5, ZENITH - 0.5, fc='royalblue', alpha=0.7)
#     ax_obj.fill_between([dusknauti.value, duskastro.value], HORIZON + 0.5, ZENITH - 0.5, fc='royalblue', alpha=0.4)
#     ax_obj.fill_between([dawnastro.value, dawnnauti.value], HORIZON + 0.5, ZENITH - 0.5, fc='royalblue', alpha=0.4)
#     # ------------------------------------------------------------------------------------------------------------- #

    # Set Plotting Tick Parameters
    # ------------------------------------------------------------------------------------------------------------- #
    ax_obj.xaxis.set_ticks_position('both')
#     ax_obj.xaxis.set_major_locator(HourLocator())
    print (ax_obj.get_xticks())
    ax_obj.yaxis.set_major_locator(FixedLocator(range(0, 91, 10)))
    ax_obj.yaxis.set_minor_locator(FixedLocator(range(0, 91, 2)))
    ax_obj.xaxis.set_major_locator(MonthLocator(bymonthday=1))
    print (ax_obj.get_xticks())
    # ax_obj.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))
    print (ax_obj.get_xticks())
    ax_obj.tick_params(axis='both', which='major', direction='in', width=1.6, length=9, labelsize=12)
    ax_obj.tick_params(axis='both', which='minor', direction='in', width=0.9, length=5, labelsize=12)
    # ------------------------------------------------------------------------------------------------------------- #

    # Plot The Y-Axis On The RHS With Airmass
    # ------------------------------------------------------------------------------------------------------------- #
    list_secz = []
    for altitude in ax_obj.get_yticks():
        if (1 / np.cos(np.radians(90 - altitude))) < 10:
            list_secz.append('%5.2f' % (1 / np.cos(np.radians(90 - altitude))))
        else:
            list_secz.append('NaN')

    ax_twin = ax_obj.twinx()
    ax_twin.set_ylim(0, 90)
    ax_twin.set_yticklabels(list_secz)
    ax_twin.set_yticks(ax_obj.get_yticks())
    ax_twin.set_yticks(ax_obj.get_yticks(minor=True), minor=True)
    ax_twin.tick_params(axis='both', which='major', direction='in', width=1.6, length=9, labelsize=12)
    ax_twin.tick_params(axis='both', which='minor', direction='in', width=0.9, length=5, labelsize=12)
    # ------------------------------------------------------------------------------------------------------------- #

    # Set Plot Global Parameters
    # ------------------------------------------------------------------------------------------------------------- #
    ax_obj.set_ylim(0, 90)
    ax_obj.set_xlim(Time(datestart).value, Time(dateend).value)
    ax_obj.legend(title='Target List', loc='center left', bbox_to_anchor=(1.05, 0.5),
                  markerscale=1.6, ncol=1, shadow=True, fancybox=True, fontsize=14)
    ax_obj.grid(True, ls='--', lw=1)
    ax_obj.set_title(display_text, fontsize=16)
    ax_obj.set_ylabel('Elevation [In Degrees]', fontsize=16)
    ax_twin.set_ylabel('Airmass', fontsize=16)

    ax_obj.set_xlabel('\nDate Range: {0} To {1}'.format(datestart, dateend), fontsize=16)
    # ------------------------------------------------------------------------------------------------------------- #

    # Show and Save the Plot
    # ------------------------------------------------------------------------------------------------------------- #
    ax_obj.autoscale_view()
    fig.autofmt_xdate()
    fig.savefig('YearlyPlan_{0}To{1}.pdf'.format(datestart, dateend), format='pdf', dpi=2000, bbox_inches='tight')
    plt.show()
    plt.close(fig)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Manual Setup - GUI Code
# ------------------------------------------------------------------------------------------------------------------- #
# Target Details
target_df = pd.read_csv(list_targets, sep='\s+', comment='#').set_index('Index')
target_df = target_df[target_df['ToPlot'].isin(['y', 'Y'])]
field_names = ['Object {0}'.format(idx) for idx in target_df.index.values]
field_values = [target_df.loc[idx, 'Name'] + ' ' + target_df.loc[idx, 'RA'] + ' ' + target_df.loc[idx, 'DEC'] for idx
                in target_df.index.values]

# Telescope Details
telescope = eg.enterbox(msg='Enter The Name of the Telescope!', title='Name of the Telescope', default='HCT')
telescope_df = pd.read_csv(list_telescopes, sep='\s+', comment='#').set_index('ShortName')

if telescope in telescope_df.index.values:
    (OBS_NAME, OBS_LONG, OBS_LAT, OBS_ALT, OBS_TIMEZONE, HORIZON, ZENITH) = telescope_df.loc[telescope].values
else:
    print("ERROR: Observatory Name '{0}' not found in the file '{1}'".format(telescope, list_telescopes))

# List of Targets
box_msg = 'Verify Name, RA, DEC of objects for Observation planning'
box_title = 'Details of Objects'
list_values = eg.multenterbox(msg=box_msg, title=box_title, fields=field_names, values=field_values)
list_values = remove_empty_values(list_values)

while len(list_values) == 0:
    err_msg = box_msg + '\n\n ERROR: Aleast 1 Object required for Observation Planning!'
    list_values = eg.multenterbox(msg=err_msg, title=box_title, fields=field_names, values=list_values)
    list_values = remove_empty_values(list_values)

# Plot wrt UTC or Local Time?
utc = eg.boolbox(msg='Plot Trajectories w.r.t UTC or Local Time?', title='UTC Or Local Time?',
                 choices=['UTC', 'Local Time'])

# Range of Dates for which Twilight Times have to be calculated
date_start = str(Time(Time.now(), format='iso', out_subfmt='date'))
date_end = str(Time(date_start, format='iso', out_subfmt='date') + 365 * u.d)
daterange = eg.multenterbox('Enter the Time Duration for which the Twilight Times need to be computed',
                            title='Time Duration', fields=['Date Start', 'Date End'], values=[date_start, date_end])
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
# Determining Time Intervals
# Compute Twilight Times and Moonrise/set Times
# ------------------------------------------------------------------------------------------------------------------- #
[datestart, dateend] = daterange
utcdate_duration = np.arange(Time(datestart), Time(dateend) + 1 * u.d, 1 * u.d)
utcdate_intervals = [utcdate.value for utcdate in utcdate_duration]
localdate_intervals = [(utcdate + OBS_TIMEZONE * u.hour).value for utcdate in utcdate_duration]

twilight_df = pd.DataFrame()
for utcdate in utcdate_duration:
    twilight_df = get_twilighttimes(twilight_df, utcdate)
twilight_df.index.name = 'Date'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The Trajectories Of Objects To Be Observed
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(18, 13))
ax = fig.add_subplot(111)

for index, value in enumerate(list_values):
    if len(value.split()) >= 3:
        ObjectToObs(name=value.split()[-3], ra=value.split()[-2], dec=value.split()[-1],
                    plot_ax=ax).plot_objtrack(utc=utc)
    elif len(value.split()) == 2:
        ObjectToObs(name='Object ' + str(int(index) + 1), ra=value.split()[-2], dec=value.split()[-1],
                    plot_ax=ax).plot_objtrack(utc=utc)
    else:
        print("ERROR : Both RA & DEC for the Object '{}' need to be specified".format(str(int(index) + 1)))
        continue

plot_obsplan(ax_obj=ax, utc=utc)
# ------------------------------------------------------------------------------------------------------------------- #
