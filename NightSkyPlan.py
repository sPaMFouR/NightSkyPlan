#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxx-------------------------NIGHT SKY PLANNER-------------------------xxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #


# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import os
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
from matplotlib.dates import DateFormatter, MinuteLocator, HourLocator

from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

plt.style.use('bmh')
plt.rc('font', family='sans-serif')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables
# Target List
# ------------------------------------------------------------------------------------------------------------------- #
list_targets = 'TargetList.dat'
list_telescopes = 'TelescopeList.dat'

dict_twilights = {'Civil': ['-6', True], 'Nautical': ['-12', True], 'Astronomical': ['-18', True],
                  'Sunset/Sunrise': ['-0.34', False], 'Moonset/Moonrise': ['-0.34', False]}
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Defaults Used In Plotting Trajectories (Including Colors)
# ------------------------------------------------------------------------------------------------------------------- #
time_offset = 0
object_count = 0

colors = ['r', 'sandybrown', 'gold', 'darkorange', 'salmon', 'deeppink', 'limegreen', 'teal', 'y', 'brown', 'c']
markers = ['o', '^', 'v', 'd', 'P', 'X', 'p', 'h', 'D', '4', '+', 's']
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function For Removing Empty Entries
# ------------------------------------------------------------------------------------------------------------------- #

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

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Manual Setup - GUI Code
# ------------------------------------------------------------------------------------------------------------------- #
# Telescope Details
telescope = eg.enterbox(msg='Enter The Name of the Telescope!', title='Name of the Telescope', default='HCT')
telescope_df = pd.read_csv(list_telescopes, sep='\s+', comment='#').set_index('ShortName')

if telescope in telescope_df.index.values:
    (OBS_NAME, OBS_LONG, OBS_LAT, OBS_ALT, OBS_TIMEZONE, HORIZON, ZENITH) = telescope_df.loc[telescope].values
else:
    print("ERROR: Observatory Name '{0}' not found in the file '{1}'".format(telescope, list_telescopes))

# List of Targets
if os.path.exists(list_targets):
    target_df = pd.read_csv(list_targets, sep='\s+', comment='#').set_index('Index')
    target_df = target_df[target_df['ToPlot'].isin(['y', 'Y'])]
    field_names = ['Object {0}'.format(idx) for idx in target_df.index.values]
    field_values = [target_df.loc[idx, 'Name'] + ' ' + target_df.loc[idx, 'RA'] + ' ' + target_df.loc[idx, 'DEC']
                    for idx in target_df.index.values]
else:
    field_names = ['Object {0}'.format(idx + 1) for idx in range(8)]
    field_values = [''] * 8

box_msg = 'Verify Name, RA, DEC of objects for Observation planning'
box_title = 'Details of Objects'
list_values = eg.multenterbox(msg=box_msg, title=box_title, fields=field_names, values=field_values)
list_values = remove_empty_values(list_values)

while len(list_values) == 0:
    err_msg = box_msg + '\n\n ERROR: Aleast 1 Object required for Observation Planning!'
    list_values = eg.multenterbox(msg=err_msg, title=box_title, fields=field_names, values=list_values)
    list_values = remove_empty_values(list_values)

# Plot Trajectories in UTC or Local Time?
utc = eg.boolbox(msg='Plot Trajectories in UTC or Local Time?', title='Time Zone', choices=['UTC', 'Local Time'])

# Current Date or Manually Entered Date?
date_obs = str(Time(Time.now(), format='iso', out_subfmt='date'))
setup_manual = eg.boolbox(msg='Manually Enter Date?', title='Manual or Current Date?', choices=['Manual', 'Current'])

if setup_manual:
    date_obs = eg.enterbox(msg='Enter The Date Of Observation!', title='Date Of Observation', default=date_obs)
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

time_midnight = (Time(date_obs) + 1 * u.day - abs(OBS_TIMEZONE) * u.hour).datetime
telescope.date = time_midnight
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

twilighttimes = [Time(datetime.strptime(str(time).split('.')[0], '%Y/%m/%d %H:%M:%S'))
                 for time in [sun_set, dusk_civil, dusk_nauti, dusk_astro, dawn_astro, 
                              dawn_nauti, dawn_civil, sun_rise, moon_rise, moon_set]]

(sunset, duskcivil, dusknauti, duskastro, dawnastro, dawnnauti,
                               dawncivil, sunrise, moonrise, moonset) = twilighttimes
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Determining Time Intervals
# ------------------------------------------------------------------------------------------------------------------- #
plot_duration = (sunrise.utc.datetime - sunset.utc.datetime).total_seconds() / 3600.
utctime_intervals = sunset + np.linspace(time_offset, time_offset + plot_duration, 140) * u.hour
localtime_intervals = utctime_intervals + OBS_TIMEZONE * u.hour
moonsep_intervals = sunset + np.linspace(time_offset, time_offset + plot_duration, 7)[1:-1] * u.hour
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Class 'ObjectToObs' For Declaring Objects To Be Observed
# ------------------------------------------------------------------------------------------------------------------- #

class ObjectToObs:
    def __init__(self, object_name, object_ra, object_dec, plot_ax):
        self.name = object_name
        self.object = ephem.FixedBody()
        self.object._epoch = ephem.J2000
        self.object._ra = object_ra
        self.object._dec = object_dec
        self.ax = plot_ax
        self.list_alt = []

    def get_altitude(self, time_obs):
        global telescope
        telescope.date = str(time_obs)
        self.object.compute(telescope)
        object_alt = Angle(str(self.object.alt) + ' degrees').degree
        return object_alt

    def get_moonsep(self, time_obs):
        global telescope
        telescope.date = str(time_obs)
        self.object.compute(telescope)
        moon_pos = ephem.Moon(str(time_obs))
        angle_ephem = ephem.separation(self.object, moon_pos)
        angle_sep = int(Angle(str(angle_ephem) + ' degrees').degree)
        return angle_sep

    def plot_objtrack(self, utc=True):
        for time_obs in list(utctime_intervals.value):
            self.list_alt.append(self.get_altitude(str(time_obs)))
        if utc:
            self.plot_in_utc()
        else:
            self.plot_in_local()

    def plot_in_utc(self):
        global object_count
        self.ax.plot(list(utctime_intervals.value), self.list_alt, c=colors[object_count],
                     marker=markers[object_count], ls='-', lw=1, ms=7, alpha=0.5,
                     label='{0} [{1:}$^\circ$]'.format(self.name, self.get_moonsep(str(moonsep_intervals[0]))))

        # plot_intervals = [time for time in moonsep_intervals if int(self.get_altitude(str(time))) > 0]
        # for time_obs in plot_intervals:
        #     self.ax.text(time_obs.value, self.get_altitude(str(time_obs)) + 0.5, self.get_moonsep(str(time_obs)),
        #                  fontsize=9, color='white', alpha=0.8)
        object_count += 1

    def plot_in_local(self):
        global object_count
        self.ax.plot(list(localtime_intervals.value), self.list_alt, color=colors[object_count], ls='-', lw=1, ms=7,
                     marker=markers[object_count], alpha=0.5, label='{0} [{1:}$^\circ$]'.format(self.name,
                     self.get_moonsep(str(moonsep_intervals[0]))))

        # plot_intervals = [time for time in moonsep_intervals if int(self.get_altitude(str(time))) > 0]
        # for time_obs in plot_intervals:
        #     local_time = time_obs + OBS_TIMEZONE * u.hour
        #     self.ax.text(local_time.value, self.get_altitude(str(time_obs)) + 0.5, self.get_moonsep(str(time_obs)),
        #                  fontsize=9, color='white', alpha=0.8)
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

    def sign(value):
        return (float(value) > 0) - (float(value) < 0)

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

    # Compute Time
    # ------------------------------------------------------------------------------------------------------------- #
    currenttime = Time.now()

    if not utc:
        sunset += OBS_TIMEZONE * u.hour
        sunrise += OBS_TIMEZONE * u.hour
        duskcivil += OBS_TIMEZONE * u.hour
        dawncivil += OBS_TIMEZONE * u.hour
        dusknauti += OBS_TIMEZONE * u.hour
        dawnnauti += OBS_TIMEZONE * u.hour
        duskastro += OBS_TIMEZONE * u.hour
        dawnastro += OBS_TIMEZONE * u.hour
        currenttime += OBS_TIMEZONE * u.hour
    # ------------------------------------------------------------------------------------------------------------- #

    # Print Text In The Plot
    # ------------------------------------------------------------------------------------------------------------- #
    ax_obj.text(sunset.value, 91, 'Sunset', rotation=50, c='navy', fontsize=9)
    ax_obj.text(sunrise.value - timedelta(minutes=10), 91, 'Sunrise', rotation=+50, c='navy', fontsize=9)
    ax_obj.text(duskcivil.value, 3, 'Civil', rotation=-90, c='navy', alpha=1, fontsize=9)
    ax_obj.text(dawncivil.value, 3, 'Civil', rotation=-90, c='navy', alpha=1, fontsize=9)
    ax_obj.text(dusknauti.value, 3, 'Nautical', rotation=-90, c='navy', alpha=1, fontsize=9)
    ax_obj.text(dawnnauti.value, 3, 'Nautical', rotation=-90, c='navy', alpha=1, fontsize=9)
    ax_obj.text(duskastro.value, 3, 'Astronomical', rotation=-90, c='navy', alpha=1, fontsize=9)
    ax_obj.text(dawnastro.value, 3, 'Astronomical', rotation=-90, c='navy', alpha=1, fontsize=9)

    nightspan = dawnnauti.value - dusknauti.value
    midnight = dusknauti.value + nightspan / 2
    printtime = datetime.strptime(str(currenttime).split('.')[0], '%Y-%m-%d %H:%M:%S')
    printnightspan = 'Night Span = {0}'.format(str(nightspan))
    printmoonphase = 'Moon Illumination = {0:.1f}%'.format(ephem.Moon(midnight).phase)
    printdate = 'Date of Observation : {0}'.format(date_obs)

    ax_obj.text(midnight - timedelta(minutes=25), HORIZON - 2, 'Telescope Horizon', c='navy', fontsize=9)
    ax_obj.text(midnight - timedelta(minutes=25), ZENITH + 1, 'Telescope Zenith', c='navy', fontsize=9)
    ax_obj.text(midnight - timedelta(hours=2), 91, s=printnightspan + ', ' + printmoonphase, c='r', fontsize=12)

    if sunset.value < currenttime.utc.datetime < sunrise.value:
        ax_obj.axvline(x=currenttime.value, ls='--', lw=1, color='red')
        ax_obj.text(currenttime.value, 3, 'Current Time', rotation=-90, c='red', fontsize=10)
    # ------------------------------------------------------------------------------------------------------------- #

    # Fill Color In Sectors of Observation/Non-Observation
    # ------------------------------------------------------------------------------------------------------------- #
    ax_obj.axvline(x=sunset.value, ls='-', lw=1, color='k')
    ax_obj.axvline(x=sunrise.value, ls='-', lw=1, color='k')
    ax_obj.axvline(x=duskcivil.value, ls='-', lw=1, color='k')
    ax_obj.axvline(x=dawncivil.value, ls='-', lw=1, color='k')
    ax_obj.axvline(x=dusknauti.value, ls='-', lw=1, color='k')
    ax_obj.axvline(x=dawnnauti.value, ls='-', lw=1, color='k')
    ax_obj.axvline(x=duskastro.value, ls='-', lw=1, color='k')
    ax_obj.axvline(x=dawnastro.value, ls='-', lw=1, color='k')

    ax_obj.set_facecolor('lightgrey')
    ax_obj.fill_between(ax_obj.get_xbound(), HORIZON - 0.5, HORIZON + 0.5, fc='dimgrey')
    ax_obj.fill_between(ax_obj.get_xbound(), ZENITH - 0.5, ZENITH + 0.5, fc='dimgrey')
    ax_obj.fill_between(ax_obj.get_xbound(), HORIZON + 0.5, ZENITH - 0.5, fc='k')

    ax_obj.fill_between([sunset.value, duskcivil.value], HORIZON + 0.5, ZENITH - 0.5, fc='royalblue', alpha=1)
    ax_obj.fill_between([dawncivil.value, sunrise.value], HORIZON + 0.5, ZENITH - 0.5, fc='royalblue', alpha=1)
    ax_obj.fill_between([duskcivil.value, dusknauti.value], HORIZON + 0.5, ZENITH - 0.5, fc='royalblue', alpha=0.7)
    ax_obj.fill_between([dawnnauti.value, dawncivil.value], HORIZON + 0.5, ZENITH - 0.5, fc='royalblue', alpha=0.7)
    ax_obj.fill_between([dusknauti.value, duskastro.value], HORIZON + 0.5, ZENITH - 0.5, fc='royalblue', alpha=0.4)
    ax_obj.fill_between([dawnastro.value, dawnnauti.value], HORIZON + 0.5, ZENITH - 0.5, fc='royalblue', alpha=0.4)
    # ------------------------------------------------------------------------------------------------------------- #

    # Set Plotting Tick Parameters
    # ------------------------------------------------------------------------------------------------------------- #
    ax_obj.xaxis.set_ticks_position('both')
    ax_obj.xaxis.set_major_locator(HourLocator())
    ax_obj.yaxis.set_major_locator(FixedLocator(range(0, 91, 10)))
    ax_obj.yaxis.set_minor_locator(FixedLocator(range(0, 91, 2)))
    ax_obj.xaxis.set_minor_locator(MinuteLocator(byminute=range(0, 60, 10)))
    ax_obj.xaxis.set_major_formatter(DateFormatter('%H:%M'))
    ax_obj.tick_params(axis='both', which='major', direction='in', width=1.6, length=8, labelsize=12)
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
    ax_twin.set_yticks(ax_obj.get_yticks())
    ax_twin.set_yticks(ax_obj.get_yticks(minor=True), minor=True)
    ax_twin.set_yticklabels(list_secz)
    ax_twin.tick_params(axis='both', which='major', direction='in', width=1.6, length=8, labelsize=12)
    ax_twin.tick_params(axis='both', which='minor', direction='in', width=0.9, length=5, labelsize=12)
    # ------------------------------------------------------------------------------------------------------------- #

    # Set Plot Global Parameters
    # ------------------------------------------------------------------------------------------------------------- #
    ax_obj.set_ylim(0, 90)
    ax_obj.set_xlim(sunset.value, sunrise.value)
    ax_obj.legend(title='Target List [Moon Angle]', loc='center left', bbox_to_anchor=(1.05, 0.5),
                  markerscale=1.6, ncol=1, shadow=True, fancybox=True, fontsize=14)
    ax_obj.grid(True, ls='--', lw=1)
    ax_obj.set_title(display_text, fontsize=16)
    ax_obj.set_ylabel('Elevation [In Degrees]', fontsize=16)
    ax_twin.set_ylabel('Airmass', fontsize=16)

    if not utc:
        ax_obj.set_xlabel('\nLocal Time [In Hours]\n' + str(printdate), fontsize=16)
    else:
        ax_obj.set_xlabel('\nUniversal Time [In Hours]\n' + str(printdate), fontsize=16)
    # ------------------------------------------------------------------------------------------------------------- #

    # Show and Save the Plot
    # ------------------------------------------------------------------------------------------------------------- #
    ax_obj.autoscale_view()
    fig.autofmt_xdate()
    fig.savefig('NightSkyPlan_{0}.pdf'.format(date_obs), format='pdf', dpi=2000, bbox_inches='tight')
    plt.show()
    plt.close(fig)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The Trajectories Of Objects To Be Observed
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(18, 13))
ax = fig.add_subplot(111)

for index, value in enumerate(list_values):
    if len(value.split()) >= 3:
        ObjectToObs(object_name=value.split()[-3], object_ra=value.split()[-2],
                    object_dec=value.split()[-1], plot_ax=ax).plot_objtrack(utc=utc)
    elif len(value.split()) == 2:
        ObjectToObs(object_name='Object ' + str(int(index) + 1), object_ra=value.split()[-2],
                    object_dec=value.split()[-1], plot_ax=ax).plot_objtrack(utc=utc)
    else:
        print("ERROR : Both RA & DEC for the Object '{}' need to be specified".format(str(int(index) + 1)))
        continue

plot_obsplan(ax_obj=ax, utc=utc)
# ------------------------------------------------------------------------------------------------------------------- #
