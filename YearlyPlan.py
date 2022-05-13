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
from datetime import datetime
from astropy.coordinates import Angle

from matplotlib import pyplot as plt
from matplotlib.ticker import FixedLocator
from matplotlib.dates import DateFormatter, MonthLocator

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

    def plot_objtrack(self, datelist):
        global object_count
        for date in list(datelist):
            self.list_alt.append(self.get_altitude(str(date)))

        self.ax.plot(datelist, self.list_alt, c=colors[object_count % len(colors)],
                     marker=markers[object_count % len(markers)], ls='-', lw=1, ms=8, alpha=0.5, label=self.name)
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

    # Fill Color In Sectors of Observation/Non-Observation
    # ------------------------------------------------------------------------------------------------------------- #

    ax_obj.text(datemid, HORIZON - 3, 'Telescope Horizon', c='navy', fontsize=10)
    ax_obj.text(datemid, ZENITH + 1, 'Telescope Zenith', c='navy', fontsize=10)

    ax_obj.set_facecolor('lightgrey')
    ax_obj.fill_between(ax_obj.get_xbound(), HORIZON - 0.5, HORIZON + 0.5, fc='dimgrey')
    ax_obj.fill_between(ax_obj.get_xbound(), ZENITH - 0.5, ZENITH + 0.5, fc='dimgrey')
    ax_obj.fill_between(ax_obj.get_xbound(), HORIZON + 0.5, ZENITH - 0.5, fc='k')
    # ------------------------------------------------------------------------------------------------------------- #

    # Set Plot Tick Parameters
    # ------------------------------------------------------------------------------------------------------------- #
    ax_obj.xaxis.set_ticks_position('both')
    ax_obj.yaxis.set_major_locator(FixedLocator(range(0, 91, 10)))
    ax_obj.yaxis.set_minor_locator(FixedLocator(range(0, 91, 2)))
    ax_obj.xaxis.set_major_locator(MonthLocator(bymonthday=1))
    ax_obj.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))

    ax_obj.tick_params(axis='both', which='major', direction='in', width=1.6, length=9, labelsize=12)
    ax_obj.tick_params(axis='both', which='minor', direction='in', width=0.9, length=5, labelsize=12)
    # ------------------------------------------------------------------------------------------------------------- #

    # Obtain the AIRMASS and Plot it on RHS of Y-Axis
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
    ax_obj.grid(True, ls='--', lw=1)
    ax_obj.legend(title='Target List', loc='center left', bbox_to_anchor=(1.05, 0.5),
                  markerscale=1.6, ncol=1, shadow=True, fancybox=True, fontsize=14)

    ax_obj.set_ylim(0, 90)
    ax_obj.set_xlim([datetime.fromisoformat(datestart), datetime.fromisoformat(dateend)])              
    ax_obj.set_ylabel('Elevation [In Degrees]', fontsize=16)
    ax_obj.set_title(display_text, fontsize=16)

    ax_twin.set_ylabel('Airmass', fontsize=16)
    ax_obj.set_xlabel('\nDate Range: {0} - {1}'.format(datestart, dateend), fontsize=16)
    # ------------------------------------------------------------------------------------------------------------- #

    # Show and Save the Plot
    # ------------------------------------------------------------------------------------------------------------- #
    ax_obj.autoscale_view()
    fig.autofmt_xdate()
    fig.savefig('YearlyPlan_{0}To{1}.pdf'.format(datestart, dateend), format='pdf', dpi=2000, bbox_inches='tight')
    # plt.show()
    plt.close(fig)

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
    target_df = pd.read_csv(list_targets, sep='\s+', comment='#')
    target_df = target_df[target_df['Plot'].isin(['y', 'Y'])]
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
datespan = Time(dateend) - Time(datestart)
datemid = (Time(datestart) + datespan * 0.45).value

date_duration = np.arange(Time(datestart), Time(dateend) + 1 * u.d, 3 * u.d)

# Date Intervals for Midnight at the Observatory
date_intervals = [(localdate - OBS_TIMEZONE * u.hour).datetime for localdate in date_duration]
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The Trajectories Of Objects To Be Observed
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(18, 13))
ax = fig.add_subplot(111)

for index, value in enumerate(list_values):
    if len(value.split()) >= 3:
        ObjectToObs(name=value.split()[-3], ra=value.split()[-2], dec=value.split()[-1],
                    plot_ax=ax).plot_objtrack(date_intervals)
    elif len(value.split()) == 2:
        ObjectToObs(name='Object ' + str(int(index) + 1), ra=value.split()[-2], dec=value.split()[-1],
                    plot_ax=ax).plot_objtrack(date_intervals)
    else:
        print("ERROR : Both RA & DEC for the Object '{}' need to be specified".format(str(int(index) + 1)))
        continue

plot_obsplan(ax_obj=ax)
# ------------------------------------------------------------------------------------------------------------------- #

