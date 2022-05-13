# NightSkyPlan: Night Sky Planner #

[![Build Status](https://img.shields.io/badge/release-0.2-red)](https://github.com/sPaMFouR/NightSkyPlan)
[![Python 3.9.5](https://img.shields.io/badge/python3.9.5-nightskyplan-green)](https://www.python.org/downloads/release/python-395/)

The directory contains python scipts to aid with planning of night sky observations, planning for proposal cycles, designing a proposal cycle with details about the duration of night and length of the night etc, and also extracting information of moon phase and moon separation for observations of photometric standards. The directory currently hosts the following codes:

## 1) NightSkyPlan.py
Required File(s): [`TelescopeList.dat`], [`TargetList.dat(optional)`]<br />
Output File(s): `NightSkyPlan_DATE.pdf`<br />

'NightSkyPlan.py' is a python script designed to assist in planning night sky observations from ground based observatories/sites. Observatory/site details can be added to the file 'TelescopeList.dat' and chosen when the script runs. The code determines observability of targets specified in 'TargetsList.dat' from the chosen observatory. The observability chart also shows altitude, airmass, moon phase and moon-separation from the specified targets.

## 2) YearlyPlan.py
Required File(s): [`TelescopeList.dat`], [`TargetList.dat(optional)`]<br />
Output File(s): `YearlyPlan_StartDATEToEndDATE.pdf`<br />

'YearlyPlan.py' is a python script designed to assist in planning long term observations (specifically with regards to observational proposals) from ground based observatories/sites. The planner can be given a custom date range to plan object long term observability. Observatory/site details can be added to the file 'TelescopeList.dat' and chosen when the script runs. The code determines the observability of targets specified in 'TargetsList.dat' from the chosen observatory. The yearly observability chart shows altitude, airmass, telescope zenith and telescope horizon.

## 3) CalcTwilightTime.py
Required File(s): [`TelescopeList.dat`]<br />
Output File(s): `NightDuration_StartDATEToEndDATE.pdf`, `TwilightTimes_StartDATEToEndDATE.asc`<br />

'CalcTwilightTime.py' is a python script designed to aid in designing a proposal cycle with details about sunset, sunrise, twilight times, duration of night and moon phase.

## 3) CalcMoonAnglePhase.py
Required File(s): [`TelescopeList.dat`], [`DateList.dat(optional)`]<br />
Output File(s): `MoonPhaseAngle.asc`<br />

'CalcMoonAnglePhase.py' is a python script designed to aid in computing details about moon phase and moon separation for observations of photometric standards.

[`TelescopeList.dat`]: https://github.com/sPaMFouR/NightSkyPlan/blob/master/TelescopeList.dat
[`TargetList.dat(optional)`]: https://github.com/sPaMFouR/NightSkyPlan/blob/master/TargetList.dat
[`DateList.dat(optional)`]: https://github.com/sPaMFouR/NightSkyPlan/blob/master/DateList.dat

The code is still under development.

Requirements:
-------

- ephem
- numpy
- pandas
- astropy
- easygui
- datetime
- matplotlib
