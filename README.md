# NightSkyPlan: Night Sky Planner #

[![Build Status](https://img.shields.io/badge/release-0.1-orange)](https://github.com/sPaMFouR/NightSkyPlan)
[![Python 3.7](https://img.shields.io/badge/python--3.7.4-nightskyplanner-brightgreen)](https://www.python.org/downloads/release/python-374/)

The directory contains python scipts to aid with planning of night sky observations, planning for proposal cycles, designing a proposal cycle with details about the duration of night and length of the night etc, and also extracting information of moon phase and moon separation for observations of photometric standards. The directory currently hosts the following codes:

## 1) NightSkyPlan.py
Required File(s): [`TelescopeList.dat`], [`TargetList.dat(optional)`]<br />
Output File(s): `NightSkyPlan_DATE.pdf`<br />

'NightSkyPlan.py' is a python script designed to assist in planning night sky observations from ground based observatories/sites. Observatory/site details can be added to the file 'TelescopeList.dat' and chosen when the script runs. The code determines observability of targets specified in 'TargetsList.dat' from the chosen observatory. The observability chart also shows altitude, airmass, moon phase and moon-separation from the specified targets.

## 2) CalcTwilightTime.py
Required File(s): [`TelescopeList.dat`]<br />
Output File(s): `NightDuration_StartDATEToEndDATE.pdf`, `TwilightTimes_StartDATEToEndDATE.asc`<br />

'CalcTwilightTime.py' is a python script designed to aid in designing a proposal cycle with details about sunset, sunrise, twilight times, duration of night and moon phase.

## 3) CalcMoonAnglePhase.py
Required File(s): [`TelescopeList.dat`], [`DateList.dat(optional)`]<br />
Output File(s): `MoonPhaseAngle.asc`<br />

'CalcMoonAnglePhase.py' is a python script designed to aid in computing details about moon phase and moon separation for observations of photometric standards.
The code is still under development.

[`TelescopeList.dat`]: https://github.com/sPaMFouR/NightSkyPlan/TelescopeList.dat
[`TargetList.dat(optional)`]: https://github.com/sPaMFouR/NightSkyPlan/TargetList.dat
[`DateList.dat(optional)`]: https://github.com/sPaMFouR/NightSkyPlan/DateList.dat


Requirements
-------

- ephem
- numpy
- astropy
- easygui
- pandas
- datetime
- matplotlib
