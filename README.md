# NightSkyPlan: Night Sky Planner #

[![Build Status](https://img.shields.io/badge/release-0.1-orange)](https://github.com/sPaMFouR/NightSkyPlan)
[![Python 3.7](https://img.shields.io/badge/python--3.7.4-nightskyplanner-brightgreen)](https://www.python.org/downloads/release/python-374/)

The directory contains python scipts to aid with planning of night sky observations, planning for proposal cycles, designing a proposal cycle with details about the duration of night and length of the night etc, and also extracting information of moon phase and moon separation for observations of photometric standards. The directory currently hosts the following codes:

## 1) NightSkyPlan.py
##### Additional Files Required: TelescopeList.dat, TargetList.dat(optional) #####
##### Output File/s: NightSkyPlan_DATE.pdf #####

'NightSkyPlan.py' is a python script designed to assist in planning night sky observations from ground based observatories/sites. Observatory/site details can be added to the file 'TelescopeList.dat' and chosen when the script runs. The code determines observability of targets specified in 'TargetsList.dat' from the chosen observatory. The observability chart also shows altitude, airmass, moon phase and moon-separation from the specified targets.

## 2) CalcTwilightTime.py
##### Additional Files Required: TelescopeList.dat #####
##### Output File/s: NightDuration_StartDATEToEndDate.pdf, TwilightTimes_StartDATEToEndDate.asc #####

'CalcTwilightTime.py' is a python script designed to aid in designing a proposal cycle with details about sunset, sunrise, twilight times, duration of night and moon phase.

## 3) CalcMoonAnglePhase.py
##### Additional Files Required: TelescopeList.dat, DateList.dat(optional) #####
##### Output File/s: MoonPhaseAngle.asc #####

'CalcMoonAnglePhase.py' is a python script designed to aid in computing details about moon phase and moon separation for observations of photometric standards.
The code is still under development.

Requirements
-------

- ephem
- numpy
- astropy
- easygui
- pandas
- datetime
- matplotlib
