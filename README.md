# NightSkyPlan: Night Sky Planner #

[![Build Status](https://img.shields.io/badge/release-0.1-orange)](https://github.com/sPaMFouR/NightSkyPlan)
[![Python 3.7](https://img.shields.io/badge/python-3.7.2-brightgreen.svg)](https://www.python.org/downloads/release/python-372/)

NightSkyPlan is a python script to help astronomers plan night sky observations from ground based obsbervatories. Observatory details can be added to the file 'TelescopeList.dat' and chosen when the script runs.

Features:
1) Determines observability of targets specified in 'TargetsList.dat' from the chosen observatory.
2) Observability chart also shows altitude, airmass, moon illumination and moon-separation from the specified targets.

The code is still under development.

Authors
-------

* **Avinash Singh** (IIA, Bengaluru)

Requirements
-------

- numpy
- astropy
- easygui
- pandas
- datetime
- matplotlib
