# ALPACAS
ALPACAS (A Library for Plotting And Computing All kinds of Sundials) is a Python library for sundials. As the name suggests, it possesses two main functionalities. On the one hand, one can use ALPACAS to plot custom-made sundials, without any expert knowledge. These sundial plots can be printed out and e.g. glued or otherwise fixated on a piece of wood or cardboard. Then one only needs to add a small vertical stick (e.g. a nail) serving as nodus and has created a fully functional sundial! This sundial can be customized to indicate the time in several different systems, as well as the date. 

On the other hand, for experts ALPACAS provides many functions for performing sundial-related computations. It allows to compute the sun position in local coordinates as function of latitude, declination and hour angle, the position of the shadow cast by a nodus onto a sundial, the angles of the hour lines cast by a polar gnomon, the hour angle of sunset/sunrise and some more quantities. Furthermore, it has a celestial mechanics module that provides the declination and equation of time for a given date (taking as input only earth's orbital parameters, which are provided in a module with constants). 

ALPACAS supports sundials of the following orientations (i.e. orientation of the dial plane): horizontal, vertical, polar, equatorial and of arbitrary orientation provided by inclination angles. 

## Requirements and Installation
ALPACAS needs only the ``numpy`` and ``matplotlib`` packages. For using the library, just download this repository and then make it known to your Python interpreter, e.g. under Linux by running 
```
export PYTHONPATH=<path to the repository>/ALPACAS:$PYTHONPATH
```
or adding this line to your ``~/.bashrc`` file to make this permanent.

## Plotting functionalities
The plotting functionalities are contained in the core module ``Alpacas``. First import this module:
```
import Alpacas as alp
```
Now we create a ``Sundial`` object:
```
sundial = alp.Sundial(latitude=50.733, longitude=7.104, orientation="horizontal")
```
The latitude and longitude of the position where the sundial will be located must be provided in degrees, with latitudes in the northern hemisphere counting as positive and in the southern hemisphere as negative, as well as longitudes east of Greenwich counted as positive and west of it as negative. The ``orientation`` parameter (which defaults to "horizontal") determines the orientation of the sundial. The following options are available:

+ ``"horizontal"`` (the dial is oriented horizontally)
+ ``"vertical"`` (the dial is oriented vertically, facing exactly south (northern hemisphere)/north (southern hemisphere))
+ ``"polar"`` (the dial is oriented parallel to the polar axis)
+ ``"equatorial_upper_side"`` (the dial is oriented parallel to the equator plane, and located on the "upper" side)
+ ``"equatorial_lower_side"`` (the dial is oriented parallel to the equator plane, and located on the "lower" side)
+ ``[alpha, beta]`` (a list of two inclination angles of the dial plane provided in degrees)

A horizontal sundial has the advantage that it can display the time during the whole day, from sunrise to sunset. It is also easy to realize and fairly common. Vertical sundials are also very widespread as they can be realized on south-facing (southern hemisphere: north-facing) walls of buildings. The other types are much rarer and their implementation in ALPACAS is provided rather for experts. 

First, we initialize the sundial plot:
```
sundial.init_dial_plot(xsize=10, ysize=10, cm_per_unit=2)
```
``xsize`` and ``ysize`` determine the *relative* size of the plotted dial, where we set the height of the nodus to 1. In many cases, the default ``xsize=10`` and ``ysize=10`` is a good choice. ``cm_per_unit`` sets how the relative units in which the height of the nodus amounts to 1 translate into centimeters. This depends on how large one wants the sundial to be.

Now we can populate the sundial with different elements. For example, before adding any lines indicating time, it is useful to add a marker for where the nodus has to be located (by default, it is placed in the center of the plot). 
```
sundial.add_nodus_pos(style={"color":"black", "ms":4})
```
The function (optionally) takes a ``style`` argument which must be provided as a dictionary setting color and markersize (ms) of the marker.

Next, we can add lines for indicating the time. For example, one of the most common times displayed on sundials is the apparent zone time, which deviates from the time displayed by clocks only by a few minutes (the equation of time). To add lines for apparent zone time, we can e.g. write:
```
sundial.add_apparent_zone_time(which="hourly", timezone=1, gnomon_type="nodus", 
                               style={"color":"blue", "lw":2, "linestyle": "-"}, 
                               which_labels="roman_numerals", 
                               labelstyle={"side":"left", "position":0.8, "position_type":"relative", "tightness":0.25, "fontsize":8, "color":"blue", "rotate": False})
```
``which`` determines for which hours to draw lines. It can be "hourly", "half-hourly", "two-hourly" or a costum list of hours. The ``timezone`` parameter must be set to the difference of the zonal time from UTC time in hours. It can also be set to "closest", in which case the closest time zone to the provided longitude is taken. This can be wrong, however, for many places on earth. ``gnomon_type`` should be "nodus" if a vertical stick (a nodus) is used, and "polar" if a polar gnomon, i.e. a stick parallel to the polar axis, is used. ``which_labels`` can be set to "roman_numerals", "arabic_numerals", "no_labels" or a custom list of labels (whose length then must match with that of the list provided for ``which``). ``style`` and ``labelstyle`` optionally provide the possibility to highly customize the appearance of the hour lines and their labels.

Other possible types of time systems supported by ALPACAS are 

+ apparent solar time (time system in which the sun reaches the zenith at 12:00): ``add_apparent_solar_time``
+ mean solar time (as mean zone time but not uniform across time zone): ``add_mean_solar_time``
+ mean zone time (time indicated by clocks): ``add_mean_zonal_time``
+ Babylonian hours (hours since sunrise): ``add_babylonian_hours``
+ Italian hours (hours since sunset): ``add_italian_hours``
+ temporal hours (period from sunset to sunrise is divided into twelve hours of equal length): ``add_temporal_hours``

These functions are well documented within the ALPACAS code.

Additionally to the lines indicating the time of the day, we can also add curves indicating a particular date. For example, in order to have the equinoxes and solstices displayed, we can write
```
sundial.add_date_line(date="equinox", label="Equinox")
sundial.add_date_line(date="summer_solstice", label="Summer solstice")
sundial.add_date_line(date="winter_solstice", label="Winter solstice")
```
Dates other than the equinox and solstices can be provided as a ``[month, day]``. If desired, the appearence of the line and the labels can again be customized by ``style`` and ``labelstyle`` arguments.

Sometimes it can be nice to add some text to the dial. This can be done as 
```
sundial.add_text("South", 0, -4, fontsize=12, color="black")
```
where the position of the text is given in coordinates in relative units.

Finally, if we are happy with the sundial we have created, we can save it to a file:
```
sundial.save_dial_plot("Sundial.png", resolution=300)
```
The ``resolution`` is provided in dpi.

## Computation functionalities
The computation functionalities are contained in the ``BasicFunctions`` module. For example, compute the position of the shadow a nodus casts onto a horizontal sundial at hour angle $\tau=20^\circ$, declination $\delta=10^\circ$ and latitude $\phi=50^\circ$,:
```
>>> import BasicFunctions as bf
>>> import numpy as np
>>> bf.nodus_horizontal(tau=np.deg2rad(20), delta=np.deg2rad(10), phi=np.deg2rad(50))
array([0.46275402, 0.82060332])
```
Note that the functions in ``BasicFunctions`` in general take angles in radian and not in degree.
The single functions provided are documented well within the code, so we refrain from a detailed documentation here.

The module ``CelestialMechanics`` provides essential functionalities such as the computation of the solar declination and equation of time for a specified date. For example, compute the declination in radian and the equation of time in minutes for 1 January:
```
>>> import CelestialMechanics as cm
>>> cm.declination_date(month=1, day=1, leap_year=False)
-0.40081282075217145
>>> cm.equation_of_time_date(month=1, day=1, leap_year=False)
-3.6196350209321593
```
Also the functions in this module are documented well within the code.
