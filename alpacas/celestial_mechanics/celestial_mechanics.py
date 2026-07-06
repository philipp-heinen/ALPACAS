"""
This module provides tools for calculating the declination and the equation of time for a specified date.
"""

import numpy as np
from numpy import cos, sin, tan, arccos, arcsin, arctan, arctan2, sqrt, pi
import datetime as dt
import Constants as const


def solve(f, a, b, tol=1e-12, args=()):
    """
    Implements a simple regula falsi solver.
    """
    A = a
    B = b
    while True:
        C = (A*f(B, *args)-B*f(A, *args))/(f(B, *args)-f(A, *args))
        if np.abs(f(C, *args)) < tol:
            break
        if np.sign(f(A, *args)) == np.sign(f(C, *args)):
            A = C
        else:
            B = C 
    return C


def myarctan2(x, y):
    temp = arctan2(x, y)
    if temp < 0:
        return temp+2*pi
    else:
        return temp
    
def mapping(x):
    if x > 0.5:
        return x-1
    elif x < -0.5:
        return x+1
    else: 
        return x

def kepler_equation(xi, M, eps):
    """
    xi: eccentric anomaly.
    M: mean anomaly.
    eps: eccentricity.
    """
    return xi-eps*sin(xi)-M

def solve_kepler_equation(M, eps, exact=True):
    """
    xi: eccentric anomaly.
    M: mean anomaly.
    eps: eccentricity.
    exact: if true, solve Kepler's equation exactly, otherwise by Taylor expansion.
    """
    if exact:
        return solve(kepler_equation, 0, 2*pi, args=(M, eps))
    else:
        return M+eps*sin(M)+0.5*eps**2*sin(2*M)
    
def kepler_trajectory(M, eps, exact=True):
    """
    Computes the position of a Keplerian body at mean anomaly M.
    The length of the long half axis of the Kepler ellipse is set to 1.
    
    Parameters
    ----------
    M: mean anomaly.
    eps : eccentricity.
    exact: if true, solve Kepler's equation exactly, otherwise by Taylor expansion.

    Returns
    -------
    The x and y coordinate of the body. The origin is at the position of the sun, x is along the long axis of the ellipse, y along the short one.

    """
    xi = solve_kepler_equation(M, eps, exact)
    return np.array([cos(xi)-eps, sqrt(1-eps**2)*sin(xi)])

def convert_true_to_mean_anomaly(phi, eps):
    """
    Helper function to convert true (phi) to mean anomaly (M).

    """
    if 0 <= phi <= pi:
        xi = arccos((eps+cos(phi))/(1+eps*cos(phi)))
    elif pi < phi <= 2*pi:
        xi = 2*pi-arccos((eps+cos(phi))/(1+eps*cos(phi)))
    else:
        xi = np.nan
        
    return xi-eps*sin(xi)

def sun_coordinate(t, exact=True):
    """
    Computes the position of the sun in the equatorial coordinate system as function of time t since the passage of the spring equinox, expressed as fraction of the total orbital period (one year).

    Parameters
    ----------
    t : time t since the passage of the spring equinox, expressed as fraction of the total orbital period (one year).
    exact: if true, solve Kepler's equation exactly, otherwise by Taylor expansion.

    Returns
    -------
    Cartsian coordinates of the sun position in the equatorial coordinate system.

    """
    kt = kepler_trajectory(2*pi*t-convert_true_to_mean_anomaly(const.LONGITUDE_OF_PERHELION/180*pi, const.ECCENTRICITY), const.ECCENTRICITY, exact)
    alpha = const.ECLIPTIC/180*pi
    beta = const.LONGITUDE_OF_PERHELION/180*pi
    v1 = np.array([cos(beta), cos(alpha)*sin(beta), sin(alpha)*sin(beta)])       #vector of long half axis in equatorial coordinate system
    v2 = np.array([-sin(beta), cos(alpha)*cos(beta), sin(alpha)*cos(beta)])       #vector of short half axis in equatorial coordinate system
    return kt[0]*v1+kt[1]*v2
    
def declination(t, exact=True):
    """
    Computes the declination of the sun as function of time t since the passage of the spring equinox, expressed as fraction of the total orbital period (one year).
    
    Parameters
    ----------
    t : time t since the passage of the spring equinox, expressed as fraction of the total orbital period (one year).
    exact: if true, solve Kepler's equation exactly, otherwise by Taylor expansion.

    Returns
    -------
    Declination.

    """

    sc = sun_coordinate(t, exact)
    return arctan(sc[2]/sqrt(sc[0]**2+sc[1]**2))

def declination_date(month, day, leap_year=False, exact=True):
    """
    Computes the declination of the sun as function of the date. 
    Sub-day precision is not provided, as this is not necessary for sundials in general. In case of need, the declination(t) function can be used.
    
    
    Parameters
    ----------
    month: month of the date.
    day: day of the date.
    leap_year: if true, calculate for leap year.
    exact: if true, solve Kepler's equation exactly, otherwise by Taylor expansion.

    Returns
    -------
    Declination.

    """
    if leap_year:
        t = (dt.datetime(4, month, day)-dt.datetime(4, const.EQUINOX[0], const.EQUINOX[1])).days/const.YEAR
    else:
        t = (dt.datetime(1, month, day)-dt.datetime(1, const.EQUINOX[0], const.EQUINOX[1])).days/const.YEAR
        
    return declination(t, exact=exact)

def equation_of_time(t, eot_at_equinox=None, exact=True):
    """
    Computes the equation of time in minutes as function of the time t since the passage of the spring equinox, expressed as fraction of the total orbital period (one year).
    Parameters
    ----------
    t : time t since the passage of the spring equinox, expressed as fraction of the total orbital period (one year).
    exact: if true, solve Kepler's equation exactly, otherwise by Taylor expansion.
    eot_at_equinox: Optionally, a value for the equation of time at the spring equinox in minutes can be set by hand.
                    If None, the position of the mean sun in the equatorial plane at t=0 is set such that it passes the equinox point at the same moment as the auxiliary mean sun moving in the ecliptic plane,
                    where the position of the latter at t=0 is set such that it passes the periapsis at the same moment as the true sun.
                    This convention makes the temporal integral over the equation of time vanish to high precision for small eccentrity and yields sufficient precision for all sundial-related purposes.

    Returns
    -------
    Equation of time in minutes.

    """

    sc = sun_coordinate(t, exact)
    if eot_at_equinox is None:
        offset = (convert_true_to_mean_anomaly(const.LONGITUDE_OF_PERHELION/180*pi, const.ECCENTRICITY)-const.LONGITUDE_OF_PERHELION/180*pi)/(2*pi)
    else:
        offset = -eot_at_equinox*60/const.SIDDAY
    return mapping(t-offset-myarctan2(sc[1], sc[0])/(2*pi))*const.SIDDAY/60


def equation_of_time_date(month, day, leap_year=False, eot_at_equinox=None, exact=True):
    """
    Computes the equation of time in minutes as function of the date.
    Sub-day precision is not provided, as this is not necessary for sundials in general. In case of need, the equation_of_time(t) function can be used.    
    Parameters
    ----------
    month: month of the date.
    day: day of the date.
    leap_year: if true, calculate for leap year.
    exact: if true, solve Kepler's equation exactly, otherwise by Taylor expansion.
    eot_at_equinox: Optionally, a value for the equation of time at the spring equinox in minutes can be set by hand.
                    If None, the position of the mean sun in the equatorial plane at t=0 is set such that it passes the equinox point at the same moment as the auxiliary mean sun moving in the ecliptic plane,
                    where the position of the latter at t=0 is set such that it passes the periapsis at the same moment as the true sun.
                    This convention makes the temporal integral over the equation of time vanish to high precision for small eccentrity and yields sufficient precision for all sundial-related purposes.

    Returns
    -------
    Equation of time in minutes.

    """

    if leap_year:
        t = (dt.datetime(4, month, day)-dt.datetime(4, const.EQUINOX[0], const.EQUINOX[1])).days/const.YEAR
    else:
        t = (dt.datetime(1, month, day)-dt.datetime(1, const.EQUINOX[0], const.EQUINOX[1])).days/const.YEAR
        
    return equation_of_time(t, eot_at_equinox=eot_at_equinox, exact=exact)


def find_solstice(which, exact=True):
    """
    Small helper function that computes the time t at which the solstices occur, expressed at fraction of the orbital period (one year), starting at the spring equinox.
    """
    eps = 0.0001
    if which == "summer":
        return solve(lambda x: (declination(x+eps/2, exact=exact)-declination(x-eps/2, exact=exact))/eps, 0, 0.5)
    elif which == "winter":
        return solve(lambda x: (declination(x+eps/2, exact=exact)-declination(x-eps/2, exact=exact))/eps, 0.5, 1)
