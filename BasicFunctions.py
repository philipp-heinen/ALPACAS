"""
This module provides useful elementary functions.
"""

import numpy as np
from numpy import sin, cos, tan, arccos, arcsin, arctan, arctan2, sqrt

def sun_vector(tau, delta, phi):
    """
    Calculates the normalized sun ray vector in the local coordinate system.
    
    Parameters
    ----------
    tau : hour angle.
    delta : declination.
    phi: latitude.

    Returns
    -------
    A 3D vector representing the sun ray in local coordinates. 
    The x coordinate is along the west-east direction, the y coordinate along the south-north direction and the z coordinate is the vertical direction.

    """
    return np.array([cos(delta)*sin(tau),
            sin(phi)*cos(delta)*cos(tau)-cos(phi)*sin(delta),
            -cos(phi)*cos(delta)*cos(tau)-sin(phi)*sin(delta)])


def nodus_horizontal(tau, delta, phi):
    """
    Calculates the position of the shadow that a nodus (a point-like gnomon) casts on a horizontal sundial.
    
    Parameters
    ----------
    tau : hour angle.
    delta : declination.
    phi: latitude.

    Returns
    -------
    A vector with the x and y coordinate of the shadow point on the dial as entries.
    The x coordinate is along the west-east (east-west on the southern hemisphere) direction and the y coordinate along the south-north (north-south on the southern hemisphere) direction.
    The origin (0,0) is placed under the nodus, whose length is set to 1.

    """
    sv=sun_vector(tau, delta, phi)
    if sv[2]<0:
        if phi>=0:
            return np.array([-sv[0]/sv[2],-sv[1]/sv[2]])
        else:
            return np.array([sv[0]/sv[2],sv[1]/sv[2]])
    else:
        return np.array([np.nan, np.nan])
    
    
def nodus_vertical(tau, delta, phi):
    """
    Calculates the position of the shadow that a nodus (a point-like gnomon) casts on a vertical sundial facing excatly to the south (or north on the southern hemisphere).
    
    Parameters
    ----------
    tau : hour angle.
    delta : declination.
    phi: latitude.

    Returns
    -------
    A vector with the x and y coordinate of the shadow point on the dial as entries.
    The x coordinate is along the west-east (east-west on the southern hemisphere) direction and the y coordinate along the vertical direction.
    The origin (0,0) is placed under the nodus (the connection between nodus and the origin is perpendicular to the dial plane), whose length is set to 1.

    """
    sv=sun_vector(tau,delta,phi)
    if phi>=0:
        if sv[1]>0 and sv[2]<0:
            return np.array([sv[0]/sv[1],sv[2]/sv[1]])
        else: 
            return np.array([np.nan, np.nan])
    else:
        if sv[1]<0 and sv[2]<0:
            return np.array([sv[0]/sv[1],-sv[2]/sv[1]])
        else: 
            return np.array([np.nan, np.nan])        


def nodus_equatorial(tau, delta, phi, upper_side=True):
    """
    Calculates the position of the shadow that a nodus (a point-like gnomon) casts on an equatorial sundial (i.e. the dial is parallel to the equatorial plane).
    
    Parameters
    ----------
    tau : hour angle.
    delta : declination.
    phi: latitude.
    upper_side: if true, calculate for the upper side of the horizontal dial (the one the sun shines onto during summer), otherwise for the lower side. 

    Returns
    -------
    A vector with the x and y coordinate of the shadow point on the dial as entries.
    The x coordinate is along the west-east (east-west on the southern hemisphere) direction and the y coordinate along the direction perpendicular to that.
    The origin (0,0) is placed under the nodus (the connection between nodus and the origin is perpendicular to the dial plane), whose length is set to 1.

    """
    sv=sun_vector(tau, delta, phi)
    if phi>=0:
        if upper_side:
            if delta>0 and sv[2]<0:
                return np.array([sin(tau)/tan(delta),cos(tau)/tan(delta)])
            else: 
                return np.array([np.nan,np.nan])
        else:
            if delta<0 and sv[2]<0:
                return np.array([-sin(tau)/tan(delta),cos(tau)/tan(delta)])
            else: 
                return np.array([np.nan,np.nan])
    else:
        if upper_side:
            if delta<0 and sv[2]<0: 
                return np.array([sin(tau)/tan(delta),-cos(tau)/tan(delta)])
            else: 
                return np.array([np.nan,np.nan])
        else:
            if delta>0 and sv[2]<0:
                return np.array([-sin(tau)/tan(delta),-cos(tau)/tan(delta)])
            else: 
                return np.array([np.nan,np.nan])


def nodus_polar(tau, delta, phi):
    """
    Calculates the position of the shadow that a nodus (a point-like gnomon) casts on a vertical sundial facing excatly to the south.
    
    Parameters
    ----------
    tau : hour angle.
    delta : declination.
    phi: latitude.

    Returns
    -------
    A vector with the x and y coordinate of the shadow point as entries.
    The x coordinate is along the west-east (east-west on the southern hemisphere) direction and the y coordinate along the direction perpendicular to that.
    The origin (0,0) is placed under the nodus (the connection between nodus and the origin is perpendicular to the dial plane), whose length is set to 1.

    """
    if -np.pi/2<tau<np.pi/2 and sun_vector(tau,delta,phi)[2]<0:
        if phi>=0:
            return np.array([tan(tau),-tan(delta)])
        else:
            return np.array([-tan(tau),tan(delta)])
    else: 
        return np.array([np.nan,np.nan])


def nodus_arbitrary_orientation(tau, delta, phi, alpha, beta):
    """
    Calculates the position of the shadow that a nodus (a point-like gnomon) casts on a sundial with arbitrary orientation.
    The normal vector vector to the dial plane is 
    n=(-sin(alpha)*sin(beta),sin(alpha)*cos(beta),cos(alpha)) 
    where the three coordinates are along the west-east, south-north and vertical direction, respectively.
    For example, a vertical dial not facing exactly to the south but being misaligned by an angle of epsilon would be represented by alpha=-pi/2 (alpha=pi/2 on the southern hemisphere), beta=epsilon.
    We define the two vectors spanning the dial plane as 
    v1=(cos(beta),sin(beta),0)
    v2=(-cos(alpha)*sin(beta),cos(alpha)*cos(beta),-sin(alpha))
    This convention is chosen such that v1 and v2 are orthogonal to each other and v1 lies in the horizontal plane.
    
    
    Parameters
    ----------
    tau : hour angle.
    delta : declination.
    phi: latitude.
    alpha: one of the orientation angles.
    betaa: one of the orientation angles.
    Returns
    -------
    A vector with the x and y coordinate of the shadow point as entries.
    The x coordinate is along v1 and the y coordinate along v2.
    The origin (0,0) is placed under the nodus (the connection between nodus and the origin is perpendicular to the dial plane, i.e. the nodus is parallel to n), whose length is set to 1.

    """
    sv=sun_vector(tau,delta,phi)
    n=np.array([-sin(alpha)*sin(beta),sin(alpha)*cos(beta),cos(alpha)])
    v1=([cos(beta),sin(beta),0])
    v2=np.array([-cos(alpha)*sin(beta),cos(alpha)*cos(beta),-sin(alpha)])
    M=np.transpose(np.array([-sv,v1,v2]))
    sol=np.linalg.solve(M,n)
    if sol[0]>0 and sv[2]<0:
        return np.array([sol[1],sol[2]])
    else: 
        return np.array([np.nan, np.nan])    
   
    
def position_gnomon_horizontal(phi):
    """
    Calculates the position where a gnomon polar intersects with the dial plane if it passes through the nodus, which is placed one unit length above the coordinate origin (in the dial coordinates).
    A horizontal orientation of the dial is assumed. Coordinates are defined as in the other functions. 

    Parameters
    ----------
    phi : latitude.

    Returns
    -------
    Coordinates of the intersection point.

    """
    if phi>=0:
        return np.array([0,-1/tan(phi)])
    else:
        return np.array([0,1/tan(phi)])


def position_gnomon_vertical(phi):
    """
    Calculates the position where a gnomon polar intersects with the dial plane if it passes through the nodus, which is placed one unit length above the coordinate origin (in the dial coordinates).
    A vertical orientation of the dial is assumed. Coordinates are defined as in the other functions. 

    Parameters
    ----------
    phi : latitude.

    Returns
    -------
    Coordinates of the intersection point.

    """
    if phi>=0:
        return np.array([0,tan(phi)])
    else:
        return np.array([0,-tan(phi)])


def position_gnomon_arbitrary_orientation(phi, alpha, beta):
    """
    Calculates the position where a gnomon polar intersects with the dial plane if it passes through the nodus, which is placed one unit length above the coordinate origin (in the dial coordinates).
    The orientation of the dial is described by the angles alpha and beta as in nodus_arbitrary_orientation(). Coordinates are defined as in the other functions. 

    Parameters
    ----------
    phi : latitude.
    alpha: one of the orientation angles.
    betaa: one of the orientation angles.
    Returns
    -------
    Coordinates of the intersection point.

    """
    pg=np.array([0,cos(phi),sin(phi)])          #vector of the polar gnomon
    n=np.array([-sin(alpha)*sin(beta),sin(alpha)*cos(beta),cos(alpha)])
    v1=([cos(beta),sin(beta),0])
    v2=np.array([-cos(alpha)*sin(beta),cos(alpha)*cos(beta),-sin(alpha)])
    M=np.transpose(np.array([-pg,v1,v2]))
    sol=np.linalg.solve(M,n)
    return np.array([sol[1],sol[2]])


def gnomon_horizontal(tau, phi):
    """ 
    Caculates the angle between the shadow line that a polar gnomon (i.e. a gnomon parallel to the polar axis) casts at hour angle tau and the south-north (north-south on the southern hemisphere) direction, on a horizontal sundial.

    Parameters
    ----------
    tau : hour angle.
    phi : latitude.
    """    
    return arctan2(sin(tau)*sin(phi), cos(tau))


def gnomon_vertical(tau, phi):
    """ 
    Caculates the angle between the shadow line that a polar gnomon (i.e. a gnomon parallel to the polar axis) casts at hour angle tau and the the vertical direction, on a vertical sundial.

    Parameters
    ----------
    tau : hour angle.
    phi : latitude.
    """     
    if -np.pi/2<tau<np.pi/2:
        if phi>=0:
            return arctan2(sin(tau)*cos(phi), -cos(tau))
        else:
            return arctan2(-sin(tau)*cos(phi), -cos(tau))
    else:
        return np.nan


def gnomon_equatorial(tau, phi, upper_side=True):
    """ 
    Caculates the angle between the shadow line that a polar gnomon (i.e. a gnomon parallel to the polar axis) casts at hour angle tau and the direction perpendicular to the west-east (east-west on the southern hemisphere) direction, on an equatorial sundial.

    Parameters
    ----------
    tau : hour angle.
    phi : latitude.
    """    
    if phi>=0:
        if upper_side:
            return tau
        else:
            if tau>0:
                return np.pi-tau
            else:
                return -np.pi-tau
    else:
        if upper_side:
            return -tau
        else:
            if tau<0:
                return np.pi+tau
            else:
                return -np.pi+tau


def gnomon_arbitrary_orientation(tau, phi, alpha, beta):
    """ 
    Caculates the angle between the shadow line that a polar gnomon (i.e. a gnomon parallel to the polar axis) casts at hour angle tau and the direction of v2=(-cos(alpha)*sin(beta),cos(alpha)*cos(beta),-sin(alpha)), on a sundial with arbitrary orientation defined by angles alpha and beta.

    Parameters
    ----------
    tau : hour angle.
    phi : latitude.
    alpha: one of the orientation angles.
    betaa: one of the orientation angles.
    """     
    sv=sun_vector(tau, 0, phi)
    pg=np.array([0,cos(phi),sin(phi)])
    aux1=np.cross(sv,pg)
    n=np.array([-sin(alpha)*sin(beta),sin(alpha)*cos(beta),cos(alpha)])
    v1=([cos(beta),sin(beta),0])
    v2=np.array([-cos(alpha)*sin(beta),cos(alpha)*cos(beta),-sin(alpha)])
    aux2=np.cross(n,aux1)
    if np.dot(aux2,sv)<0:
        aux2=-aux2
    return arctan2(np.dot(aux2,v1), np.dot(aux2,v2))


def sunset(delta, phi):
    """
    Computes the hour angle at which sunset occurs (without taking into account atmospherical effects).

    Parameters
    ----------
    delta : declination.
    phi: latitude.
    
    Returns
    -------
    Hour angle of sunset.

    """
    if -tan(phi)*tan(delta)<=-1:
        return np.pi
    elif -tan(phi)*tan(delta)>=1:
        return 0
    else:
        return arccos(-tan(phi)*tan(delta))


def sunset_declination(tau, phi):
    """
    Helper function that computes the declination that is necessary for sunset to occur at a certain hour angle tau
    """
    return arctan(-cos(tau)/tan(phi))
