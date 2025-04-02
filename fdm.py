"""Module Finite Difference Method describe derivatives used in numerical method (FDM)

Note:
    The file contains a set of procedures allowing calculations using the finite difference method in cylindrical space (2D)
    temperature change at any point in time. The procedures were derived for an asymmetric mesh.

    .. code-block:: none

         t[n] - temperature in node n, d[n] - distance in the coordinates Z and R, between node 4 and node n,
         a(k) - temperature equalization factor [m^2/s] for material k
        
                     t[0]                         Z ^
                      O                             |
                      |     (0)                     |
                (3)   |d[0]                         |
                      |                             0----------> R
         t[3] O------t[4]-----------O t[1]
                d[3]  |        d[1]
                      |     (1)
                (2)   | d[2]
                      |
                      O
                     t[2]
"""

def d2tdr2(t, dx, a):
    """Second derivative of temperature by radius

    Args:
        t: temperature vector [°C]
        dx: distance from the central node [m]
        a: temperature equalization factor vector [m^2/s] 

    Returns: 
        second derivative of temperature by radius

    """
    out = (
        2.0
        * (
            (
                (t[1] - t[4]) * ((dx[0] * a[0] + dx[2] * a[1]))
                - (dx[1] * (t[4] - t[3]) * (dx[0] * a[3] + dx[2] * a[2]))
                / dx[3]
            )
        )
        / ((dx[0] + dx[2]) * ((dx[1]) ** 2.0 + dx[1] * dx[3]))
    )
    return out


def dtdr(t, dx, a, r):
    """First derivative of temperature by radius

    Args:
        t: temperature vector [°C]
        dx: distance from the central node [m]
        a: temperature equalization factor vector [m^2/s]
        r: radius counted from the axis of symmetry [m]

    Returns: 
        first derivative of temperature by radius
    """
    out = d2tdr2(t, dx, a) * r + (
        ((t[4] - t[3]) * (dx[0] * a[3] + dx[2] * a[2]))
        + ((dx[1] * dx[3]) / (dx[3] * dx[1] + (dx[1]) ** 2.0))
        * (
            (t[1] - t[4]) * (dx[0] * a[0] + dx[2] * a[1])
            + (t[3] - t[4]) * (dx[0] * a[3] + dx[2] * a[2])
        )
    ) / (dx[3] * (dx[0] + dx[2]))
    return out


def d2tdz2(t, dx, a, r):
    """Second derivative of temperature by depth

    Args:
        t: temperature vector [°C]
        dx: distance from the central node [m]
        a: temperature equalization factor vector [m^2/s]
        r: radius counted from the axis of symmetry [m]

    Returns:
        Second derivative of temperature by depth

    """
    count1 = (
        (t[0] - t[4])
        * (
            ((r**2.0 - (r - dx[3]) ** 2.0) * a[3])
            + ((r + dx[1]) ** 2.0 - r**2.0) * a[0]
        )
        / ((r + dx[1]) ** 2.0 - (r - dx[3]) ** 2.0)
    )
    count2 = (
        (dx[0] / dx[2])
        * (t[2] - t[4])
        * (
            (
                (r**2.0 - (r - dx[3]) ** 2.0) * a[2]
                + ((r + dx[1]) ** 2.0 - r**2.0) * a[1]
            )
            / ((r + dx[1]) ** 2.0 - (r - dx[3]) ** 2.0)
        )
    )
    out = 2.0 * (count1 + count2) / ((dx[2]) ** 2.0 + dx[2] * dx[0])
    return out


def dtdtime(t, dx, a, r):
    """Derivative of temperature by time

    Args:
        t: temperature vector [°C]
        dx: distance from the central node [m]
        a: temperature equalization factor vector [m^2/s]
        r: radius counted from the axis of symmetry [m]

    Returns:
        Derivative of temperature by time
    """
    out = d2tdr2(t, dx, a) + dtdr(t, dx, a, r) / r + d2tdz2(t, dx, a, r)
    return out
