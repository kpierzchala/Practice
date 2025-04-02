"""
This modules contains functions of brines properties depending on temperature, pressure and salinity
"""

import math
import numpy as np
import CoolProp.CoolProp as cp


def c_heat(t, pressure, salt_mass_fraction):
    """Specific heat of brine.

    Args:
        t: temperature in range from 0 to 1000°C
        pressure: presure in range from 1 to 5000 bar
        salt_mass_fraction: molar fraction in range from 0 to 1 NaCl

    Return:
        specific heat of brine [J/kg*K]

    Note:
        source: T. Driesner 2007, Procedure for binary H2O-NaCl solutions

    """
    temperature = t - 273.15
    pbar = pressure / 1e5
    molar_mass_H2O = 18.015
    molar_mass_NaCl = 58.443
    n_NaCl = salt_mass_fraction * 1000 / molar_mass_NaCl
    n_H2O = 1000 / molar_mass_H2O
    xmol = n_NaCl / (n_NaCl + n_H2O)
    xmol1 = 1 - xmol
    xmol12 = xmol1 * xmol1

    q11 = -32.1724 + 0.0621255 * pbar
    q21 = -1.69513 - (4.52781e-4) * pbar - (6.04279e-8) * (pbar**2)
    q22 = 0.0612567 + (1.88082e-5) * pbar

    q1x1 = 47.9048 - (9.36994e-3) * pbar + (6.51059e-6) * (pbar**2)
    q2x1 = 0.241022 + (3.45087e-5) * pbar - (4.28356e-9) * (pbar**2)

    q10 = q1x1
    q20 = 1 - q21 * np.sqrt(q22)
    q12 = -q11 - q10
    q23 = q2x1 - q20 - q21 * np.sqrt(1 + q22)

    q1 = q10 + q11 * xmol1 + q12 * xmol12
    q2 = q20 + (q21 * np.sqrt(xmol + q22)) + q23 * xmol

    tstar_h = q1 + q2 * temperature

    c_heat = q2 * cp.PropsSI(
        "C", "T", tstar_h + 273.15, "P", pressure, "Water"
    )
    return c_heat


def density(t, pressure, salt_mass_fraction):
    """Density of brine

    Args:
        t: temperature in range from 0 to 1000°C
        pressure: pressure in range from 1 to 5000 bar
        salt_mass_fraction: molar fraction in range from 0 to 1 NaCl

    Returns
        density of brine [kg/(m^3)]

    Note:
        source: T.Driesner & C.A. Heinrich 2007
    """
    temperature = t - 273.15
    try:
        salt_molecular_weight = 58.44
        water_molecular_weight = 18.01528

        pbar = pressure / 1.0e5

        xmol = ((salt_mass_fraction * 1000) / salt_molecular_weight) / (
            ((salt_mass_fraction * 1000) / salt_molecular_weight)
            + (((1 - salt_mass_fraction) * 1000) / water_molecular_weight)
        )
        xmol1 = 1.0 - xmol
        xmol12 = xmol1 * xmol1
        brine_molecular_weight = (xmol * salt_molecular_weight) + (
            xmol1 * water_molecular_weight
        )

        n11 = -54.2958 - 45.7623 * np.exp((-9.44785e-4) * pbar)
        n21 = -2.6142 - 0.000239092 * pbar
        n22 = 0.0356828 + (4.37235e-6) * pbar + (2.0566e-9) * (pbar**2)

        n1x1 = (
            330.47
            + 0.942876 * np.sqrt(pbar)
            + 0.0817193 * pbar
            - ((2.47556e-8) * (pbar**2))
            + ((3.45052e-10) * (pbar**3))
        )
        n2x1 = (
            -0.0370751
            + 0.00237723 * np.sqrt(pbar)
            + (5.42049e-5) * pbar
            + (5.84709e-9) * (pbar**2)
            - (5.99373e-13) * (pbar**3)
        )

        n10 = n1x1
        n20 = 1.0 - n21 * np.sqrt(n22)
        n12 = -n11 - n10
        n23 = n2x1 - n20 - n21 * np.sqrt(1.0 + n22)

        n1 = n10 + n11 * xmol1 + n12 * xmol12
        n2 = n20 + n21 * np.sqrt(xmol + n22) + n23 * xmol

        n300 = (7.60664e6) / ((pbar + 472.051) ** 2)
        n301 = -50 - 86.1446 * math.exp((-6.21128e-4) * pbar)
        n302 = 294.318 * math.exp((-5.66735e-3) * pbar)
        n310 = -0.0732761 * math.exp((-2.3772e-3) * pbar) - (5.2948e-5) * pbar
        n311 = -47.2747 + 24.3653 * math.exp((-1.25533e-3) * pbar)
        n312 = -0.278529 - 0.00081381 * pbar

        n31 = n310 * math.exp(n311 * xmol) + n312 * xmol
        n30 = n300 * (math.exp(n301 * xmol) - 1) + n302 * xmol
        D = n30 * math.exp(n31 * temperature)
        tstar_v = n1 + n2 * temperature + D

        if pressure <= cp.PropsSI("Water", "pcrit"):
            ts = cp.PropsSI("T", "P", pressure, "Q", 0, "Water") - 273.15
            extrapolate = tstar_v > ts
        else:
            extrapolate = False

        def polynomial(coeffs, x):
            return sum(c * x**i for i, c in enumerate(coeffs))

        def extrapolation(pressure, pbar, ts, tstar_v, brine_molecular_weight):
            dt = 0.2
            water_molecular_weight = 18.01528

            try:
                props_s = cp.PropsSI(
                    "D", "T", ts + 273.15, "P", pressure, "Water"
                )
                dws = props_s
                vws = 1e3 * water_molecular_weight / dws

                props_s1 = cp.PropsSI(
                    "D", "T", (ts - dt) + 273.15, "P", pressure, "Water"
                )
                dws1 = props_s1
                vws1 = 1e3 * water_molecular_weight / dws1

                dvdt = (vws - vws1) / dt
                logp = np.log(pbar)

                o2 = polynomial(
                    [
                        2.0125e-7 + 3.29977e-9 * np.exp(-4.31279 * logp),
                        -1.17748e-7,
                        7.58009e-8,
                    ],
                    logp,
                )
                ts2 = ts * ts
                o1 = dvdt - 3 * o2 * ts2
                o0 = vws - ts * (o1 + o2 * ts2)
                vb = polynomial([o0, o1, 0.0, o2], tstar_v)

                return 1e3 * brine_molecular_weight / vb

            except Exception as e:
                print(f"An error occurred in extrapolation: {e}")
                return None

        if extrapolate:
            brine_density = extrapolation(
                pressure, pbar, ts, tstar_v, brine_molecular_weight
            )
        else:
            d_H20 = cp.PropsSI(
                "D", "T", tstar_v + 273.15, "P", pressure, "Water"
            )
            vm_H20 = (water_molecular_weight / 1000) / d_H20
            vm_sol = vm_H20
            rho = (brine_molecular_weight / 1000) / vm_sol
        return rho

    except Exception as e:
        print(f"An error occurred: {e}")
        return None


def dyn_viscosity(t, p, s):
    """Dynamic viscosity of brine

    Args:
        t: Temperature in range from 0 to 350°C
        p: Presure in range from 0.1 to 100 MPa
        s: Molar fraction in range from 0 to 5 mol/lg NaCl

    Returns:
        dynamic viscosity of brine [Pa*s]
    
    Note:
        source: Mao S. & Duan Z. 2009
    """
    rho_H2O = (cp.PropsSI("D", "T", t, "P", p, "Water")) / 1000
    p = p / 1e5
    m = (s * 1000) / 58.44
    d1 = 0.28853170e7
    d2 = -0.11072577e5
    d3 = -0.90834095e1
    d4 = 0.30925651e-1
    d5 = -1 * 0.27407100e-4
    d6 = -0.19283851e7
    d7 = 0.56216046e4
    d8 = 0.13827250e2
    d9 = -0.47609523e-1
    d10 = 0.35545041e-4
    a0 = -0.21319213
    a1 = 0.13651589e-2
    a2 = -0.12191756e-5
    b0 = 0.69161945e-1
    b1 = -0.27292263e-3
    b2 = 0.20852448e-6
    c0 = -0.25988855e-2
    c1 = 0.77989227e-5
    d = np.array([d1, d2, d3, d4, d5, d6, d7, d8, d9, d10])
    a = a0 + a1 * t + a2 * (t**2)
    b = b0 + b1 * t + b2 * (t**2)
    c = c0 + c1 * t
    ln_nr = a * m + b * (m**2) + c * (m**3)
    nr = np.exp(ln_nr)
    sum1 = 0
    for i in range(1, 6):
        term = d[i - 1] * t ** (i - 3)
        sum1 += term

    sum2 = 0
    for i in range(6, 11):
        term = d[i - 1] * rho_H2O * t ** (i - 8)
        sum2 += term
    ln_eta_H2O = sum1 + sum2
    n_H2O = np.exp(ln_eta_H2O)
    n_sol = n_H2O * nr
    return n_sol


def brine_compressibility(t, p, s):
    """Brine compressibility

    Args:
        t: temperature 
        p: pressure 
        s: mineralisation
    
    Returns:
        brine compressibility  [1/Pa]

    """
    out = np.abs(
        15.5
        / (0.485 * p / 1.0e5 + 0.5415 * (s * 1000.0) - 966.6 * t - 156.28e3)
    )
    return out


def compressibility(tr, pr, s):
    """Compressibility

    Args:
        t: temperature 
        p: pressure 
        s: mineralisation

    Returns:
        compressibility  [1/Pa]

    """
    p = pr / 6894.757
    t = (9.0 / 5.0) * (tr - 273.15) + 32.0
    dVwt = -1.0001e-2 + 1.33391e-4 * t + 5.50654e-7 * t**2.0
    dVwp = (
        -1.95301e-9 * p * t
        - 1.72834e-13 * p**2.0 * t
        - 3.58922e-7 * p
        - 2.25341e-10 * p**2.0
    )
    b = (1 - dVwp) * (1 + dVwt)
    return b


def heat_convection(t, p, s, d, w, lch, t_wall):
    """Heat convection

    Args:
        t: temperature 
        p: pressure 
        s: mineralisation
        d: diameter
        lch: characteristic length
        t_wall: temperature of the wall of well

    Returns:
        heat convection  [W/(m^2*K)]

    """
    g = 9.81
    Re = (density(t, p, s) * w * d) / dyn_viscosity(t, p, s)
    Pr = c_heat(t, p, s) * dyn_viscosity(t, p, s) / heat_conduction(t, p, s)

    if Re < 2000.0:
        Pr_wsf = (
            c_heat(t_wall, p, s)
            * dyn_viscosity(t_wall, p, s)
            / heat_conduction(t_wall, p, s)
        )
        Gr = (
            g
            * lch**3.0
            * brine_compressibility(t, p, s)
            * (np.abs(t_wall - t))
        ) / (dyn_viscosity(t, p, s) / density(t, p, s)) ** 2.0
        alfa = (
            (0.15 * Re**0.33 * Gr**0.1 * Pr**0.43 * (Pr / Pr_wsf) ** 0.25)
            * heat_conduction(t, p, s)
            / d
        )
    else:
        alfa = (0.0155 * Re**0.83 * Pr**0.5) * (heat_conduction(t, p, s) / d)

    return alfa


def heat_conduction(t, pr, s):
    """Heat conduction

    Args:
        t: temperature 
        p: pressure 
        s: mineralisation

    Returns:
        heat conduction  [W/(m*K)]
    """
    p = pr / 1e6
    k_H20 = (
        ((7e-9) * (t**3))
        - ((1.5113e-5) * (t**2))
        + ((8.801e-3) * t)
        - (0.8624)
        + ((1.57e-6) * p * t)
    )
    salt_molecular_weight = 58.44
    water_molecular_weight = 18.015
    # Na_molecular_weight = 22.99
    # Cl_molecular_weight = 35.34
    salt_moles = (s * 1000) / salt_molecular_weight
    water_moles = 1000 / water_molecular_weight
    X_Na = salt_moles / (salt_moles + water_moles)
    X_Cl = salt_moles / (salt_moles + water_moles)
    a1Na = 0
    a2Na = 0
    a_Na = a1Na + a2Na * math.exp(-0.023 * (t - 273.15))
    a1Cl = -0.360439
    a2Cl = 0.006076
    a_Cl = a1Cl + a2Cl * math.exp(-0.023 * (t - 273.15))
    k_ions_solvent = a_Na * X_Na + a_Cl * X_Cl
    k = k_H20 + k_ions_solvent
    return k


def dp_flow(t, p, s, fi, w, l, d):
    """Flow pressure drops

    Args:
        t: temperature 
        p: pressure 
        s: mineralisation
        fi: coefficient of pipe roughness
        w: velocity of brine
        l: characteristic length
        d: diameter

    Returns
        pressure drop  [Pa/m]

    """

    Re = (density(t, p, s) * w * d) / dyn_viscosity(t, p, s)

    if Re == 0:
        beta = 0
    elif (Re > 0) and (Re < 3.0e3):
        beta = 64.0 / Re
    elif (Re >= 3.0e3) and (Re <= 1.0e4):
        beta = 0.3164 / Re**0.25
    elif Re > 1.0e4:
        beta = (0.221 / Re**0.237) + 0.0032

    dp = (fi * beta * (w**2.0) * density(t, p, s) * l) / (2.0 * d)

    return dp
