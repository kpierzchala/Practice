"""Modules defines main functions connected to wells"""
import math
import numpy as np
import pandas as pd
import brine_prop as bp


def third_type_boundary(
    nlayer,
    alfa,
    tbrine,
    tnext,
    tside0,
    dtime,
    n2e,
    twell,
    nc,
    z,
    mat,
    am,
    nodes_in_elements,
):
    """Collecting all the data set neccessary for temperature of rock in a side of the well estimation
    
    Args:
        nlayer: number of layer in a well
        alfa: heat convection coefficient [W/(m2 K)]
        tbrine: the brine (fluid) temperature [°C]
        tnext: temperature of the next node, linked with an element belonging to a well side, next by radius
        dtime: time step [s]
        n2c: nodes in element
        twell: matrix includes information in rows - number of zones in a well, columns: 0-node number, 1-radius (coordinate R),2-hight counted from the bottom of the well (Z coordinate, 0. = bottom), 3-temperature distribution
        nc: matrix nodes coordinate
        z: coordinates of nodes by Z
        mat: matrix includes information how elements are filled by materials
        am: matrix of temperature equalization factor [m2/s]
        nodes_in_elements: matrix includes information how nodes are linked with elements, in rows - number of elements, in columns - number of nodes in an element

    Returns:
        time on a well side after the time step dtime

    """

    if 1:
        nw = int(twell[nlayer, 0])
        r = nc[nw, 1]
        dr = nc[nw + len(z), 1] - nc[nw, 1]

        if (nodes_in_elements[nw, 2] >= 0) and (nodes_in_elements[nw, 3] >= 0):
            x1 = nodes_in_elements[nw, 2]
            x2 = nodes_in_elements[nw, 3]
            a = (
                (nc[nw + 1, 2] - nc[nw, 2]) * am[mat[x1], 1]
                + (nc[nw, 2] - nc[nw - 1, 2]) * am[mat[x2], 1]
            ) / (nc[nw + 1, 2] - nc[nw - 1, 2])
            lbd = (
                (nc[nw + 1, 2] - nc[nw, 2]) * am[mat[x1], 2]
                + (nc[nw, 2] - nc[nw - 1, 2]) * am[mat[x2], 2]
            ) / (nc[nw + 1, 2] - nc[nw - 1, 2])

        if (nodes_in_elements[nw, 2] >= 0) and (nodes_in_elements[nw, 3] < 0):
            x = nodes_in_elements[nw, 2]
            a = am[mat[x], 1]
            lbd = am[mat[x], 2]

        if (nodes_in_elements[nw, 2] < 0) and (nodes_in_elements[nw, 3] >= 0):
            x = nodes_in_elements[nw, 3]
            a = am[mat[x], 1]
            lbd = am[mat[x], 2]

        m = dr**2.0 / (a * dtime)
        c = (alfa * dr) / lbd
        tsn = (1 / m) * (
            2.0 * c * tbrine
            + (m - 2.0 * (c + 1.0 - dr / (2.0 * r))) * tside0
            + 2.0 * tnext * (1.0 - dr / (2.0 * r))
        )
    else:
        tsn = tbrine

    return tsn


def prod_well(
    m,
    t,
    p,
    s,
    tr,
    trock_matrix,
    dtime,
    fi,
    acc,
    n2e,
    nc,
    z,
    mat,
    am,
    nodes_in_elements,
):
    """Description of production well

    Args:
        m: mass flux [kg/s]
        t: inlet temperature [°C]
        p: pressure at the depth of a well liner [Pa]
        S: TDS [kg/kg]
        tr: matrix of rock temperature in nodes belonged to a well (only), nw is a node number tr[n,0] - node number nw where is the side of a well, tr[n,1] - radius of a side of a well [m], tr[n,2] - Z coordinate [m], tr[n,3] - temperature [°C]
        trock_matrix: the matrix with temperature in nodes at the following time step [°C]. 1D matrix = vector
        dtime: time step value [s]
        fi: pipe wall roughness coefficient
        acc: required accuracy
        n2e: matrix describing how nodes belonges to elements
        nc: matrix with nodes coordinates R,Z
        z: matrix with the Z coordinate
        mat: the matrix materials linked to elements in a grid
        am: matrix of temperature equalization factor [m^2/s]
        nodes_in_elements: matrix includes information how nodes are linked with elements, in rows - number of elements, in columns - number of nodes in an element
        
    Returns:
        matrix with parameters at the begining and at the end of n zone. Col 0-rock temperature in nodes belonged to a well, 1-temperature of rock in a side of the well ,2-temperature of brine entering zone n, 3- temperature of brine leaving zone n, 4- pressure of brine entering zone n, 5-pressure of brine leaving zone n, 6- heat convection coefficient in zone n  
    
    Note:
        gravity value g=9.80655
    
    """
    g = 9.80655
    out = np.zeros([len(tr), 7], dtype=float)
    trock = np.zeros([len(tr)], dtype=float)
    trock2 = np.zeros([len(tr)], dtype=float)
    tbrine1 = np.zeros([len(tr)], dtype=float)
    tbrine2 = np.zeros([len(tr)], dtype=float)
    pbrine1 = np.zeros([len(tr)], dtype=float)
    pbrine2 = np.zeros([len(tr)], dtype=float)

    for i in range(len(tr), 0, -1):
        n = i - 1
        trock[n] = tr[n, 3]

        if n == len(tr) - 1:
            tbrine1[n] = t
            pbrine1[n] = p

        else:
            tbrine1[n] = tbrine2[n + 1]
            pbrine1[n] = pbrine2[n + 1]

        w = m / (
            bp.density(tbrine1[n] + 273.15, pbrine1[n], s)
            * 3.141592654
            * tr[n, 1] ** 2.0
        )
        pDynamic = (
            0.5 * bp.density(tbrine1[n] + 273.15, pbrine1[n], s) * w**2.0
        )
        roB = bp.density(tbrine1[n] + 273.15, pbrine1[n] - pDynamic, s)
        cB = bp.c_heat(tbrine1[n] + 273.15, pbrine1[n] - pDynamic, s)
        if n > 0:
            hzone = tr[n, 2] - tr[n - 1, 2]
        else:
            hzone = 0.0
        dp_flow_res = bp.dp_flow(
            tbrine1[n] + 273.15, pbrine1[n], s, fi, w, hzone, tr[n, 1] * 2.0
        )
        alfaB = bp.heat_convection(
            tbrine1[n] + 273.15,
            pbrine1[n] - pDynamic,
            s,
            2.0 * tr[n, 1],
            w,
            hzone,
            tr[n, 3] + 273.15,
        )
        A = alfaB * (math.pi * 2.0 * tr[n, 1] * hzone) / ((m + 1e-3) * cB)
        if hzone > 0.0:
            tbrine2[n] = (
                2.0 * tbrine1[n] + A * (tr[n, 3] + tr[n - 1, 3] - tbrine1[n])
            ) / (A + 2.0)
        else:
            tbrine2[n] = tbrine1[n]
        pbrine2[n] = (
            pbrine1[n] - (g * roB * (tr[n, 2] - tr[n - 1, 2])) - dp_flow_res
        )

        if n == 0:
            pbrine2[n] = pbrine1[n]

        if n == 1:
            print(
                f"tbrine1[n]:{n}, {round(tbrine1[n],3)}, {round(pbrine1[n]/1e5,3)}, {dp_flow_res}"
            )

        trock2[n] = third_type_boundary(
            n,
            alfaB,
            tbrine1[n],
            trock_matrix[int(tr[n, 0] + len(z))],
            trock[n],
            dtime,
            n2e,
            tr,
            nc,
            z,
            mat,
            am,
            nodes_in_elements,
        )

        if i == 1:
            PR = [
                tbrine1[0],
                pbrine1[0],
                w[0],
                pDynamic[0],
                roB[0],
                cB[0],
                dp_flow_res[0],
                m,
                alfaB[0],
            ]

        out[n, 0] = trock[n]
        if n == len(tr) - 1:
            out[n, 1] = t
        else:
            out[n, 1] = trock2[n]

        out[n, 2] = tbrine1[n]
        out[n, 3] = tbrine2[n]
        out[n, 4] = pbrine1[n]
        out[n, 5] = pbrine2[n]
        out[n, 6] = alfaB

    df = pd.DataFrame(out)
    df.to_excel(
        "prod_well_output.xlsx",
        index=False,
        header=[
            "trock",
            "trock2",
            "tbrine1",
            "tbrine2",
            "pbrine1",
            "pbrine2",
            "alfaB",
        ],
    )

    return [out, PR]


def inj_well(
    m,
    t,
    p,
    s,
    tr,
    trock_matrix,
    dtime,
    fi,
    acc,
    n2e,
    nc,
    z,
    mat,
    am,
    nodes_top,
    nodes_in_elements,
):
    """Description of injection well

    Args:
        m: mass flux [kg/s]
        t: inlet temperature [°C]
        p: pressure at the depth of a well liner [Pa]
        s: TDS [kg/kg]
        tr: matrix of rock temperature in nodes belonged to a well (only), nw is a node number tr[n,0] - node number nw where is the side of a well, tr[n,1] - radius of a side of a well [m], tr[n,2] - Z coordinate [m], tr[n,3] - temperature [°C]
        trock_matrix: the matrix with temperature in nodes at the following time step [°C]. 1D matrix = vector
        dtime: time step value [s]
        fi: pipe wall roughness coefficient, usually has the values:1 –for smooth, new pipes made of steel, copper, brass, 1.1 – for new cast iron pipes, 1.56 – for cast iron and riveted steel, 1.78 – for old steel pipes (this value can be considered typical for drill pipes), 2.22 – for old riveted steel pipes,
        acc: required accuracy
        n2e: matrix describing how nodes belonges to elements
        nc: matrix with nodes coordinates R,Z
        z: matrix with the Z coordinate
        mat: the matrix materials linked to elements in a grid
        am: matrix of temperature equalization factor [m2/s]
        nodes_top: nodes in
        nodes_in_elements: matrix includes information how nodes are linked with elements, in rows - number of elements, in columns - number of nodes in an element

    Returns:
        matrix with parameters at the begining and at the end of n zone. Col 0-rock temperature in nodes belonged to a well, 1-temperature of rock in a side of the well ,2-temperature of brine entering zone n, 3- temperature of brine leaving zone n, 4- pressure of brine entering zone n, 5-pressure of brine leaving zone n, 6- heat convection coefficient in zone n
        
    Note:
        gravity value g=9.80655
    
    """
    g = 9.80655
    out = np.zeros([len(tr), 7], dtype=float)
    trock = np.zeros([len(tr)], dtype=float)
    trock2 = np.zeros([len(tr)], dtype=float)
    tbrine1 = np.zeros([len(tr)], dtype=float)
    tbrine2 = np.zeros([len(tr)], dtype=float)
    pbrine1 = np.zeros([len(tr)], dtype=float)
    pbrine2 = np.zeros([len(tr)], dtype=float)

    for i in range(0, len(tr)):
        n = i

        trock[n] = tr[n, 3]

        if n == 0:
            tbrine1[n] = t
            pbrine1[n] = p
        else:
            tbrine1[n] = tbrine2[n - 1]
            pbrine1[n] = pbrine2[n - 1]

        w = m / (
            bp.density(tbrine1[n] + 273.15, pbrine1[n], s)
            * (math.pi * tr[n, 1] ** 2.0)
        )
        pDynamic = (
            0.5 * bp.density(tbrine1[n] + 273.15, pbrine1[n], s) * w**2.0
        )
        roB = bp.density(tbrine1[n] + 273.15, pbrine1[n] - pDynamic, s)
        cB = bp.c_heat(tbrine1[n] + 273.15, pbrine1[n] - pDynamic, s)
        if n > 0:
            hzone = tr[n, 2] - tr[n - 1, 2]
        else:
            hzone = 0.0
        dp_flow_res = bp.dp_flow(
            tbrine1[n] + 273.15, pbrine1[n], s, fi, w, hzone, tr[n, 1] * 2.0
        )
        alfaB = bp.heat_convection(
            tbrine1[n] + 273.15,
            pbrine1[n] - pDynamic,
            s,
            2.0 * tr[n, 1],
            w,
            hzone,
            tr[n, 3] + 273.15,
        )

        A = alfaB * (math.pi * 2.0 * tr[n, 1] * hzone) / (m * cB)

        if n < (len(tr) - 1):
            tbrine2[n] = (
                2.0 * tbrine1[n] + A * (tr[n, 3] + tr[n + 1, 3] - tbrine1[n])
            ) / (A + 2.0)
        else:
            tbrine2[n] = (
                2.0 * tbrine1[n] + A * (2.0 * tr[n, 3] - tbrine1[n])
            ) / (A + 2.0)

        if n == (len(tr) - 1):
            pbrine2[n] = pbrine1[n] + (g * roB * hzone) - dp_flow_res
        else:
            pbrine2[n] = pbrine1[n] + (g * roB * hzone) - dp_flow_res

        trock2[n] = third_type_boundary(
            n,
            alfaB,
            tbrine1[n],
            trock_matrix[int(tr[n, 0]) + len(z)],
            trock[n],
            dtime,
            n2e,
            tr,
            nc,
            z,
            mat,
            am,
            nodes_in_elements,
        )

        out[n, 0] = trock[n]

        if n == 0:
            out[n, 1] = t
        else:
            out[n, 1] = trock2[n]

        out[n, 2] = tbrine1[n]
        out[n, 3] = tbrine2[n]
        out[n, 4] = pbrine1[n]
        out[n, 5] = pbrine2[n]
        out[n, 6] = alfaB

    return out


def twell_convert_twell0(twell0, z):
    """Function which convert initial parameters of well
    """
    X = np.zeros([len(twell0), 4])
    for nw in range(0, len(twell0)):
        X[nw, 0] = twell0[nw, 0]
        X[nw, 1] = twell0[nw, 1]
        X[nw, 2] = twell0[nw, 2]
        X[nw, 3] = z[nw - len(twell0), 1]
    return X
