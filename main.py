"""This is main programs's executable file

    Note:
        Authors: L.Pająk, K.Pierzchała, M.Miecznik\n
        Affiliation: Mineral and Energy Economy Research Institute\n
        Polish Academy of Sciences (MEERI PAS), Kraków, Poland\n
        Date: March 2025\n
        Version: 1.0\n
        Welbore flow simulator (More information in the user's guide)\n
        Calculator developed in the GeoModel project (https://geomodel.pl/en/)\n

    Important:
        Change argument for 1 to modell production well, 0 for injection well in first if condition i main module.
        All crucial functions are exceuted here. They are described in detail in particular modules.       
"""
import time
import pandas as pd
import major_lib as ml


def main():
    start_time = time.perf_counter()

    if 1:
        file_name = "well_prod_input.xlsx"
    else:
        file_name = "well_inj_input.xlsx"

    print("The name of the input file: ", file_name)

    well_geom = ml.well_geometry_in(file_name)
    wws = ml.working_scheme(file_name)
    others = ml.read_others(file_name)
    in_temp_lack = ml.prof_temp_start(file_name, well_geom)
    grid = ml.grid_gener(in_temp_lack, float(others[3, 0]), 1.25)
    r = grid[0]
    z = grid[1]
    nc = ml.nodes_coordinates(r, z)
    t0 = ml.initial_temp(nc, in_temp_lack)
    n2e = ml.nodes_to_elements(nc, r, z)
    nodes_in_elements = ml.nodes_in_elements(n2e)
    am = ml.a_mat_fix(file_name)
    mat = ml.check_mat(in_temp_lack, grid, n2e, nc)
    nodes_first_boundary = ml.first_type_boundary(nc, n2e, mat)
    active_nodes = ml.seek_actv_nodes(nc, nodes_first_boundary)
    nodes_bottom = ml.nodes_bottom(nc, z, r)
    nodes_top = ml.nodes_top(nc, z, r)
    nodes_in_well = ml.nodes_zones(z, nc, n2e, mat)
    nodes_with_brine = ml.nodes_filled_with_brine(nodes_in_well, n2e, mat)
    dtime_max = ml.max_dtime(nc, n2e, mat, am, 4.0)
    dtime = 2.0
    t0 = ml.eq_temp_in_zone(t0, nodes_in_well)
    twell0 = ml.temp_in_well(nodes_in_well, t0, nc)

    print("The suggested length of time step is : ", round(dtime, 2), " sec")

    if wws[0, 1] >= 0.0:
        result = ml.prod_calculations(
            wws,
            dtime,
            dtime_max,
            t0,
            twell0,
            n2e,
            nc,
            z,
            r,
            mat,
            am,
            nodes_in_well,
            nodes_with_brine,
            active_nodes,
            others,
            nodes_top,
            nodes_in_elements,
        )
    else:
        result = ml.inj_calculations(
            wws,
            dtime,
            dtime_max,
            t0,
            twell0,
            n2e,
            nc,
            z,
            r,
            mat,
            am,
            nodes_in_well,
            nodes_with_brine,
            active_nodes,
            others,
            nodes_top,
            nodes_bottom,
            nodes_in_elements,
        )

    if True:
        out = []
        for i, t in enumerate(t0):
            hlp = [nc[i, 1], nc[i, 2], t, result[i]]
            out.append(hlp)
        df = pd.DataFrame(out)
        head = ["Radius r [m]", "Depth z [m]", "T_init [°C]", "T_end [°C]"]
        df.to_excel("out.xlsx", index=False, header=head)

    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print(f"Simulatio duration: {elapsed_time:.4f} seconds")


main()
