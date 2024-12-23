import openmc
import math
import sys
sys.path.append('../')
from constants import core_height, rod_pitch, turnkey_size, ring_number
from constants import r_fuel, r_hole, r_fr, delta_shell, fr_number, n_fa, cr_pos

sys.path.append('../materials/')
from materials import shell_110, gaz, coolant_arr, fuel_uo2_arr, grey_fuel_arr, s_ring
from materials import shell_temperature, gaz_gap_temperature, central_gaz_temperature
from fuel_assemblies import fa_types, find_name

def make_fuel_rod(num, j, fr_num, r_hole, r_fuel, delta_shell, r_fr, helium, fuel, clad, water, central_temp, gaz_gap_temp, shell_temp):
    hole_c = openmc.ZCylinder(surface_id=int(2E11 + num * 1E9 + j*1E6 + fr_num * 1E3 + 99), r=r_hole)
    fuel_c = openmc.ZCylinder(surface_id=int(2E11 + (num + 1) * 1E9 + j*1E6 + fr_num * 1E3 + 99), r=r_fuel)
    clad_c = openmc.ZCylinder(surface_id=int(2E11 + (num + 2) * 1E9 + j*1E6 + fr_num * 1E3 + 99), r=r_fr)
    inner_clad_c = openmc.ZCylinder(surface_id=int(2E11 + (num + 3) * 1E9 + j*1E6 + fr_num * 1E3 + 99), r=r_fr - delta_shell)
    central_hole_cell = openmc.Cell(cell_id=int(3E11 + num * 1E9 + j*1E6 + fr_num * 1E3 + 99), fill=helium, region=-hole_c)
    central_hole_cell.temperature = central_temp
    fuel_cell = openmc.Cell(cell_id=int(3E11 + (num + 1) * 1E9 + j*1E6 + fr_num * 1E3 + 99), fill=fuel, region=+hole_c & -fuel_c)
    gap_cell = openmc.Cell(cell_id=int(3E11 + (num + 2) * 1E9 + j*1E6 + fr_num * 1E3 + 99), fill=helium, region=+fuel_c & -inner_clad_c)
    gap_cell.temperature = gaz_gap_temp
    clad_cell = openmc.Cell(cell_id=int(3E11 + (num + 3) * 1E9 + j*1E6 + fr_num * 1E3 + 99), fill=clad, region=+inner_clad_c & -clad_c)
    clad_cell.temperature = shell_temp
    water_cell = openmc.Cell(cell_id=int(3E11 + (num + 4) * 1E9 + j*1E6 + fr_num * 1E3 + 99), fill=water, region=+clad_c)
    pin_universe = openmc.Universe(universe_id = int(4E11 + num * 1E9 + j*1E6 + fr_num * 1E3 + 99),
                                   cells=[central_hole_cell, fuel_cell, gap_cell, clad_cell, water_cell],
                                   name=f"fa_number_{j}__fr_number_{fr_num}")
    return pin_universe

def grey_fuel_rod(num, j, fr_num, r_hole, delta_shell, r_fr, helium, grey_fuel, clad, water, central_temp, gaz_gap_temp, shell_temp):
    hole_c = openmc.ZCylinder(surface_id=int(2E11 + num * 1E9 + j*1E6 + fr_num * 1E3 + 99), r=r_hole)
    clad_c = openmc.ZCylinder(surface_id=int(2E11 + (num + 1) * 1E9 + j*1E6 + fr_num * 1E3 + 99), r=r_fr)
    inner_clad_c = openmc.ZCylinder(surface_id=int(2E11 + (num + 2) * 1E9 + j*1E6 + fr_num * 1E3 + 99), r=r_fr - delta_shell)
    r_rings = []
    ring_cells = []
    r_rings.append(hole_c)
    for i in range(0, ring_number):
        r_ring = math.sqrt(s_ring/math.pi + r_rings[i].r * r_rings[i].r)
        r_rings.append(openmc.ZCylinder(surface_id=int(2E11 + (num + 4) * 1E9 + j*1E6 + fr_num * 1E3 + i), r=r_ring))
        ring_cells.append(openmc.Cell(cell_id=int(3E11 + num * 1E9 + j*1E6 + fr_num * 1E3 + i), fill=grey_fuel[i], region= -r_rings[i+1] & +r_rings[i]))
    central_hole_cell = openmc.Cell(cell_id=int(3E11 + (num + 1)* 1E9 + j*1E6 + fr_num * 1E3 + 99), fill=helium, region=-hole_c)
    central_hole_cell.temperature = central_temp
    gap_cell = openmc.Cell(cell_id=int(3E11 + (num + 2) * 1E9 + j*1E6 + fr_num * 1E3 + 99), fill=helium, region=+r_rings[-1] & -inner_clad_c)
    gap_cell.temperature = gaz_gap_temp
    clad_cell = openmc.Cell(cell_id=int(3E11 + (num + 3) * 1E9 + j*1E6 + fr_num * 1E3 + 99), fill=clad, region=+inner_clad_c & -clad_c)
    clad_cell.temperature = shell_temp
    water_cell = openmc.Cell(cell_id=int(3E11 + (num + 4) * 1E9 + j*1E6 + fr_num * 1E3 + 99), fill=water, region=+clad_c)
    grey_pin_universe = openmc.Universe(universe_id = int(4E11 + num * 1E9 + j*1E6 + fr_num * 1E3 + 99),
                                   cells=[central_hole_cell, gap_cell, clad_cell, water_cell] + ring_cells,
                                   name=f"grey_fa_number_{j}__fr_number_{fr_num}")
    return grey_pin_universe

def fuel_assembly(i, shell, gaz, coolant, uo2_fa_arr, grey_fa_arr,
                  shell_temp, gaz_gap_temp, central_gaz_temp,
                  r_hole, r_fuel, delta_shell, r_fr, type,
                  fr_number, rod_pitch, cr_pos, turnkey_size, boundary):
    #central tube
    central_tube_out_cylinder = openmc.ZCylinder(surface_id=int(2E11 + 1 * 1E9 + i*1E6 + 312 * 1E3 + 99), r=0.65)
    central_tube_in_cylinder = openmc.ZCylinder(surface_id=int(2E11 + 2 * 1E9 + i*1E6 + 321* 1E3 + 99), r=0.55)
    #control sistem channel
    csc_out_cylinder = openmc.ZCylinder(surface_id=int(2E11 + 3 * 1E9 + i*1E6 + 312* 1E3 + 99), r=0.63)
    csc_in_cylinder = openmc.ZCylinder(surface_id=int(2E11 + 4 * 1E9 + i*1E6 + 312* 1E3 + 99), r=0.545)

    #create central tube
    central_tube_in_cell = openmc.Cell(cell_id = int(3E11 + 1 * 1E9 + i*1E6 + 312* 1E3 + 99),
                                       fill = coolant, region = -central_tube_in_cylinder)
    central_tube_cell = openmc.Cell(cell_id = int(3E11 + 2 * 1E9 + i*1E6 + 312* 1E3 + 99),
                                    fill = shell, region = +central_tube_in_cylinder & -central_tube_out_cylinder)
    central_tube_cell.temperature = coolant.temperature
    water2_cell = openmc.Cell(cell_id = int(3E11 + 3 * 1E9 + i*1E6 + 312* 1E3 + 99),
                              fill = coolant, region = +central_tube_out_cylinder)
    ct = openmc.Universe(universe_id = int(4E11 + 1 * 1E9 + i*1E6 + 312* 1E3 + 99),
                         cells=[central_tube_in_cell, central_tube_cell, water2_cell])

    #create control sistem channel
    csc_cell = openmc.Cell(cell_id = int(3E11 + 4 * 1E9 + i*1E6 + 312* 1E3 + 99),
                           fill = shell, region = +csc_in_cylinder & -csc_out_cylinder)
    csc_cell.temperature = coolant.temperature
    water3_cell = openmc.Cell(cell_id = int(3E11 + 5 * 1E9 + i*1E6 + 312* 1E3 + 99),
                              fill = coolant, region = +csc_out_cylinder)
    csc_in_cell = openmc.Cell(cell_id = int(3E11 + 6 * 1E9 + i*1E6 + 312* 1E3 + 99),
                              fill = coolant, region = -csc_in_cylinder)
    csc = openmc.Universe(universe_id = int(4E11 + 2* 1E9 + i*1E6 + 312* 1E3 + 99),
                          cells = [csc_in_cell, csc_cell, water3_cell])

    rods_arr = []
    n = 0
    n_ = 0
    for j in range(0, fr_number + 18):
        if j == cr_pos[n_]:
            rods_arr.append(csc)
            n_ += 1
        else:
            if len(type["grey_pos"]) == 0:
                fr = make_fuel_rod(6, i, j, r_hole, r_fuel, delta_shell, r_fr,
                                    gaz, uo2_fa_arr[j], shell, coolant,
                                    central_gaz_temp, gaz_gap_temp, shell_temp)
            else:
                if j == type["grey_pos"][n]:
                    fr = grey_fuel_rod(11, i, j, r_hole, delta_shell, r_fr, gaz, grey_fa_arr[n], shell, coolant,
                                        central_gaz_temp, gaz_gap_temp, shell_temp)
                    n += 1
                else:
                    fr = make_fuel_rod(6, i, j, r_hole, r_fuel, delta_shell, r_fr,
                                    gaz, uo2_fa_arr[j], shell, coolant,
                                    central_gaz_temp, gaz_gap_temp, shell_temp)
            rods_arr.append(fr)

    all_water_out=openmc.Cell(cell_id = int(3E11 + 16 * 1E9 + i*1E6 + 312* 1E3 + 99), fill=coolant)
    all_water_out_u=openmc.Universe(universe_id = int(4E11 + 16 * 1E9 + i*1E6 + 312* 1E3 + 99), cells=[all_water_out])

    lat=openmc.HexLattice(lattice_id = int(5E11 + 1E9 + i*1E6 + 312* 1E3 + 99), name=f'assembly_{i}')
    lat.center = (0.0, 0.0)
    lat.pitch = (rod_pitch,)
    lat.outer=all_water_out_u
    lat.orientation = 'x'

    ring10 = rods_arr[0:60]
    ring9 = rods_arr[60:114]
    ring8 = rods_arr[114:162]
    ring7 = rods_arr[162:204]
    ring6 = rods_arr[204:240]
    ring5 = rods_arr[240:270]
    ring4 = rods_arr[270:294]
    ring3 = rods_arr[294:312]
    ring2 = rods_arr[312:324]
    ring1 = rods_arr[324:330]

    ring0 = [ct]
    rings = [ring10, ring9, ring8, ring7, ring6, ring5, ring4, ring3, ring2, ring1, ring0]
    lat.universes = rings
    assembly_cell = openmc.Cell(cell_id = int(3E11 + 17 * 1E9 + i*1E6 + 312* 1E3 + 99), name=f'cell_assembly_{i}')
    hex_prizm = openmc.model.HexagonalPrism(edge_length = turnkey_size/math.sqrt(3), orientation = 'x', boundary_type = boundary)
    assembly_cell.region = -hex_prizm
    assembly_cell.fill = lat

    fa_un = openmc.Universe(universe_id = int(4E11 + 17 * 1E9 + i*1E6 + 312* 1E3 + 99), cells=[assembly_cell])

    return fa_un
