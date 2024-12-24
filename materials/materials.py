import csv
import numpy as np
import openmc
import sys
from math import log
from math import pi
from math import sqrt
from scipy.integrate import quad

import openmc.model
from fuel_assemblies import fa_types, find_name

sys.path.append('../')
from constants import temp_path
from constants import b_conc, dif_fu_cart
from constants import numbers, fr_number, ring_number
from constants import r_fuel, r_hole

#materials ID definition
#1???????? - first 1
#1__?????? - material number
#1??__???? - fuel asssembly number
#1????___? - fuel rod number or 999
#1???????_ - ring in grey rod number or 9

def cr_steel(num):
    _08x18h10t = openmc.Material(material_id = int(1E8 + num * 1E6 + 99*1E4 + 999 * 1E1 + 9), name = "08x18h10t")
    _08x18h10t.add_element('Fe', 67.665,'wo')
    _08x18h10t.add_element('Cr', 18.0,'wo')
    _08x18h10t.add_element('Ni', 10.0,'wo')
    _08x18h10t.add_element('Ti', 0.6,'wo')
    _08x18h10t.add_element('C', 0.08,'wo')
    _08x18h10t.add_element('Si', 0.8,'wo')
    _08x18h10t.add_element('Mn', 2.0,'wo')
    _08x18h10t.add_element('Mo', 0.5,'wo')
    _08x18h10t.add_element('S', 0.02,'wo')
    _08x18h10t.add_element('P', 0.035,'wo')
    _08x18h10t.add_element('Cu', 0.3,'wo')
    _08x18h10t.set_density('g/cm3', 7.85)
    #_08x18h10t.temperature = temp
    return _08x18h10t

def cr_helium(j, num):
    helium=openmc.Material(material_id = int(1E8 + num * 1E6 + j*1E4 + 999 * 1E1 + 9), name = "He")
    helium.add_element('He', 1.0)
    #helium.temperature = temp
    helium.set_density('g/cm3', 3.24e-3)
    return helium

def cr_uo2_fuel(j, fr_num, num, temp, enrich, ring):
    #fu = openmc.Material(material_id = int(1E11 + num * 1E9 + j * 1E6 + fr_num * 1E3 + ring), name = f"UO2_{enrich}")
    fu = openmc.Material(material_id = int(1E8 + num * 1E6 + j*1E4 + fr_num * 1E1 + ring), name = f"UO2_{enrich}")
    fu.add_element('U', 1.0, enrichment = enrich, enrichment_type='wo')
    fu.add_element('O', 2.0)
    fu.set_density('g/cm3', 10.48)
    fu.temperature = temp + 273.15
    fu.depletable = True
    fu.volume = pi * (r_fuel * r_fuel - r_hole * r_hole)
    return fu

def cr_uo2_gdo2(j, fr_num, num, ring, temp, enrich, gdo2_pt, square):
    uo2 = cr_uo2_fuel(j, fr_num, num, temp, enrich, ring)
    gdo2 = openmc.Material(material_id = int(1E8 + (num + 1) * 1E6 + j*1E4 + fr_num * 1E1 + ring), name = f'GdO2_{j}_{fr_num}_{ring}')
    #gdo2 = openmc.Material(material_id = int(1E6 + (num + 1) * 1E5 + j*1E4 + fr_num * 1E3 + ring), name = 'GdO2')
    gdo2.add_element('Gd', 2.0)
    gdo2.add_element('O', 3.0)
    gdo2.set_density('g/cm3', 7.407)
    gdo2_uo2 = openmc.Material.mix_materials([uo2, gdo2], [1-gdo2_pt*1E-2, gdo2_pt*1E-2], 'wo')
    gdo2_uo2.id = int(1E8 + (num + 2) * 1E6 + j*1E4 + fr_num * 1E1 + ring)
    #gdo2_uo2.id = int(1E4+ (num + 2) * 1E4 + j * 1E3 + fr_num *1E3 + ring)
    gdo2_uo2.name = 'GdO2_UO2'
    gdo2_uo2.temperature = temp + 273.15
    gdo2_uo2.depletable = True
    gdo2_uo2.volume = square
    return gdo2_uo2

def cr_water(j, num, temp, density, b_conc):
    if b_conc > 0.00001:
        b_ppm = 1/(1 + 61.83/18 * (1/(b_conc*1E-3)-1)) * 1E6
        water = openmc.model.borated_water(boron_ppm = b_ppm, density=density*1E-3)
        water.id = int(1E8 + num * 1E6 + j*1E4 + 999 * 1E1 + 9)
        water.temperature = temp + 273.15
        water.name = 'H2O_b'
    else:
        water = openmc.Material(material_id = int(1E8 + (num + 1) * 1E6 + j*1E4 + 999 * 1E1 + 9), name = "H2O")
        water.add_element('H', 2.0)
        water.add_element('O', 1.0)
        water.set_density('g/cm3', density*1E-3)
        water.temperature = temp + 273.15
        water.add_s_alpha_beta('c_H_in_H2O')
    return water
def cr_shell_110(j, num):
    shell_alloy = openmc.Material(material_id = int(1E8 + num * 1E6 + j*1E4 + 999 * 1E1 + 9), name = "110")
    shell_alloy.add_element('Zr', 0.99, percent_type='wo')
    shell_alloy.add_element('Nb',0.1, percent_type='wo')
    #shell_alloy.temperature = temp + 273.15
    shell_alloy.set_density('g/cm3',6.5)
    return shell_alloy
file_path = temp_path + 'burnup_temp.txt'
hc_density = []
hc_temperature = []
shell_temperature = []
gaz_gap_temperature = []
fuel_temperature = []
central_gaz_temperature = []
with open(file_path, 'r') as file:
    for line in file:
        numbers = list(map(float, line.strip().split('\t')))
        hc_density.append(numbers[0])
        hc_temperature.append(numbers[1])
        shell_temperature.append(numbers[2])
        gaz_gap_temperature.append(numbers[3])
        fuel_temperature.append(numbers[4])
        central_gaz_temperature.append(numbers[5])
shell_110 = cr_shell_110(18, 1)
gaz = cr_helium(18, 2)
coolant_arr = []
fuel_uo2_arr = []
grey_fuel_arr = []
s_ring = pi * (r_fuel * r_fuel - r_hole * r_hole)/ring_number
for i in range(0, len(dif_fu_cart)):
    fuel_fa = []
    grey_fuel_fa = []
    type = find_name(dif_fu_cart[i], fa_types)
    coolant = cr_water(i, 3, hc_temperature[i], hc_density[i], b_conc)
    coolant_arr.append(coolant)
    for j in range (0, fr_number - len(type["grey_pos"])):
        fuel_uo2 = cr_uo2_fuel(i, j, 5, fuel_temperature[i], type["enrichment"], ring_number + 1)
        fuel_fa.append(fuel_uo2)
    fuel_uo2_arr.append(fuel_fa)
    for j in range (0, len(type["grey_pos"])):
        grey_rod = []
        for k in range(0, ring_number):
            grey_fuel = cr_uo2_gdo2(i, j, 6, k, fuel_temperature[i], type["grey_enrichment"], type["gdo2_wo"], s_ring)
            grey_rod.append(grey_fuel)
        grey_fuel_fa.append(grey_rod)
    grey_fuel_arr.append(grey_fuel_fa)
print("materials created!")
