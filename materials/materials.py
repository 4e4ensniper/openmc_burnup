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
from constants import core_height, csv_path
from constants import b_conc, dif_fu_cart
from constants import numbers, fr_number
from constants import r_fr, delta_shell, l_shell, r_fuel, l_g, l_fuel, r_hole

#materials ID definition
#1????? - first 1
#1__??? - material number
#1??___ - fuel asssembly number

def cr_steel(j, num, temp):
    _08x18h10t = openmc.Material(material_id = int(1E7 + num*1E5 + j), name = "08x18h10t")
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
    _08x18h10t.temperature = temp
    return _08x18h10t

def cr_helium_(j, num, temp):
    helium=openmc.Material(material_id = int(1E5 + num * 1E3 + j), name = "He")
    helium.add_element('He', 1.0)
    helium.temperature = temp
    helium.set_density('g/cm3', 3.24e-3)
    return helium

def cr_uo2_fuel(j, num, temp, enrich):
    fu = openmc.Material(material_id = int(1E7 + num * 1E5 + j * 1E2 + ring), name = f"UO2_{enrich}")
    fu.add_element('U', 1.0, enrichment = enrich, enrichment_type='wo')
    fu.add_element('O', 2.0)
    fu.set_density('g/cm3', 10.48)
    fu.temperature = temp
    return fu

def cr_uo2_gdo2(j, num, ring, temp, enrich, gdo2_pt):
    uo2 = cr_uo2_fuel(j, num, temp, enrich)
    gdo2 = openmc.Material(material_id = int(1E7 + (num + 1) * 1E5 + j*1E2 + ring), name = 'GdO2')
    gdo2.add_element('Gd', 2.0)
    gdo2.add_element('O', 3.0)
    gdo2.set_density('g/cm3', 7.407)
    gdo2_uo2 = openmc.Material.mix_materials([uo2, gdo2], [1-gdo2_pt*1E-2, gdo2_pt*1E-2], 'wo')
    gdo2_uo2.id = int(1E7 + (num + 2) * 1E5 + j * 1E2 + ring)
    gdo2_uo2.name = 'GdO2_UO2'
    gdo2_uo2.temperature = temp
    return gdo2_uo2

def cr_water(j, num, temp, density, b_conc):
    if b_conc > 0.00001:
        b_ppm = 1/(1 + 61.83/18 * (1/(b_conc*1E-3)-1)) * 1E6
        water = openmc.model.borated_water(boron_ppm = b_ppm, density=density*1E-3)
        water.id = int(1E5 + num * 1E3 + j)
        water.temperature = temp + 273.15
        water.name = 'H2O_b'
    else:
        water = openmc.Material(material_id = int(1E5 + (num+1) * 1E3 + j), name = "H2O")
        water.add_element('H', 2.0)
        water.add_element('O', 1.0)
        water.set_density('g/cm3', density*1E-3)
        water.temperature = temp + 273.15
        water.add_s_alpha_beta('c_H_in_H2O')
    return water

file_path = 'burnup_temp.txt'
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

