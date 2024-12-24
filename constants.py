core_height = 3.53 #m
n_fa = 151
r_fuel = 7.57/20 #sm
r_hole = 2.35/20 #sm
r_fr = 9.1/20 #sm
delta_shell = 0.69/10 #sm
l_shell = 21.2 #W/m*K
l_g  = 0.25 #W/m*K
l_fuel = 2.9 #W/m*K
fr_number = 312
rod_pitch = 1.275 #sm
turnkey_size = 23.4 #sm
core_barrel_in_r = 339.72/2 #sm
core_barrel_out_r = 349.88/2 #sm
b_conc = 4 #g/kg
ring_number = 7
#The number of fuel assemblies per row, starting from the bottom of the core map.
line = [4, 7, 10, 11, 12, 13, 12, 13, 12, 13, 12, 11, 10, 7, 4]
#Path to the location of files with temperature distributions.
temp_path = "/home/ubuntu24/Desktop/openmc_burnup/materials/"
#csv_path = "/home/adminsrv/projects/EfanovKS/openmc_ap1000/materials/temperature_distributions/"
#Calculation parameters
batches = 200
inactive = 10
particles = 100000
#Number of different fuel assemblies defined as different universes
dif_fu_cart = ['Z49A2']
'''
dif_fu_cart = ['Z49A2', 'Z49A2',
               'Z40', 'Z24', 'Z33Z2', 'Z24',
               'Z24', 'Z33Z9', 'Z13',
               'Z13', 'Z33Z2', 'Z24',
               'Z24', 'Z13',
               'Z33Z9', 'Z24',
               'Z13',
               'Z33Z9']
'''
#An array that contains the numbers of the various universes of fuel assemblies on the cartogram.
half_numbers = [1, 2, 2, 1,
                3, 4, 5, 6, 5, 4, 3,
                1, 4, 7, 8, 9, 9, 8, 7, 4, 1,
                2, 5, 8, 10, 11, 12, 11, 10, 8, 5, 2,
                2, 6, 9, 11, 13, 14, 14, 13, 11, 9, 6, 2,
                1, 5, 9, 12, 14, 15, 16, 15, 14, 12, 9, 5, 1,
                4, 8, 11, 14, 16, 17, 17, 16, 14, 11, 8, 4,
                3, 7, 10, 13, 15, 17]

r_half_numbers = list(reversed(half_numbers))
half_numbers.append(len(dif_fu_cart))
numbers = half_numbers + r_half_numbers

cr_pos = [207, 213, 219, 225, 231, 237,
          240, 245, 250, 255, 260, 265,
          295, 298, 301, 304, 307, 310]

