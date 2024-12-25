import openmc
import openmc.deplete
import os
import sys
import numpy as np

from depletion_steps import depletion_steps

sys.path.append('../../')
from constants import batches, inactive, particles, q_r, n_fa, fr_number, core_height

if __name__ == '__main__':
    chain_file = "/home/ubuntu24/Desktop/libs/chain_endfb71_pwr.xml"
    settings = openmc.Settings()
    settings.batches = batches
    settings.inactive = inactive
    settings.particles = particles

    settings.temperature = {'method': 'interpolation'}
    settings.export_to_xml()

    chain = openmc.deplete.Chain.from_xml(chain_file)

    geometry = os.path.join('..', 'geometry.xml')
    materials = os.path.join('..', 'materials.xml')
    model = openmc.model.Model.from_xml(geometry = geometry, materials = materials, settings='settings.xml')

    op = openmc.deplete.CoupledOperator(model, chain_file)

    sec_in_day = 24 * 60 * 60
    #depletion_steps = [0.01, 0.01]
    #timesteps = [ds * sec_in_day for ds in depletion_steps]

    timesteps = depletion_steps

    powr = q_r/(n_fa * core_height*100) * 1.31

    integrator = openmc.deplete.CELIIntegrator(operator = op, timesteps = timesteps, power = powr, timestep_units='MWd/Kg')
    integrator.integrate()

    results = openmc.deplete.Results("depletion_results.h5")
    time, keff = results.get_keff()
    time = time / (sec_in_day)
    k_eff = []
    st_dev = []
    for i in range(0, len(time)):
        k_eff.append(keff[i][0])
        st_dev.append(keff[i][1])
    res = np.column_stack((time, np.array(k_eff), np.array(st_dev)))
    np.savetxt('results.txt', res, delimiter = '\t', fmt = '%.7f' )
    res = np.column_stack((time, [0.0] + depletion_steps, np.array(k_eff), np.array(st_dev)))
    np.savetxt('results.txt', res, delimiter = '\t', fmt = '%.7f' )
