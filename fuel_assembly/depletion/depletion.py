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
    gramms_hm = op.heavy_metal

    sec_in_day = 24 * 60 * 60

    timesteps = np.diff(depletion_steps)
    powr = q_r/(n_fa * core_height*100) * 1.31

    '''
    MWd_per_kgHM_convert_to_days = 1*gramms_hm / powr * 1E3
    print("\"Maximum\" depletion step: {:5.3} [d] in 1MWd/kgHM".format(MWd_per_kgHM_convert_to_days))

    MWd_per_kgHM_dt=np.diff(timesteps)

    timesteps_unk=timesteps * MWd_per_kgHM_convert_to_days
    print(timesteps_unk)
    '''
    integrator = openmc.deplete.CELIIntegrator(operator = op, timesteps = timesteps, power = powr, timestep_units='MWd/Kg')
    integrator.integrate()

    results = openmc.deplete.Results("depletion_results.h5")
    time, keff = results.get_keff()
    time = time / sec_in_day

    k_eff = []
    st_dev = []
    for i in range(0, len(time)):
        k_eff.append(keff[i][0])
        st_dev.append(keff[i][1])
    res = np.column_stack((time, depletion_steps, np.array(k_eff), np.array(st_dev)))
    np.savetxt('results.txt', res, delimiter = '\t', fmt = '%.7f' )
