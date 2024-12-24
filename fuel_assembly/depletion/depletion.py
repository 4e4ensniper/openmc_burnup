import openmc
import openmc.deplete
import os
import sys


sys.path.append('../../')
from constants import batches, inactive, particles

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

    op = openmc.deplete.Operator(model, chain_file)

    sec_in_day = 24 * 60 * 60
    depletion_steps = [0.01, 0.01, 0.01, 0.05]
    timesteps = [ds * sec_in_day for ds in depletion_steps]

    integrator = openmc.deplete.CELIIntegrator(op, timesteps, power_density=52.4)
    integrator.integrate()
