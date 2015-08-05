from psa import *
from MDAnalysis import Universe
from MDAnalysis.analysis.align import rotation_matrix
import numpy as np
import os, sys
WORKDIR = '/nfs/homes/sseyler/Repositories/python/psanalysis'
# WORKDIR = '/Users/sseyler/Repositories/python/psanalysis'
sys.path.append(WORKDIR)



print("Generating AdK CORE C-alpha reference coords + structure...")
cref_filename = '%s/structs/1ake_a_ca_core.pdb' % WORKDIR
oref_filename = '%s/structs/4ake_a_ca_core.pdb' % WORKDIR

c_ref = MDAnalysis.Universe(cref_filename)
o_ref = MDAnalysis.Universe(oref_filename)
u_ref = MDAnalysis.Universe(cref_filename)

c_ref_ca = c_ref.selectAtoms('name CA')
o_ref_ca = o_ref.selectAtoms('name CA')

adkCORE_resids = "(resid 1:29 or resid 60:121 or resid 160:214)"
c_ref_CORE_ca = c_ref_ca.selectAtoms(adkCORE_resids).coordinates() \
        - c_ref_ca.selectAtoms(adkCORE_resids).centerOfMass()
o_ref_CORE_ca = o_ref_ca.selectAtoms(adkCORE_resids).coordinates() \
        - o_ref_ca.selectAtoms(adkCORE_resids).centerOfMass()
ref_coords = 0.5*(c_ref_CORE_ca + o_ref_CORE_ca)

u_ref.atoms.translate(-c_ref_ca.selectAtoms(adkCORE_resids).centerOfMass())
o_ref.atoms.translate(-o_ref_ca.selectAtoms(adkCORE_resids).centerOfMass())
u_ref.selectAtoms(adkCORE_resids).CA.set_positions(ref_coords)



print("Building collection of simulations...")
# List of method names (same as directory names)
# method_names = ['DIMS', 'FRODA', 'MAP']
method_names = ['DIMS', 'TMD', 'FRODA', 'MAP', 'iENM', 'MENM-SP', 'MENM-SD',       \
                'MDdMD', 'GOdMD', 'Morph', 'ANMP', 'LinInt']
labels = [] # Heat map labels
simulations = [] # List of simulation topology/trajectory filename pairs

# Build list of simulations, each represented by a pair of filenames
# ([topology filename], [trajectory filename]). Generate corresponding label
# list.
for method in method_names:
    # Note: DIMS uses the PSF topology format
    topname = 'top.psf' if ('DIMS' in method or 'TMD' in method) else 'top.pdb'
    pathname = 'path.dcd'
    method_dir = '{}/methods/{}'.format(WORKDIR, method)
    if method is not 'LinInt':
        nruns = 3 if 'TMD' not in method else 6
        if 'TMD' in method:
            for run in xrange(1, nruns+1): # 3 runs per method
                run_dir = '{}/{:03n}'.format(method_dir, run)
                topology = '{}/{}'.format(method_dir, topname)
                trajectory = '{}/{}'.format(run_dir, pathname)
                labels.append(method + '(' + str(run) + ')')
                simulations.append((topology, trajectory))
        else:
            for run in xrange(1, nruns+1): # 3 runs per method
                run_dir = '{}/{:03n}'.format(method_dir, run)
                topology = '{}/{}'.format(method_dir, topname)
                trajectory = '{}/{}'.format(run_dir, pathname)
                labels.append(method + '(' + str(run) + ')')
                simulations.append((topology, trajectory))
    else: # only one LinInt trajectory
        topology = '{}/{}'.format(method_dir, topname)
        trajectory = '{}/{}'.format(method_dir, pathname)
        labels.append(method)
        simulations.append((topology, trajectory))



# universes = [] # List of MDAnalysis Universes representing simulations
# for sim in simulations:
#     universes.append(Universe(*sim))


import parallelizer as parallel
p = parallel.Parallelizer(workers=1)
p.prun(simulations)