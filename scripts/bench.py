import numpy as np
import pathgroup3 as pg

if __name__ == "__main__":
    ntrj = 100
    nproc = 11
    selection = 'all'
    outputname = 'dh_mat_ws_adkaa1-100gp.dat'

    workdir = '/nfs/homes/sseyler/Desktop/trj_RADICAL/'
    top = workdir + 'adk4ake.psf'
    top2 = workdir + '1ake_a_core.pdb'

    # trj_list = ['%sraw_trj/dims%04i_fit-core.dcd' % (workdir, i) for i in xrange(1,ntrj+1)]
    trj_list2 = ['%sraw_trj/pathway%i_fit-core.dcd' % (workdir, i) for i in xrange(1,ntrj+1)]

    # group = pg.PathGroup(trj_list,topology=top,selection=selection)
    group2 = pg.PathGroup(trj_list2,topology=top2,selection=selection)
    newgroup = group2 #+ group2

    print newgroup[0].getCoordinates().shape

    dh = newgroup.compareTo(metric='dh',processes=nproc)

    np.savetxt(outputname,dh,fmt='%5.3f')