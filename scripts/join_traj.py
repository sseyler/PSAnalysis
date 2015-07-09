import MDAnalysis as mda

u = mda.Universe("open_final.psf", "adk-tmd_1.dcd")

w = mda.Writer("path.dcd", u.trajectory.numatoms)

for i in xrange(99):
    trajname = "adk-tmd_%i.dcd" % (i+1)
    u = mda.Universe("open_final.psf", trajname)
    for ts in u.trajectory:
        w.write(ts)