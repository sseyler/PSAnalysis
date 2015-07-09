import numpy as np
from itertools import izip

# def write_xyz(filename, atoms, coordinates, box, title='simulation'):
#     """
#       Write atoms and coordinates to XYZ file *filename*.

#       :Arguments:
#         *filename*
#             name of the output file
#         *atoms*
#             list of the N atom names
#     """
#     from itertools import izip
#     with open(filename, "w") as xyz:
#         for frame in xrange(coordinates.shape[0]):
#             snapshot = coordinates[frame,:,:]
#             xyz.write("%d\n" % len(atoms))
#             xyz.write("frame  %s\n" % (frame+1))
#             # xyz.write("box %g %g %g  %s\n" % (box[0], box[1], box[2], frame))
#             for atom, (x,y,z) in izip(atoms, snapshot):
#                 xyz.write("%8s   %10.5f %10.5f %10.5f\n" % (atom, x, y, z))

def write_frame(xyz, atoms, snapshot, frame):
    xyz.write("%d\n" % len(atoms))
    xyz.write("frame  %s\n" % (frame+1))
    for atom, (x,y,z) in izip(atoms, snapshot):
        xyz.write("%8s   %10.5f %10.5f %10.5f\n" % (atom, x, y, z))

def write_xyz(filename, atoms, coordinates, title='simulation'):
    """
      Write atoms and coordinates to XYZ file *filename*.

      :Arguments:
        *filename*
            name of the output file
        *atoms*
            list of the N atom names
    """
    with open(filename, "w") as xyz:
        # Iterate through each frame
        for frame in xrange(coordinates.shape[0]):
            snapshot = coordinates[frame,:,:]
            write_frame(xyz, atoms, snapshot, frame)

# def write_frame(filename, atoms, coordinates, time, box):

#     from itertools import izip
#     with open(filename, "a") as xyz:
#         xyz.write("%d\n" % len(atoms))
#         xyz.write("box %g %g %g frame %0.1f\n" % (box[0], box[1], box[2], time))
#         for atom, (x,y,z) in izip(atoms, coordinates):
#             (x,y,z) = minimumImage((x,y,z), box)
#             xyz.write("%8s   %10.5f %10.5f %10.5f\n" % (atom, x, y, z))


# def minimumImage(x, box):
#     """
#     Returns a modified list of coordinates after applying periodic
#     boundary conditions specified by box.

#     Variables:
#         x - an Nx3 numpy array of coordinates
#         box - a 1x3 numpy array specifying the box size that enforces
#                 periodic boundary conditions, centered on the origin
#     """
#     return x - box*np.round(x/box)