import numpy as np
from itertools import izip

import MDAnalysis
from MDAnalysis.core.log import ProgressMeter

import sys
import os.path


class PDBToTraj(object):

    def __init__(self, universe, filename=None, infix='_conv',
                 targetdir=os.path.curdir, format='dcd'):
        self.universe = universe
        self.frames = self.universe.trajectory
        self.targetdir = targetdir
        self.ext = '.' + format

        if filename is None:
            head, tail = os.path.split(self.frames.filename)
            basename, ext = os.path.splitext(tail)
            filename = basename
        self.top_new = os.path.join(self.targetdir, 'top' + infix + '.pdb')
        self.trj_new = os.path.join(self.targetdir, filename + infix + self.ext)


    def convert(self, gentop=True):
        natoms = self.frames.numatoms
        nframes = self.frames.numframes
        if gentop:
            self.universe.atoms.write(self.top_new) # write first frame as topology
        self.frames.rewind()
        writer = MDAnalysis.Writer(self.trj_new, natoms)

        percentage = ProgressMeter(nframes, interval=10,
                format="Fitted frame %(step)5d/%(numsteps)d  [%(percentage)5.1f%%]\r")

        for ts in self.frames:
            writer.write(self.universe.atoms)
            percentage.echo(ts.frame)
        writer.close_trajectory()

    def get_top_and_traj(self):
        return self.top_new, self.trj_new


class Coord2Traj(object):

    def __init__(self, topology, timeseries, atoms, filename=None,
                 formats=['pdb'], save_xyz=False):
        """
          Write *timeseries* of *atoms* to new trajectory *format*.

          :Arguments:
            *timeseries*
               :class:`numpy.ndarray` of the time series of coordinates
            *atoms*
               list of *N* atom (or collective variable) names
            *filename*
               string, name of trajectory output file excluding extension
            *formats*
               list of strings, file extensions to set trajectory formats
            *save_xyz*
               boolean, if ``True`` keep intermediate xyz trajectory; delete
               otherwise [``False``]
        """
        self.topology = topology
        self.timeseries = timeseries
        self.atoms = atoms
        self.filename = filename or 'new_traj'
        self.formats = formats
        self.save_xyz = save_xyz

        self.nframes = self.timeseries.shape[0]
        self.N = len(self.atoms) # number of atoms

        if len(self.atoms) != self.timeseries.shape[1]:
            raise ValueError

        self._write_xyz_traj(self.filename, self.atoms, self.timeseries)


    def _write_xyz_frame(self, xyz_writer, atoms, coords, frame):
        xyz_writer.write("{0:d}\n".format(len(atoms)))
        xyz_writer.write("frame  %s\n" % (frame + 1))
        for atom, (x,y,z) in izip(atoms, coords):
            xyz_writer.write("%8s   %10.5f %10.5f %10.5f\n" % (atom, x, y, z))


    def _write_xyz_traj(self, filename, atoms, timeseries):
        """
          Write *atoms* and *timeseries* to XYZ file *filename*.

          :Arguments:
            *filename*
               name of the output file
            *atoms*
               list of the N atom names
            *timeseries*
               :class:`numpy.ndarray` for the time series to write
        """
        if '.xyz' not in filename:
            filename = filename + '.xyz'
        with open(filename, "w") as xyz:
            for frame, coords in enumerate(timeseries):
                self._write_xyz_frame(xyz, atoms, coords, frame)


    def _generate_topology(self, filename, atoms, coords):
        filename = (filename or self.filename) + '_top'
        topname_xyz = filename + '.xyz'
        topname_pdb = filename + '.pdb'

        with open(topname_xyz, "w") as xyz_top:
            self._write_xyz_frame(xyz_top, atoms, coords, 0)

        u = MDAnalysis.Universe(topname_xyz, multiframe=True)
        u.atoms.write(topname_pdb)

        return topname_pdb


    def convert(self, formats=None):
        formats = formats or self.formats
        coords_top = self.timeseries[0]
        # topname_pdb = self._generate_topology(self.filename, self.atoms, coords_top)

        trajname_xyz = self.filename + '.xyz'
        trajname_pdb = self.filename + '.pdb'

        u_xyz = MDAnalysis.Universe(self.topology, trajname_xyz, multiframe=True)
        pdb_writer = MDAnalysis.Writer(trajname_pdb, self.N)
        for ts in u_xyz.trajectory:
            pdb_writer.write(u_xyz.atoms)
        pdb_writer.close_trajectory()


        trajnames = []
        for format in formats:
            sys.stderr.write('Writing to {} format...'.format(format))
            converter = PDBToTraj(u_xyz, self.filename, infix='', format=format)
            converter.convert(gentop=False)
            topname, trajname = converter.get_top_and_traj()
            trajnames.append(trajname)
        sys.stderr.write('Finished trajectory conversion.')