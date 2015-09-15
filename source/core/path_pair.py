# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see KATHRYN for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

def hausdorffNN(P, Q, N=None, indices=False):
    r"""Calculate the Hausdorff distance between two paths.

    .. |3Dp| replace:: :math:`N_p \times N \times 3`
    .. |2Dp| replace:: :math:`N_p \times (3N)`
    .. |3Dq| replace:: :math:`N_q \times N \times 3`
    .. |2Dq| replace:: :math:`N_q \times (3N)`

    *P* (*Q*) is a :class:`numpy.ndarray` of :math:`N_p` (:math:`N_q`) time
    steps, :math:`N` atoms, and :math:`3N` coordinates (e.g.,
    :meth:`MDAnalysis.core.AtomGroup.AtomGroup.coordinates`). *P* (*Q*) has
    either shape |3Dp| (|3Dq|), or |2Dp| (|2Dq|) in flattened form.

    :Arguments:
      *P*
         :class:`numpy.ndarray` representing a path
      *Q*
         :class:`numpy.ndarray` representing a path
      *N*
         int, number of atoms; if ``None`` then *P* and *Q* are each assumed to
         be a 3D :class:`numpy.ndarray`; else, they are 2D (flattened)
      *indices*
         generate indices of nearest neighbors rather than distances [``False``]

    :Returns:
      list of two pairs of numpy arrays, the first containing Hausdorff nearest
      neighbor indices for *P* and *Q*, the second containing Hausdorff
      nearest neighbor distances for *P* and *Q*

    Example::
     >>> u = Universe(PSF,DCD)
     >>> mid = len(u.trajectory)/2
     >>> ca = u.selectAtoms('name CA')
     >>> P = numpy.array([
     ...                ca.positions for _ in u.trajectory[:mid:]
     ...              ]) # first half of trajectory
     >>> Q = numpy.array([
     ...                ca.positions for _ in u.trajectory[mid::]
     ...              ]) # second half of trajectory
     >>> hausdorff(P,Q)
     4.7786639840135905
     >>> hausdorff(P,Q[::-1]) # hausdorff distance w/ reversed 2nd trajectory
     4.7786639840135905
    """
    if N is None:
        N = P.shape[1] # number of atoms from 2nd dim of P
        axis = (1,2)
    else:
        axis = 1
    d = numpy.array([sqnorm(p - Q, axis) for p in P])

    NN_indices_dists = []
    NN_indices_dists.append( numpy.argmin(d, axis=1), numpy.argmin(d, axis=0) )
    NN_indices_dists.append(
        numpy.amin(d, axis=1)/N)**0.5, (numpy.amin(d, axis=0)/N)**0.5 )
    return NN_indices_dists


class PSAPair(object):

    def __init__(self, idx, i, j, npaths):
        self.npaths = npaths
        self.matrixIdx = (i,j)
        self.pairIdx = idx #_dvec_idx(npaths,i,j)

        # Set by calling hausdorff_nn
        self.NNOfP = None# numpy array of length len(P) of NN (indices) in Q
        self.NNOfQ = None# numpy array of length len(Q) of NN (indices) in P
        self.NNDistancesOfP = None# numpy array of length len(P) of NN (distances) in P
        self.NNDistancesOfQ = None# numpy array of length len(P) of NN (distances) in Q

        # Set by self.getHausdorffPair
        self.hausdorffPair = (None, None)
        self.hausdorffDist = None


    def _dvec_idx(self, i, j):
        return self.npaths*i + j - (i+2)*(i+1)/2


    def runHausdorffNN(self, P, Q):
        r"""Convenience function to generate Hausdorff nearest neighbor lists
        of indices and distances for *this* pair of paths using `hausdorffNN`
        and setting instance variables of the PSAPair object.
        """
        NN_indices_dists = hausdorffNN(P, Q)
        self.NNOfP, self.NNOfQ = NN_indices_dists[0]
        self.NNDistancesOfP, self.NNDistancesOfQ = NN_indices_dists[1]
        return self


    def getHausdorffPair(self):
        maxNNDistP = max(self.NNDistancesOfP)
        maxNNDistQ = max(self.NNDistancesOfQ)
        if maxNNDistP > maxNNDistP:
            maxNNP = numpy.argmax(self.NNDistancesOfP)
            self.hausdorffPair = maxNNP, self.NNOfP[maxNNP]
            self.hausdorffDist = maxNNDistP
        else:
            maxNNQ = numpy.argmax(self.NNDistancesOfQ)
            self.hausdorffPair = self.NNOfQ[maxNNQ], maxNNQ
            self.hausdorffDist = maxNNDistQ

