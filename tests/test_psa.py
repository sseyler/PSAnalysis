from __future__ import print_function
import MDAnalysis.analysis.psa
from MDAnalysis import SelectionError, SelectionWarning, FinishTimeException

from numpy.testing import *
import numpy as np
import nose
from nose.plugins.attrib import attr

import os
import shutil
import errno
import tempfile
import itertools
import warnings

from MDAnalysisTests.datafiles import PSF, DCD, DCD2
from MDAnalysisTests import executable_not_found_runtime

import unittest

class TestPSAnalysis(TestCase):
    def setUp(self):
        self.iu1 = np.triu_indices(3, k=1)
        self.universe1 = MDAnalysis.Universe(PSF, DCD)
        self.universe2 = MDAnalysis.Universe(PSF, DCD2)
        self.universe_rev = MDAnalysis.Universe(PSF, DCD)
        self.universes = [self.universe1, self.universe2, self.universe_rev]
        self.psa = MDAnalysis.analysis.psa.PSAnalysis(self.universes,           \
                                               path_select='name CA')
        self.psa.generate_paths(align=True)
        self.psa.paths[-1] = self.psa.paths[-1][::-1,:,:] # reverse third path
        self._run()

    def _run(self):
        self.psa.run(metric='hausdorff')
        self.hausd_matrix = self.psa.get_pairwise_distances()
        self.psa.run(metric='discrete_frechet')
        self.frech_matrix = self.psa.get_pairwise_distances()
        self.hausd_dists = self.hausd_matrix[self.iu1]
        self.frech_dists = self.frech_matrix[self.iu1]

    def tearDown(self):
        del self.universe1
        del self.universe2
        del self.universe_rev
        del self.psa
        shutil.rmtree('psadata') # remove psa data directory

    def test_hausdorff_bound(self):
        err_msg = "Some Frechet distances are smaller than corresponding "      \
                + "Hausdorff distances"
        assert_array_less(self.hausd_dists, self.frech_dists, err_msg)

    def test_reversal_hausdorff(self):
        err_msg = "Hausdorff distances changed after path reversal"
        assert_array_almost_equal(self.hausd_matrix[1,2],                       \
                                  self.hausd_matrix[0,1],                       \
                                  decimal=3, err_msg=err_msg)

    def test_reversal_frechet(self):
        err_msg = "Frechet distances did not increase after path reversal"
        assert self.frech_matrix[1,2] >= self.frech_matrix[0,1], err_msg

    def test_avg_metric_bound(self):
        pass

    def test_alignment(self):
        pass

    def test_neighbor_sequence(self):
        pass

    def test_clustering(self):
        pass

if __name__ == '__main__':
    unittest.main()