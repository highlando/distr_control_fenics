import unittest
import dolfin
import scipy.sparse.linalg as spsla
import numpy as np

import distributed_control_fenics.cont_obs_utils as cou


import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=UserWarning)


# unittests for the input and output operators


class TestInoutOpas(unittest.TestCase):

    # @unittest.skip("for now")
    def test_outopa_workingconfig(self):
        """ The innerproducts that assemble the output operator

        are accurately sampled for this parameter set NV=40 and NYx=NYy=2
        """

        NV = 40
        # NYx, NYy = 2, 2
        NYx, NYy = 2, 2
        NY = NYx*NYy

        mesh = dolfin.UnitSquareMesh(NV, NV)
        V = dolfin.VectorFunctionSpace(mesh, "CG", 2)

        xmin, xmax, ymin, ymax = .45, .55, .6, .8
        odcoo = dict(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
        MyC, My = cou.get_mout_opa(odcoo=odcoo, V=V, mfgrid=(NYy, NYx))

        exv = dolfin.Expression(('x[0]', 'x[0]*x[1]'),
                                element=V.ufl_element())
        testv = dolfin.interpolate(exv, V)

        testvi = testv.vector().get_local()
        testy = spsla.spsolve(My, MyC * testvi)

        # for linear funcs the average is the value at the midpoints
        xspan = xmax-xmin
        yspan = ymax-ymin
        mpxone, mpyone = xmin+.25*xspan, ymax-.25*yspan
        mpxfour, mpyfour = xmin+.75*xspan, ymax-.75*yspan
        vmpxone, vmpyone = exv((mpxone, mpyone))
        vmpxfour, vmpyfour = exv((mpxfour, mpyfour))

        self.assertTrue(np.allclose(vmpxone, testy[0]) and
                        np.allclose(vmpyone, testy[NY]) and
                        np.allclose(vmpxfour, testy[NY-1]) and
                        np.allclose(vmpyfour, testy[-1]))

suite = unittest.TestLoader().loadTestsFromTestCase(TestInoutOpas)
unittest.TextTestRunner(verbosity=2).run(suite)
