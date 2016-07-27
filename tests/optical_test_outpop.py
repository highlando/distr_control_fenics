# this is rather optical checking
import dolfin
import scipy.sparse.linalg as spsla
import numpy as np
import matplotlib.pyplot as plt

import dolfin_navier_scipy.dolfin_to_sparrays as dts
import sadptprj_riclyap_adi.lin_alg_utils as lau
import distr_control_fenics.cont_obs_utils as cou

dolfin.parameters.linear_algebra_backend = "Eigen"


def check_output_opa(NY, femp=None):

    if femp is None:
        # from dolfin_navier_scipy.problem_setups import cyl_fems
        # femp = cyl_fems(2)
        from dolfin_navier_scipy.problem_setups import drivcav_fems
        femp = drivcav_fems(20)

    V = femp['V']
    Q = femp['Q']

    odcoo = femp['odcoo']

    testcase = 2  # 1,2
    # testvelocities
    if testcase == 1:
        """case 1 -- not div free"""
        exv = dolfin.Expression(('x[1]', 'x[1]'), element=V.ufl_element())
    if testcase == 2:
        """case 2 -- disc div free"""
        exv = dolfin.Expression(('1', '1'), element=V.ufl_element())

    testv = dolfin.interpolate(exv, V)

    odcoo = femp['odcoo']

    stokesmats = dts.get_stokessysmats(V, Q, nu=1)
    # remove the freedom in the pressure
    stokesmats['J'] = stokesmats['J'][:-1, :][:, :]
    stokesmats['JT'] = stokesmats['JT'][:, :-1][:, :]

    bc = dolfin.DirichletBC(V, exv, 'on_boundary')

    # reduce the matrices by resolving the BCs
    (stokesmatsc,
     rhsd_stbc,
     invinds,
     bcinds,
     bcvals) = dts.condense_sysmatsbybcs(stokesmats, [bc])

    # check the C
    MyC, My = cou.get_mout_opa(odcoo=odcoo, V=V, NY=NY)
    # MyC = MyC[:, invinds][:, :]

    # signal space
    ymesh = dolfin.IntervalMesh(NY - 1, odcoo['ymin'], odcoo['ymax'])
    Y = dolfin.FunctionSpace(ymesh, 'CG', 1)

    y1 = dolfin.Function(Y)
    y2 = dolfin.Function(Y)
    # y3 = dolfin.Function(Y)

    # dolfin.plot(V.mesh())

    ptmct = lau.app_prj_via_sadpnt(amat=stokesmats['M'],
                                   jmat=stokesmats['J'],
                                   rhsv=MyC.T,
                                   transposedprj=True)

    testvi = testv.vector().array()
    testvi0 = np.atleast_2d(testv.vector().array()).T
    testvi0 = lau.app_prj_via_sadpnt(amat=stokesmats['M'],
                                     jmat=stokesmats['J'],
                                     rhsv=testvi0)

    print "||J*v|| = {0}".format(np.linalg.norm(stokesmats['J'] * testvi))
    print "||J* v_df|| = {0}".format(np.linalg.norm(stokesmats['J'] * testvi0))

    # # testsignals from the test velocities
    testy = spsla.spsolve(My, MyC * testvi)
    testyv0 = spsla.spsolve(My, MyC * testvi0)
    # testyg = spsla.spsolve(My, MyC * (testvi.flatten() - testvi0))
    testry = spsla.spsolve(My, np.dot(ptmct.T, testvi))

    print "||C v_df - C_df v|| = {0}".format(np.linalg.norm(testyv0 - testry))

    plt.figure(1)
    y1.vector().set_local(testy[:NY])
    dolfin.plot(y1, title="x-comp of C*v")

    plt.figure(2)
    y2.vector().set_local(testy[NY:])
    dolfin.plot(y2, title="y-comp of C*v")

    # y2.vector().set_local(testyv0[:NY])
    # dolfin.plot(y2, title="x-comp of $C*(P_{df}v)$")

    # y3.vector().set_local(testyg[:NY])
    # dolfin.plot(y3, title="x-comp of $C*(v - P_{df}v)$")

    plt.show(block=False)

if __name__ == '__main__':
    check_output_opa(NY=4)
