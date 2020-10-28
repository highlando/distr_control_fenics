# this is rather optical checking
import dolfin
import scipy.sparse.linalg as spsla
import numpy as np
import matplotlib.pyplot as plt

import dolfin_navier_scipy.dolfin_to_sparrays as dts
import sadptprj_riclyap_adi.lin_alg_utils as lau
import distributed_control_fenics.cont_obs_utils as cou

N = 40


def check_output_opa(NYx=1, NYy=1, femp=None):

    NY = NYx*NYy
    if femp is None:
        # from dolfin_navier_scipy.problem_setups import cyl_fems
        # femp = cyl_fems(2)
        from dolfin_navier_scipy.problem_setups import drivcav_fems
        femp = drivcav_fems(N)

    V = femp['V']
    Q = femp['Q']

    odcoo = femp['odcoo']

    testcase = 1  # 1,2
    testcase = 2  # 1,2
    # testvelocities
    if testcase == 1:
        """case 1 -- not div free"""
        exv = dolfin.Expression(('x[1]', 'x[1]'), element=V.ufl_element())
        exv = dolfin.Expression(('x[0]*x[0]', 'x[0]*x[1]'),
                                element=V.ufl_element())
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
    MyC, My = cou.get_mout_opa(odcoo=odcoo, V=V, mfgrid=(NYy, NYx))

    ptmct = lau.app_prj_via_sadpnt(amat=stokesmats['M'],
                                   jmat=stokesmats['J'],
                                   rhsv=MyC.T,
                                   transposedprj=True)

    NV = V.dim()
    testvi = testv.vector().get_local().reshape((NV, 1))
    testvi0 = np.copy(testvi)
    testvi0 = lau.app_prj_via_sadpnt(amat=stokesmats['M'],
                                     jmat=stokesmats['J'],
                                     rhsv=testvi0)

    print("||J*v|| = {0}".format(np.linalg.norm(stokesmats['J']*testvi)))
    print("||J* v_df|| = {0}".format(np.linalg.norm(stokesmats['J']*testvi0)))

    # # testsignals from the test velocities
    testy = spsla.spsolve(My, MyC * testvi)
    testyv0 = spsla.spsolve(My, MyC * testvi0)
    testry = spsla.spsolve(My, np.dot(ptmct.T, testvi))

    print("||C v_df - C_df v|| = {0}".format(np.linalg.norm(testyv0 - testry)))

    yx = testy[:NY]
    yy = testy[NY:]

    plt.figure(1)
    plt.plot(yx)

    plt.figure(2)
    plt.plot(yy)

    plt.figure(3)
    dolfin.plot(testv, title='test velocity field')

    print('yx={0}'.format(yx))
    print('yy={0}'.format(yy))

    print('obsdom: [{0},{1}]x[{2},{3}]'.format(odcoo['xmin'], odcoo['xmax'],
                                               odcoo['ymin'], odcoo['ymax']))

    plt.show()


if __name__ == '__main__':
    check_output_opa(NYx=2, NYy=4)
