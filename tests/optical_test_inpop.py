# this is rather optical checking
import dolfin
import numpy as np
import scipy.sparse.linalg as spsla
import matplotlib.pyplot as plt

import dolfin_navier_scipy.dolfin_to_sparrays as dts
import distr_control_fenics.cont_obs_utils as cou

dolfin.parameters.linear_algebra_backend = "Eigen"


def check_input_opa(NU, femp=None):

    if femp is None:
        # from dolfin_navier_scipy.problem_setups import cyl_fems
        # femp = cyl_fems(2)
        from dolfin_navier_scipy.problem_setups import drivcav_fems
        femp = drivcav_fems(20)

    V = femp['V']
    Q = femp['Q']

    cdcoo = femp['cdcoo']

    # get the system matrices
    stokesmats = dts.get_stokessysmats(V, Q)

    # check the B
    B1, Mu = cou.get_inp_opa(cdcoo=cdcoo, V=V, NU=NU, xcomp=0)
    B2, Mu = cou.get_inp_opa(cdcoo=cdcoo, V=V, NU=NU, xcomp=1)

    # get the rhs expression of Bu
    Bu1 = spsla.spsolve(stokesmats['M'],
                        B1*np.vstack([np.linspace(0, 1, NU).reshape((NU, 1)),
                                      np.linspace(0, 1, NU).reshape((NU, 1))]))

    Bu2 = spsla.spsolve(stokesmats['M'],
                        B2*np.vstack([np.linspace(0, 1, NU).reshape((NU, 1)),
                                      np.linspace(0, 1, NU).reshape((NU, 1))]))
    Bu3 = spsla.spsolve(stokesmats['M'], B1*np.vstack([1*np.ones((NU, 1)),
                                                       0.2*np.ones((NU, 1))]))

    bu1 = dolfin.Function(V)
    bu1.vector().set_local(Bu1)
    bu1.vector()[2] = 1  # for scaling and illustration purposes

    bu2 = dolfin.Function(V)
    bu2.vector().set_local(Bu2)
    bu2.vector()[2] = 1  # for scaling and illustration purposes

    bu3 = dolfin.Function(V)
    bu3.vector().set_local(Bu3)
    bu3.vector()[2] = 1  # for scaling and illustration purposes

    plt.figure(1)
    dolfin.plot(bu1, title='plot of Bu - extending in x')
    plt.figure(2)
    dolfin.plot(bu2, title='plot of Bu - extending in y')
    plt.figure(3)
    dolfin.plot(bu3, title='plot of Bu - extending in y')
    # dolfin.plot(V.mesh())

    dolfin.interactive()
    plt.show(block=False)

if __name__ == '__main__':
    check_input_opa(NU=4)
