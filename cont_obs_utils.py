import dolfin
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spsla

import sadptprj_riclyap_adi.lin_alg_utils as lau

from dolfin import dx, inner

__all__ = ['get_inp_opa',
           'get_mout_opa',
           'app_difffreeproj',
           'get_regularized_c',
           'Cast1Dto2D',
           'get_rightinv',
           'extract_output',
           'CharactFun',
           'get_pavrg_onsubd']


def get_inp_opa(cdcoo=None, NU=8, V=None, xcomp=0):
    """dolfin.assemble the 'B' matrix

    the findim array representation
    of the input operator

    Parameters
    ----------
    cdcoo : dictionary with xmin, xmax, ymin, ymax of the control domain
    NU : dimension of the input space
    V : FEM space of the velocity
    xcomp : spatial component that is preserved, the other is averaged out

    """

    cdom = ContDomain(cdcoo)

    v = dolfin.TestFunction(V)
    v_one = dolfin.Expression(('1', '1'), element=V.ufl_element())
    v_one = dolfin.interpolate(v_one, V)

    BX, BY = [], []

    for nbf in range(NU):
        ubf = L2abLinBas(nbf, NU)
        bux = Cast1Dto2D(ubf, cdom, vcomp=0, xcomp=xcomp)
        buy = Cast1Dto2D(ubf, cdom, vcomp=1, xcomp=xcomp)
        bx = inner(v, bux) * dx
        by = inner(v, buy) * dx
        Bx = dolfin.assemble(bx)
        By = dolfin.assemble(by)
        Bx = Bx.get_local()
        By = By.get_local()
        Bx = Bx.reshape(len(Bx), 1)
        By = By.reshape(len(By), 1)
        BX.append(sps.csc_matrix(By))
        BY.append(sps.csc_matrix(Bx))

    Mu = ubf.massmat()

    try:
        return (
            sps.hstack([sps.hstack(BX), sps.hstack(BY)], format='csc'),
            sps.block_diag([Mu, Mu]))
    except AttributeError:  # e.g. in scipy <= 0.9
        return (
            sps.hstack([sps.hstack(BX), sps.hstack(BY)], format='csc'),
            sps.hstack([sps.vstack([Mu, sps.csc_matrix((NU, NU))]),
                        sps.vstack([sps.csc_matrix((NU, NU)), Mu])]))


def get_mout_opa(odcoo=None, NY=None, V=None,
                 msrfncs='pc', mfgrid=(1, 1)):
    """dolfin.assemble the 'MyC' matrix

    the find an array representation
    of the output operator

    the considered output is

    .. math::
        y(t) = M_y^{-1} \\int_{\\Omega_o} \\Psi v(t) dx

    where :math:`\\Psi` is the formal vector of measurement functions
    that have support on the domain of observation :math:`\\Omega_0`

    Arrangement of the `mfgrid(3,2)` vector array is matrix like

    s11 s21
    s12 s22
    s13 s23   in physical space

    and then output will be `[y(s11), y(s12), ..., y(s23)].T`


    Parameters
    ----------
    msrfncs : string, optional
        type of the measurement functions, defaults to `'pc'` -- piecewise
        constant (which measures averages over subdomains)
    mfgrid : tuple, optional
        grid for the definition of the `pc` measurments functions
    odcoo : dictonary
        with keys `'xmin'`, `'xmax'`, `'ymin'`, `'ymax'` -- the coordinates
        of the domain of observation
    """

    if NY is not None:
        raise NotImplementedError('by now only piecewise constant defined ' +
                                  ' via mfgrid')
    NV = V.dim()

    # define the subdomains == sensor locations
    xpartitions = mfgrid[1]
    xmin, xspan = odcoo['xmin'], odcoo['xmax']-odcoo['xmin']
    xspspan = xspan / xpartitions

    ypartitions = mfgrid[0]
    ymax, yspan = odcoo['ymax'], odcoo['ymax']-odcoo['ymin']
    yspspan = yspan / ypartitions

    spcoolist = []  # list of subdomain coordinates
    for nxp in range(xpartitions):
        spxmin = xmin + nxp*xspspan
        spxmax = xmin + (nxp+1)*xspspan
        for nyp in range(ypartitions):
            spymax = ymax - nyp*yspspan
            spymin = ymax - (nyp+1)*yspspan
            spcoo = dict(xmin=spxmin, xmax=spxmax,
                         ymin=spymin, ymax=spymax)
            spcoolist.append(spcoo)

    # filter functions that filter x, y components
    voney = dolfin.Expression(('0', '1'), element=V.ufl_element())
    vonex = dolfin.Expression(('1', '0'), element=V.ufl_element())
    voney = dolfin.interpolate(voney, V)
    vonex = dolfin.interpolate(vonex, V)

    u = dolfin.TrialFunction(V)
    cellmarkers = dolfin.MeshFunction('size_t', V.mesh(),
                                      V.mesh().topology().dim())

    cxll, cyll = [], []
    strtdx = 101
    for spcoo in spcoolist:
        strtdx += 1
        spodom = ContDomain(spcoo)
        spodom.mark(cellmarkers, strtdx)  # just 0 didnt work
        dx = dolfin.Measure('dx', subdomain_data=cellmarkers)
        cxl = dolfin.assemble(inner(vonex, u) * dx(strtdx)).get_local()
        cxll.append(sps.csr_matrix(cxl.reshape((1, NV))))
        cyl = dolfin.assemble(inner(voney, u) * dx(strtdx)).get_local()
        cyll.append(sps.csr_matrix(cyl.reshape((1, NV))))

    mycx = sps.vstack(cxll)
    mycy = sps.vstack(cyll)
    MyC = sps.vstack([mycx, mycy])

    My = (xspspan*yspspan)*sps.eye(2*xpartitions*ypartitions)  # the volumes

    return MyC, My


def app_difffreeproj(v=None, J=None, M=None):
    """apply the regularization (projection to divfree vels)

    i.e. compute v = [I-M^-1*J.T*S^-1*J]v
    """

    vg = lau.app_schurc_inv(M, J, np.atleast_2d(J * v).T)

    vg = spsla.spsolve(M, J.T * vg)

    return v - vg


def get_regularized_c(Ct=None, J=None, Mt=None):
    """apply the regularization (projection to divfree vels)

    i.e. compute rC = C*[I-M^-1*J.T*S^-1*J] as
    rCT = [I - J.T*S.-T*J*M.-T]*C.T
    """

    raise UserWarning('deprecated - use more explicit approach to proj via ' +
                      'sadpoints systems as implemented in linalg_utils')

    Nv, NY = Mt.shape[0], Ct.shape[1]
    try:
        rCt = np.load('data/regCNY{0}vdim{1}.npy'.format(NY, Nv))
    except IOError:
        print('no data/regCNY{0}vdim{1}.npy'.format(NY, Nv))
        MTlu = spsla.factorized(Mt)
        auCt = np.zeros(Ct.shape)
        # M.-T*C.T
        for ccol in range(NY):
            auCt[:, ccol] = MTlu(np.array(Ct[:, ccol].todense())[:, 0])
        # J*M.-T*C.T
        auCt = J * auCt
        # S.-T*J*M.-T*C.T
        auCt = lau.app_schurc_inv(MTlu, J, auCt)
        rCt = Ct - J.T * auCt
        np.save('data/regCNY{0}vdim{1}.npy'.format(NY, Nv), rCt)

    return np.array(rCt)


# Subdomains of Control and Observation
class ContDomain(dolfin.SubDomain):

    def __init__(self, ddict):
        dolfin.SubDomain.__init__(self)
        self.minxy = [ddict['xmin'], ddict['ymin']]
        self.maxxy = [ddict['xmax'], ddict['ymax']]

    def inside(self, x, on_boundary):
        insi = (dolfin.between(x[0], (self.minxy[0], self.maxxy[0]))
                and
                dolfin.between(x[1], (self.minxy[1], self.maxxy[1])))
        # print(x, insi)
        return insi


class L2abLinBas():
    """ return the hat function related to the num-th vertex

    from the interval [a=0, b=1] with an equispaced grid
    of N vertices

    """

    def __init__(self, num, N, a=0.0, b=1.0):
        self.dist = (b - a) / (N - 1)
        self.vertex = a + num * self.dist
        self.num, self.N = num, N
        self.a, self.b = a, b

    def evaluate(self, s):
        # print s
        if max(self.a, self.vertex - self.dist) <= s <= self.vertex:
            sval = 1.0 - 1.0 / self.dist * (self.vertex - s)
        elif self.vertex <= s <= min(self.b, self.vertex + self.dist):
            sval = 1.0 - 1.0 / self.dist * (s - self.vertex)
        else:
            sval = 0
        return sval

    def massmat(self):
        """ return the mass matrix
        """
        mesh = dolfin.IntervalMesh(self.N - 1, self.a, self.b)
        Y = dolfin.FunctionSpace(mesh, 'CG', 1)
        yv = dolfin.TestFunction(Y)
        yu = dolfin.TrialFunction(Y)
        my = yv * yu * dx
        my = dolfin.assemble(my)
        return dolfin.as_backend_type(my).sparray()


def Cast1Dto2D(u, cdom, vcomp=None, xcomp=0, degree=2):
    """ casts a function u defined on [u.a, u.b]

    into the f[comp] of an expression
    defined on a 2D domain cdom by
    by scaling to fit the xcomp extension
    and simply extruding into the other direction
    """

    # control 1D basis function
    u = u
    # domain of control
    cdom = cdom
    # component of the value to be set as u(s)
    vcomp = vcomp
    # component of x to be considered as s coordinate
    xcomp = xcomp
    # transformation of the intervals [cd.xmin, cd.xmax] -> [a, b]
    # via s = m*x + d
    m = (u.b - u.a) / \
        (cdom.maxxy[xcomp] - cdom.minxy[xcomp])
    d = u.b - m * cdom.maxxy[xcomp]

    class IDtoIIDExpr(dolfin.UserExpression):

        def eval(self, value, x):
            if cdom.inside(x, False):
                if xcomp is None:
                    value[:] = u.evaluate(m * x[xcomp] + d)
                else:
                    value[:] = 0
                    value[vcomp] = u.evaluate(
                        m * x[xcomp] + d)
            else:
                value[:] = 0

        def value_shape(self):
            return (2,)

    return IDtoIIDExpr(degree=degree)


def get_rightinv(C):
    """compute the rightinverse bmo SVD

    """
    # use numpy routine for dense matrices
    try:
        u, s, vt = np.linalg.svd(np.array(C.todense()), full_matrices=0)
    except AttributeError:
        u, s, vt = np.linalg.svd(C, full_matrices=0)

    return np.dot(vt.T, np.dot(np.diag(1.0 / s), u.T))


def get_vstar(C, ystar, odcoo, NY):

    ystarvec = get_ystarvec(ystar, odcoo, NY)
    Cgeninv = get_rightinv(C)

    return np.dot(Cgeninv, ystarvec)


def get_ystarvec(ystar, odcoo, NY):
    """get the vector of the current target signal

    """
    ymesh = dolfin.IntervalMesh(NY - 1, odcoo['ymin'], odcoo['ymax'])
    Y = dolfin.FunctionSpace(ymesh, 'CG', 1)

    ystarvec = np.zeros((NY * len(ystar), 1))
    for k, ysc in enumerate(ystar):
        cyv = dolfin.interpolate(ysc, Y)
        ystarvec[k * NY:(k + 1) * NY, 0] = cyv.vector().get_local()

    return ystarvec


def extract_output(strdict=None, tmesh=None, c_mat=None,
                   ystarvec=None, load_data=None):
    """extract the output `y` by applying `C` to the data

    returns lists of lists to be plotted pickled in json files

    Parameters
    ----------
    strdict : dictionary
        with time as the keys and path to data as values
    tmesh : iterable list or ndarray
        values of the time instances
    c_mat : (K,N) sparse matrix
        output operator
    ystarvec : callable f(t), optional
        returns a (K,) or (K,1) vector at time `t`
    load_data : callable f(string)
        returns the data as ndarray

    Returns
    -------
    yscomplist : list
        of the outputs at `tmesh`, where `yscomplist[0] = y(trange[0])`
    ystarlist : list, optional
        of ystarvec values as yscomplist

    """

    cur_v = load_data(strdict[tmesh[0]])
    yn = c_mat.dot(cur_v)
    yscomplist = [yn.flatten().tolist()]

    for t in tmesh[1:]:
        cur_v = load_data(strdict[t])
        yn = c_mat.dot(cur_v)
        yscomplist.append(yn.flatten().tolist())

    if ystarvec is None:
        return yscomplist

    else:
        ystarlist = [ystarvec(0).flatten().tolist()]
        for t in tmesh[1:]:
            ystarlist.append(ystarvec(t).flatten().tolist())

    return yscomplist, ystarlist


def CharactFun(subdom, degree=2):
    """ characteristic function of subdomain """
    class xifunexp(dolfin.UserExpression):

        def eval(self, value, x):
            if subdom.inside(x, False):
                value[:] = 1
            else:
                value[:] = 0

        def value_shape(self):
            return (1,)

    return xifunexp(degree=degree)


def get_pavrg_onsubd(odcoo=None, Q=None, ppin=None):
    """assemble matrix that returns the pressure average over a subdomain

    TODO: deprecate this -- use get_mout_opa instead
    """

    prsodom = ContDomain(odcoo)
    q = dolfin.TrialFunction(Q)
    # p = dolfin.TestFunction(Q)

    # factor to compute the average via \bar u = 1/h \int_0^h u(x) dx
    Ci = 1.0 / ((odcoo['xmax'] - odcoo['xmin']) *
                (odcoo['ymax'] - odcoo['ymin']))

    cellmarkers = dolfin.MeshFunction('size_t', Q.mesh(),
                                      Q.mesh().topology().dim())

    prsdx = 303
    prsodom.mark(cellmarkers, prsdx)  # just 0 didnt work
    dx = dolfin.Measure('dx', subdomain_data=cellmarkers)
    cp = dolfin.assemble(Ci*q*dx(prsdx)).get_local().reshape((1, Q.dim()))
    ccp = sps.csc_matrix(cp)

    if ppin is None:
        return ccp  # np.atleast_2d(cp.array())
    else:
        raise UserWarning('Need to implement/specify the pinned pressure')
