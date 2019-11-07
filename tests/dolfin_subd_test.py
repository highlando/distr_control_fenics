import dolfin

N = 10
mesh = dolfin.UnitSquareMesh(N, N)
V = dolfin.FunctionSpace(mesh, 'CG', 1)
u = dolfin.TrialFunction(V)
cellmarkers = dolfin.MeshFunction('size_t', mesh, mesh.topology().dim())


class Omega0(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        return x[1] <= 0.5


class Omega1(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        return x[1] >= 0.5


lowerhalf = Omega0()
upperhalf = Omega1()

lowerhalf.mark(cellmarkers, 22)  # just 0 didnt work
dx = dolfin.Measure('dx', subdomain_data=cellmarkers)
alh = dolfin.assemble(u*dx(22)).get_local().sum()
print('area lower half: {0}'.format(alh))

upperhalf.mark(cellmarkers, 22)
dx = dolfin.Measure('dx', subdomain_data=cellmarkers)
auh = dolfin.assemble(u*dx(22)).get_local().sum()
print('area upperhalf: {0}'.format(auh))
