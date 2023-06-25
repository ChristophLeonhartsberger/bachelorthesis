from netgen.geom2d import unit_square
from ngsolve import *
import netgen.geom2d as g2d
from netgen.meshing import Pnt

geo = g2d.CSG2d()
geo.Add(g2d.Circle(Pnt((0,0)), 1.0))
ngmesh = geo.GenerateMesh(maxh=0.3)
mesh = Mesh(ngmesh)
fes_ho = Discontinuous(H1(mesh, order=10))
fes_lo = H1(mesh, order=1, dirichlet=".*")
fes_lam = Discontinuous(H1(mesh, order=1))
fes = fes_ho*fes_lo*fes_lam
uho, ulo, lam = fes.TrialFunction()

a = BilinearForm(fes)
a += Variation(0.5 * grad(uho)*grad(uho)*dx
               - 1*uho*dx
               + (uho-ulo)*lam*dx(element_vb=BBND))
gfu = GridFunction(fes)
solvers.Newton(a=a, u=gfu)
Draw(gfu.components[0],deformation=True)