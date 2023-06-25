from ngsolve import *
import netgen.geom2d as g2d
from netgen.meshing import Pnt

geo = g2d.CSG2d()
geo.Add(g2d.Circle(Pnt((0,0)), 1.0))
ngmesh = geo.GenerateMesh(maxh=1.0)
mesh = Mesh(ngmesh)
mesh.ngmesh.SecondOrder()
#mesh.Curve(2)
Draw(mesh)

V = VectorH1(mesh, order=2)
u,v = V.TnT()

ndscal = V.ndof//mesh.dim
nv = mesh.nv

eps = lambda U : Grad(U)+Grad(U).trans

V.SetCouplingType(IntRange(0, V.ndof), COUPLING_TYPE.INTERFACE_DOF)
# for k,v in enumerate(mesh.vertices):
#     V.SetCouplingType(k, WIREBASKET_DOF)
#     V.SetCouplingType(ndscal+k, WIREBASKET_DOF)
V.SetCouplingType(IntRange(0, nv), COUPLING_TYPE.WIREBASKET_DOF)
V.SetCouplingType(IntRange(ndscal, ndscal+nv), COUPLING_TYPE.WIREBASKET_DOF)

a = BilinearForm(V, eliminate_internal=True)
a += InnerProduct(eps(u), eps(v))*dx
a += u*v *dx
c1 = Preconditioner(a, "local")
c2 = Preconditioner(a, "h1amg")
c3 = Preconditioner(a, "bddc")
a.Assemble()

gfu = GridFunction(V)

# gfu.vec[:] = 1.0

ndscal = V.ndof//mesh.dim
nv = mesh.nv
for k,v in enumerate(mesh.vertices):
    x,y = v.point
    gfu.vec[k] = -y
    gfu.vec[ndscal+k] = x

from ngsolve.la import EigenValues_Preconditioner
lams = EigenValues_Preconditioner(mat=a.mat, pre=c1)
print("jacobi")
print(lams[0:3], "...\n", lams[-3:])
print("precond1", max(lams) / min(lams))


lams = EigenValues_Preconditioner(mat=a.mat, pre=c2)
print("multigrid")
print(lams[0:3], "...\n", lams[-3:])
print("precond2", max(lams) / min(lams))

lams = EigenValues_Preconditioner(mat=a.mat, pre=c3)
print("bddc")
print(lams[0:3], "...\n", lams[-3:])
print("precond3", max(lams) / min(lams))


print("\nA-norm**2 = {}\n".format(InnerProduct(a.mat*gfu.vec,gfu.vec)))
print("\nc1-norm**2 = {}\n".format(InnerProduct(c1.mat*gfu.vec,gfu.vec)))
#print("\nc2-norm**2 = {}\n".format(InnerProduct(c2.mat*gfu.vec,gfu.vec)))
print("\nc3-norm**2 = {}\n".format(InnerProduct(c3.mat*gfu.vec,gfu.vec)))
Draw(gfu)
