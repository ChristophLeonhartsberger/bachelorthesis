from ngsolve import *
import netgen.geom2d as g2d
from netgen.meshing import Pnt

#from netgen.csg import *

#cube = OrthoBrick( Pnt(0,0,0), Pnt(1,1,1) )
#geo = CSGeometry()
#geo.Add (cube)
#ngmesh = geo.GenerateMesh(maxh=0.1)
#mesh = Mesh(ngmesh)

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
ne = mesh.ne
nedges = len(mesh.edges)
print("nv:", nv)
print("ne:", ne)
print("nedges:", nedges)
print("ndof:",V.ndof)
print("ndscal:",ndscal)

#for i in range(V.ndof):
#    print(V.CouplingType(i))
#
#print("+++++++++++++++++++++")

V.SetCouplingType(IntRange(0, V.ndof), COUPLING_TYPE.LOCAL_DOF)
# for k,v in enumerate(mesh.vertices):
#     V.SetCouplingType(k, WIREBASKET_DOF)
#     V.SetCouplingType(ndscal+k, WIREBASKET_DOF)
V.SetCouplingType(IntRange(0, nv), COUPLING_TYPE.WIREBASKET_DOF)
V.SetCouplingType(IntRange(ndscal, ndscal+nv), COUPLING_TYPE.WIREBASKET_DOF)

#for i in range(V.ndof):
#    print(V.CouplingType(i))

eps = lambda U : Grad(U)+Grad(U).trans
a = BilinearForm(V, eliminate_internal=True)
a += InnerProduct(eps(u), eps(v))*dx
a.Assemble()

gfu = GridFunction(V)


# gfu.vec[:] = 1.0

ndscal = V.ndof//mesh.dim
nv = mesh.nv
for k,v in enumerate(mesh.vertices):
    x,y = v.point
    gfu.vec[k] = -y
    gfu.vec[ndscal+k] = x


print("\nA-norm**2 = {}\n".format(InnerProduct(a.mat*gfu.vec,gfu.vec)))
# print("\nc-norm**2 = {}\n".format(InnerProduct(c.mat*gfu.vec,gfu.vec)))
Draw(gfu)
