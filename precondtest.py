from ngsolve import *
from netgen.csg import *

cube = OrthoBrick( Pnt(0,0,0), Pnt(1,1,1) )
geo = CSGeometry()
geo.Add (cube)
ngmesh = geo.GenerateMesh(maxh=0.1)
mesh = Mesh(ngmesh)
meshOrder = 2
if meshOrder == 2:
    mesh.ngmesh.SecondOrder()
Draw(mesh)

if meshOrder == 1:
    fes = VectorH1(mesh, order=meshOrder) #add interleaved = True to get xyzxyz sorting
else:
    fes = VectorH1(mesh, order=meshOrder)

u = fes.TrialFunction()
v = fes.TestFunction()
bfK = BilinearForm(fes)
bfM = BilinearForm(fes)

def sigma(eps, mu, lam):
    return 2*mu*eps + lam*Trace(eps) * Id(eps.dims[0])

#für crankshaft
density = 7850
youngsModulus = 2.1e11 *1e-1
poissonsRatio = 0.3

E = youngsModulus
nu = poissonsRatio
rho = density

mu  = E / 2 / (1+nu) #Lame parameters
lam = E * nu / ((1+nu)*(1-2*nu))

bfK += InnerProduct(sigma(Sym(Grad(u)),mu,lam), Sym(Grad(v)))*dx
bfM += rho*u*v * dx

bfK.Assemble()
bfM.Assemble()

ndscal = fes.ndof // mesh.dim
nv = mesh.nv

fes.SetCouplingType(IntRange(0, fes.ndof), COUPLING_TYPE.INTERFACE_DOF)
fes.SetCouplingType(IntRange(0, nv), COUPLING_TYPE.WIREBASKET_DOF)
fes.SetCouplingType(IntRange(ndscal, ndscal + nv), COUPLING_TYPE.WIREBASKET_DOF)
fes.SetCouplingType(IntRange(2*ndscal, 2*ndscal + nv), COUPLING_TYPE.WIREBASKET_DOF)

with TaskManager():
    bfKbddc = BilinearForm(fes, eliminate_internal=True)
    bfKbddc += InnerProduct(sigma(Sym(Grad(u)), mu, lam), Sym(Grad(v))) *dx
    bfKbddc += 1e6*rho*u*v *dx
    prebddc = Preconditioner(bfKbddc, "bddc")
    bfKbddc.Assemble()

from ngsolve.la import EigenValues_Preconditioner
lams = EigenValues_Preconditioner(bfKbddc.mat, prebddc.mat)
print("bddc")
print("precond", max(lams) / min(lams))
print(lams[0:3], "...\n", lams[-3:])

if meshOrder > 1:
    fesLo = VectorH1(mesh, order=1)

    # create finite element spaces for mass matrix and stiffness matrix for lowest order
    uLo = fesLo.TrialFunction()
    vLo = fesLo.TestFunction()
    bfKLo = BilinearForm(fesLo)
    bfMLo = BilinearForm(fesLo)

    # setup (linear) mechanical FE-space for lowest order
    bfKLo += InnerProduct(sigma(Sym(Grad(uLo)), mu, lam), Sym(Grad(vLo))) *dx
    bfMLo += rho * uLo * vLo *dx

    bfKLo.Assemble()
    bfMLo.Assemble()

    # setup blocks for blockjacobi
def VertexPatchBlocks(fes):
    mesh = fes.mesh
    blocks = []
    freedofs = fes.FreeDofs()
    for v in mesh.vertices:
        vdofs = set()
        for el in mesh[v].elements:
            vdofs |= set(d for d in fes.GetDofNrs(el)
                            if freedofs[d])
        blocks.append(vdofs)
    return blocks

if meshOrder == 1:
    blocks = VertexPatchBlocks(fes)
    pre = bfK.mat.CreateBlockSmoother(blocks)
else:
    KM = bfKLo.mat.CreateMatrix()
    KM.AsVector().data = bfKLo.mat.AsVector() + 1e6* bfMLo.mat.AsVector()
    preLo = KM.Inverse(inverse='sparsecholesky')

    blocks = VertexPatchBlocks(fes)
    preHo = bfKbddc.mat.CreateBlockSmoother(blocks)

    emblo = Embedding(fes.ndof, IntRange(0, fesLo.ndof))
    pre = emblo @ preLo @ emblo.T
    pre = pre + preHo
    #pre = preHo

lams = EigenValues_Preconditioner(bfKbddc.mat, pre)
print("inv+hojac")
print("precond", max(lams) / min(lams))
print(lams[0:3], "...\n", lams[-3:])

print("jacobi")
prejac = bfKbddc.mat.CreateSmoother()
lams = EigenValues_Preconditioner(bfKbddc.mat, prejac)
print("precond", max(lams) / min(lams))
print(lams[0:3], "...\n", lams[-3:])