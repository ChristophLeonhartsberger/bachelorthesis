from ngsolve import *
from netgen.csg import *
import numpy as np

#parameters
meshSize = 0.0075
#crank:
b1 = 0.012 #width of journal bearing
r1 = 0.012 #radius of journal bearing
dk = 0.015 #crank arm width (z)
bk = 0.032 #crank arm size (y)

l3 = 0.030
l4 = 0.040
#l4x= 0.005 #offset of counterweight
lk = 0.030 #l4*0.5+l3 #crank arm length (x)
bm = 0.065
dBevel = dk*0.5
#shaft:
r0 = 0.012 #0.012
d0 = 0.020 #shaft length at left/right support
d1 = 0.012 #shaft length at intermediate support

#distance rings:
db = 0.002          #width of distance ring
rdb0 = r0+db        #total radius of distance ring, shaft
rdb1 = r1+db        #total radius of distance ring, crank

#conrod:
bc = 0.024      #height of conrod
dc = 0.012      #width of conrod
lc = 0.080      #length of conrod (axis-axis)
r1o= r1+0.006   #outer radius of conrod at crank joint
r2 = 0.008      #radius of piston journal bearing
r2o= r2+0.006   #outer radius of conrod at piston joint

cylOffZ=0.010  #z-offset of cylinder cut out of conrod
cylR = 0.008    #radius of cylinder cut out of conrod

angC = 4*np.pi/180

#piston:
dpb = r2o-0.000   #axis inside piston
r2p = r2o+0.004   #0.018
lp = 0.034
bp = 0.050
lpAxis = dc+2*db
lOffCut = 0.011 #offset for cutout of big cylinder

#total length of one segment:
lTotal = db+dk+db+b1+db+dk+db+d1

#eps
eps = 5e-4 #added to faces, to avoid CSG-problems
def NormL2(vector):
    value = 0
    for x in vector:
        value += x**2
    return value**0.5

def Normalize(v):
    v2=[0]*len(v)

    fact = NormL2(v)
    fact = 1./fact
    for i in range(len(v2)): 
        v2[i]=fact*v[i]
    return v2

def VAdd(v0, v1):
    if len(v0) != len(v1): print("ERROR in VAdd: incompatible vectors!")
    n = len(v0)
    v = [0]*n
    for i in range(n):
        v[i] = v0[i]+v1[i]
    return v

def VSub(v0, v1):
    if len(v0) != len(v1): print("ERROR in VSub: incompatible vectors!")
    n = len(v0)
    v = [0]*n
    for i in range(n):
        v[i] = v0[i]-v1[i]
    return v

def CSGcube(pCenter,size):
    s2 = [0.5*size[0],0.5*size[1],0.5*size[2]]
    p0 = VSub(pCenter,s2)
    p1 = VAdd(pCenter,s2)
    brick = OrthoBrick(Pnt(p0[0],p0[1],p0[2]),Pnt(p1[0],p1[1],p1[2]))
    return brick

def CSGcylinder(p0,p1,r):
    v = VSub(p1,p0)
    v = Normalize(v)
    cyl = Cylinder(Pnt(p0[0],p0[1],p0[2]), Pnt(p1[0],p1[1],p1[2]), 
                   r) * Plane(Pnt(p0[0],p0[1],p0[2]), Vec(-v[0],-v[1],-v[2])) * Plane(Pnt(p1[0],p1[1],p1[2]), Vec(v[0],v[1],v[2])) 
    return cyl

def GeneratePiston(zOff):
    p0 = [-dpb,0,zOff]
    p1 = [-dpb+lp,0,zOff]
    cylPo   = CSGcylinder(p0, p1, 0.5*bp) #piston outside
    cylPaxis= CSGcylinder([0,0,-0.5*lpAxis-eps+zOff],     [0,0, 0.5*lpAxis+eps+zOff], r2) #piston axis
    cylPaxis0= CSGcylinder([0,0,-0.5*lpAxis-eps+zOff],    [0,0,-0.5*lpAxis+db+zOff], r2+db) #piston axis
    cylPaxis1= CSGcylinder([0,0, 0.5*lpAxis-db+zOff], [0,0, 0.5*lpAxis+eps+zOff], r2+db) #piston axis
    cylPin  = CSGcylinder([0,0,-0.5*lpAxis+zOff], [0,0, 0.5*lpAxis+zOff], r2p) #piston inner cutout

    #box = CSGcube([0,0,zOff], [dpb+r2p,2*(r2p),lpAxis])
    box = CSGcube([-0.5*dpb,0,zOff], [dpb,2*(r2p)-0.002,lpAxis-0.000])

    cylCut  = CSGcylinder([-(l4+l3+lOffCut),0,-bp+zOff], [-(l4+l3+lOffCut),0, bp+zOff], l4+l3) #piston inner cutout

    return (cylPo-box-cylCut-cylPin)+cylPaxis+cylPaxis0+cylPaxis1

geoPiston = CSGeometry()
piston = GeneratePiston(0)#db+dk+db+0.5*b1
geoPiston.Add(piston)

meshPiston = Mesh( geoPiston.GenerateMesh(maxh=meshSize+0.001*0))
Draw(meshPiston)