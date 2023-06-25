from ngsolve import *
from netgen.csg import *
import numpy as np

#parameters

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

def RotationMatrixZ(angleRad):
    return np.array([ [np.cos(angleRad),-np.sin(angleRad), 0],
                      [np.sin(angleRad), np.cos(angleRad), 0],
                      [0,	    0,        1] ]);

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

def GenerateConrod(zOff):
    ey0 = [0,1,0] #top/bottom face vector of conrod
    ey1 = [0,-1,0]

    ex0 = [1,0,0] #top/bottom face vector of conrod
    ex1 = [1,0,0]
    
    ey0 = RotationMatrixZ(-angC)@ey0
    ey1 = RotationMatrixZ(angC)@ey1
    ex0 = RotationMatrixZ(-angC)@ex0
    ex1 = RotationMatrixZ(angC)@ex1


    pl1 = Plane(Pnt(0, 0.5*bc,0),Vec(ey0[0],ey0[1],ey0[2]))
    pl2 = Plane(Pnt(0,-0.5*bc,0),Vec(ey1[0],ey1[1],ey1[2]))

    pl3 = Plane(Pnt(-0.5*lc,0,0),Vec(-1,0,0))
    pl4 = Plane(Pnt( 0.5*lc,0,0),Vec( 1,0,0))

    pl5 = Plane(Pnt( 0,0,-0.5*dc+zOff),Vec( 0,0,-1))
    pl6 = Plane(Pnt( 0,0, 0.5*dc+zOff),Vec( 0,0, 1))

    
    cylC1 = Cylinder(Pnt(-0.5*lc,0,-1), Pnt(-0.5*lc,0,1), r1)
    #cylC1o = Cylinder(Pnt(-0.5*lc,0,-1), Pnt(-0.5*lc,0,1), r1o)
    cylC1o = Sphere(Pnt(-0.5*lc,0,zOff), r1o) #in fact is a sphere

    cylC2 = Cylinder(Pnt( 0.5*lc,0,-1), Pnt( 0.5*lc,0,1), r2)
    #cylC2o = Cylinder(Pnt(0.5*lc,0,-1), Pnt( 0.5*lc,0,1), r2o)
    cylC2o = Sphere(Pnt(0.5*lc,0,zOff), r2o) #in fact is a sphere

    cylSideA = (Cylinder(Pnt(-0.5*lc+r1o,0,cylOffZ+zOff), Pnt(0.5*lc-r2o,0,cylOffZ+zOff), cylR)*
                Plane(Pnt(-0.5*lc+r1o-0.002,0,0),Vec(-1,0,0))*
                Plane(Pnt( 0.5*lc-r2o+0.002,0,0),Vec( 1,0,0)))

    cylSideB = (Cylinder(Pnt(-0.5*lc+r1o,0,-cylOffZ+zOff), Pnt(0.5*lc-r2o,0,-cylOffZ+zOff), cylR)*
                Plane(Pnt(-0.5*lc+r1o-0.002,0,0),Vec(-1,0,0))*
                Plane(Pnt( 0.5*lc-r2o+0.002,0,0),Vec( 1,0,0)))


    return ((pl1*pl2*pl3*pl4+cylC1o+cylC2o)-cylC1-cylC2)*pl5*pl6-cylSideA-cylSideB
    #return pl1*pl2*pl3*pl4*pl5*pl6

geoConrod = CSGeometry()
conrod = GenerateConrod(0)#db+dk+db+0.5*b1
geoConrod.Add(conrod)

meshSize = 0.0075
mesh = Mesh( geoConrod.GenerateMesh(maxh=meshSize+0.001*0))
Draw(mesh)