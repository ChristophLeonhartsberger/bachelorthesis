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

#points
pLB = [0 ,0,-d0]
p0B = [0 ,0,0]
p1B = [0 ,0,db]
#p2B = [0, 0,db+dk]
p21B =[lk,0,db+dk]
p31B = [lk,0,db+dk+db]
p41B = [lk,0,db+dk+db+b1]
p51B =[lk,0,db+dk+db+b1+db]
p6B = [0 ,0,db+dk+db+b1+db+dk]
p7B = [0 ,0,db+dk+db+b1+db+dk+db]
p8B = [0 ,0,lTotal]

def CSGcylinder(p0,p1,r):
    v = VSub(p1,p0)
    v = Normalize(v)
    cyl = Cylinder(Pnt(p0[0],p0[1],p0[2]), Pnt(p1[0],p1[1],p1[2]), 
                   r) * Plane(Pnt(p0[0],p0[1],p0[2]), Vec(-v[0],-v[1],-v[2])) * Plane(Pnt(p1[0],p1[1],p1[2]), Vec(v[0],v[1],v[2])) 
    return cyl

def CSGcube(pCenter,size):
    s2 = [0.5*size[0],0.5*size[1],0.5*size[2]]
    p0 = VSub(pCenter,s2)
    p1 = VAdd(pCenter,s2)
    brick = OrthoBrick(Pnt(p0[0],p0[1],p0[2]),Pnt(p1[0],p1[1],p1[2]))
    return brick


#transform points
def TransformCrank(p, zOff, zRot):
    p2 = RotationMatrixZ(zRot) @ p
    pOff=[0,0,zOff]
    return VAdd(p2,pOff)

#cube only in XY-plane, z infinite
def CSGcubeXY(pCenter,sizeX,sizeY,ex,ey):
    #print("pCenter=",pCenter)
    pl1 = Plane(Pnt(pCenter[0]-0.5*sizeX*ex[0],pCenter[1]-0.5*sizeX*ex[1],0),Vec(-ex[0],-ex[1],-ex[2]))
    pl2 = Plane(Pnt(pCenter[0]+0.5*sizeX*ex[0],pCenter[1]+0.5*sizeX*ex[1],0),Vec( ex[0], ex[1], ex[2]))

    pl3 = Plane(Pnt(pCenter[0]-0.5*sizeY*ey[0],pCenter[1]-0.5*sizeY*ey[1],0),Vec(-ey[0],-ey[1],-ey[2]))
    pl4 = Plane(Pnt(pCenter[0]+0.5*sizeY*ey[0],pCenter[1]+0.5*sizeY*ey[1],0),Vec( ey[0], ey[1], ey[2]))

    return pl1*pl2*pl3*pl4
    

#create one crank face at certain z-offset and rotation; side=1: left, side=-1: right
def GetCrankFace(zOff, zRot, side=1):
    ex = RotationMatrixZ(zRot) @ [1,0,0]
    ey = RotationMatrixZ(zRot) @ [0,1,0]
    #print("zOff=",zOff, "zRot=", zRot, "side=", side,"ex=", ex)
    pLeft = [0,0,zOff]
    pRight = [0,0,zOff+dk]
    pMid = [0,0,zOff+0.5*dk]

    pcLeft=VAdd(pLeft,lk*ex)
    pcRight=VAdd(pRight,lk*ex)
    f=0.5**0.5
    cyl1pl = Plane(Pnt(pcLeft[0],pcLeft[1],pcLeft[2]+0.5*dk-side*dk),Vec(f*ex[0],f*ex[1],f*ex[2]-side*f))        
    cyl1 = Cylinder(Pnt(pcLeft[0],pcLeft[1],pcLeft[2]-1), Pnt(pcRight[0],pcRight[1],pcRight[2]+1), 0.5*bk)*cyl1pl

    #cone2 = Cylinder(Pnt(pcLeft[0],pcLeft[1],pcLeft[2]-1), Pnt(pcRight[0],pcRight[1],pcRight[2]+1), lk+l4)
    cone2 = Cone(Pnt(pcLeft[0],pcLeft[1],pcLeft[2]-side*dBevel+0.5*dk), Pnt(pcLeft[0],pcLeft[1],pcLeft[2]+side*dBevel+0.5*dk), lk+l4-1.5*dBevel, lk+l4-0.5*dBevel)
    cube1 = CSGcubeXY(VAdd(pMid,0.49*l3*ex),1.02*l3,bk,ex,ey) #make l3 a little longer, to avoid bad edges
    cube2 = CSGcubeXY(VAdd(pMid,-0.5*l4*ex),1.0*l4,bm,ex,ey)*cone2

    pc3a = VAdd(pLeft,0.*l3*ex+(0.5*bk+0.4*l3)*ey)
    cyl3a = Cylinder(Pnt(pc3a[0],pc3a[1],pc3a[2]-1), Pnt(pc3a[0],pc3a[1],pc3a[2]+1), 0.42*l3)
    pc3b = VAdd(pLeft,0.*l3*ex+(-0.5*bk-0.4*l3)*ey)
    cyl3b = Cylinder(Pnt(pc3b[0],pc3b[1],pc3b[2]-1), Pnt(pc3b[0],pc3b[1],pc3b[2]+1), 0.42*l3)
    #cube3a = (CSGcubeXY(VAdd(pMid,0.26*l3*ex+(0.5*bk+0.26*l3)*ey),0.5*l3,0.5*l3,ex,ey)-cyl3a)
    
    return ((cube1+cube2+cyl1)-(cyl3a+cyl3b))*Plane(Pnt(0,0,pLeft[2]),Vec(0,0,-1))*Plane(Pnt(0,0,pRight[2]),Vec(0,0,1))
    #return (cube1+cube2+cyl1)*Plane(Pnt(0,0,pLeft[2]),Vec(0,0,-1))*Plane(Pnt(0,0,pRight[2]),Vec(0,0,1))

#generate one crank, rotated around z-axis in radiant
def GenerateCrank(zOff, zRot):
    pL = TransformCrank(pLB,zOff, zRot)
    p0 = TransformCrank(p0B,zOff, zRot)
    p1 = TransformCrank(p1B,zOff, zRot)

    p21 = TransformCrank(p21B,zOff, zRot)
    p31 = TransformCrank(p31B,zOff, zRot)
    p41 = TransformCrank(p41B,zOff, zRot)
    p51 = TransformCrank(p51B,zOff, zRot)

    p6 = TransformCrank(p6B,zOff, zRot)
    p7 = TransformCrank(p7B,zOff, zRot)
    p8 = TransformCrank(p8B,zOff, zRot)
    
    crank0 = CSGcylinder(pL,[p0[0],p0[1],p0[2]+eps],r0)
    crank1 = CSGcylinder(p0,[p1[0],p1[1],p1[2]+eps],rdb0)

    #conrod bearing:
    crank3 = CSGcylinder([p21[0],p21[1],p21[2]-eps],p31,rdb1)
    crank7 = CSGcylinder(p31,p41,r1)
    crank8 = CSGcylinder(p41,[p51[0],p51[1],p51[2]+eps],rdb1)
    
    crank9 = CSGcylinder([p6[0],p6[1],p6[2]-eps],p7,rdb0)
    crank10 = CSGcylinder([p7[0],p7[1],p7[2]-eps],p8,r0)

    #return crank0+crank1+crank3+crank4+crank5+crank6+crank7+crank8+crank4b+crank5b+crank6b+crank9+crank10
    if zOff==0:#add first shaft
        crank1 = crank1+crank0
    return crank1+GetCrankFace(db+zOff,zRot,1)+crank3+crank7+crank8+GetCrankFace(db+2*db+dk+b1+zOff,zRot,-1)+crank10+crank9

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#choose configuration for crankshaft:
#crankConfig = [0] #1-piston
#crankConfig = [np.pi/2] #1-piston
#crankConfig = [0,np.pi] #2-piston
#crankConfig = [0,np.pi*2./3.,2.*np.pi*2./3.] #3-piston
#crankConfig = [0,np.pi,np.pi,0] #4-piston
crankConfig = [0,np.pi*2./3.,2.*np.pi*2./3.,2.*np.pi*2./3.,np.pi*2./3.,0] #6-piston

nPistons = len(crankConfig)
crank = GenerateCrank(0, crankConfig[0])

zPos = lTotal
for i in range(len(crankConfig)-1):
    angle = crankConfig[i+1]
    crank += GenerateCrank(zPos, angle)
    zPos += lTotal

# crank = (GenerateCrank(0, 0) + GenerateCrank(lTotal, np.pi*2./3.) + GenerateCrank(2*lTotal, np.pi*2.*2./3.)+
#           GenerateCrank(3*lTotal, np.pi*2.*2./3.) + GenerateCrank(4*lTotal, np.pi*2./3.))

geoCrank = CSGeometry()
geoCrank.Add(crank)

meshSize = 0.0075
mesh = Mesh( geoCrank.GenerateMesh(maxh=meshSize+0.001*0))
Draw(mesh)