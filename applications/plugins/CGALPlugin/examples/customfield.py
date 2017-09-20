import pyximport; pyximport.install()
import dfgeom
from dfgeom import *
import math

import accordion
import primitives


class Difference2(PythonShape):
        def __init__(self, a, b):
            self.a = a
            self.b = b

        def evalField(self, p):
                return  max(-self.a.evalField(p), self.b.evalField(p))


r = 0.
dd = -0.0
s1 = Sphere([1.0,0.0,0.0],1.0)
s2 = Sphere([dd,0.0,0.0],1.3)
s3 = Sphere([-1.0,0.0,0.0],1.0)
u = Union(s1,s3)
d = Difference(u, s2)
u2 = Sphere([dd,1.0,0],r)
u3 = Sphere([dd,-1.0,0],r)
u4 = Sphere([dd,0,1.0],r)
u5 = Sphere([dd,0,-1.0],r)
u3= Union(u2, u3)
u4= Union(u4, u5)
u5= Union(u3,u4)
d = Difference(u5, d)


#def getCythonFunction():
#    return dfgeom.getCythonRawFunction(d)

ellipsoid1=primitives.Ellipsoid("+",1.0,1.0,5.0,0.0,0.0,primitives.Point(0.0,0.0,5.0))
ellipsoid2=primitives.Ellipsoid("+",0.5,0.5,4.5,0.0,0.0,primitives.Point(0.0,0.0,5.0))
ellipsoidcreux=primitives.Difference(ellipsoid1,ellipsoid2)

cylinder1=primitives.Cylinder("+",1.0,1.0,5.0,0.0,0.0,primitives.Point(0.0,0.0,5.0))
cylinder2=primitives.Cylinder("+",0.5,0.5,4.5,0.0,0.0,primitives.Point(0.0,0.0,5.0))
cylindercreux=primitives.Difference(cylinder1,cylinder2)
#shape1=primitives.Union(cylindercreux,ellipsoidcreux)


#formepleine=primitives.Union(ellipsoid1,cylinder1)
#formecreuse=primitives.Union(ellipsoid2,cylinder2)
#shape2=primitives.Difference(formepleine,formecreuse)

#(shape,shapeInt)=accordion.accordionUniform(7.0,2.0,0.5, "ellipsoid",1,5.0,5.0,0.0)
#(shape,shapeInt)=accordion.accordionFreeDimension(7.0,3.0,1.0,[[5.0,"ellipsoid",5.0,5.0,5.0]])

shape=primitives.Difference(cylinder1,cylinder2)

def evalField(x, y, z):
    point=primitives.Point(x,y,z)
    return shape.eval(point) #dfgeom.evalField(d, x, y, z)
