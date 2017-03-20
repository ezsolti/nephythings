# -*- coding: utf-8 -*-
"""
plane equation calculator from three points
implemented theory: http://www.math.cornell.edu/~froh/231f08e1a.pdf
ax+by+cz=d

The script is originally made for mcnp input deck creation, hence the "-" and "+" 
sign for the plane is related to the cell definition.
"""

import numpy as np
#three known points on the plane 
p1=np.array([-4.1,60,-0.25])
p2=np.array([-11.7,50,-0.25])
p3=np.array([-4.1,60,0.25])

#test point in space which is on the "preferred side" of the plane
p0=np.array([0,50,0])
v1=p1-p2
v2=p3-p2

ov=np.cross(v1,v2) #direction of normal vector

a=ov[0]
b=ov[1]
c=ov[2]

d=a*p1[0]+b*p1[1]+c*p1[2]

print("A: %f, B: %f, C: %f, D: %f" %(a,b,c,d))


p0Dir=a*p0[0]+b*p0[1]+c*p0[2]-d
if p0Dir>0:
    print("The surface has + sign")
elif p0Dir<0:
    print("The surface has - sign")
else:
    print("This point is on the surface")
