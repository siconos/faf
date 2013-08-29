#!/usr/bin/env python

import random

from math import sqrt, cos, sin, pi
from numpy.linalg import norm

from Siconos.Mechanics.ContactDetection.Bullet import btQuaternion

from Siconos.Mechanics.ContactDetection.Bullet.BulletWrap import __mul__ as mul


theta1 = pi/2
a1 = 1
b1 = 0
c1 = 0
n1 = sin(theta1/2)/norm((a1,b1,c1))

r1 = btQuaternion(a1*n1, b1*n1, c1*n1, cos(theta1/2))

def s_v(v):
    return ' '.join('{0}'.format(iv) for iv in v)

alpha = pi/6

with open('input.dat','w') as f:
    #    f.write('0 0 10 50 0 20 1 0 0 0 -100. 0 0 10 10 10\n')
    f.write('1 -1 0 0  0 -.5 1 0 0 0 0    0 0 0 0 0\n')
    for k in range(0,2):

        for i in range(0,12):

            theta = (i+3 + (k%2) * 0.5) * pi / 6
            a = 0
            b = 0
            c = 1
            n = sin(theta / 2) / norm((a, b, c))

            r = btQuaternion(a*n, b*n, c*n, cos(theta/2))
            r = mul(r, r1)

            r.normalize()

            q = (23.48*cos(alpha*(i + (k%2)* 0.5)), 23.48*sin(alpha*(i + (k%2)* 0.5)), (2.348 + 1.) * k + 2.360 /2. + 2)
            o = (r.w(), r.x(), r.y(), r.z())

            f.write('2 0 1 {0} {1} 0 0 0 0 0 0\n'.format(s_v(q),s_v(o)))


