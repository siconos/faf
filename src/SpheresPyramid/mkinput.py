#!/usr/bin/env python

from math import cos, sin, pi, sqrt
from numpy.linalg import norm

from Siconos.Mechanics.ContactDetection.Bullet import btQuaternion

from Siconos.Mechanics.ContactDetection.Bullet.BulletWrap import __mul__ as mul

from random import random


N = 3

theta1 = pi / 2
a1 = 1
b1 = 0
c1 = 0
n1 = sin(theta1 / 2) / norm((a1, b1, c1))

r1 = btQuaternion(a1 * n1, b1 * n1, c1 * n1, cos(theta1 / 2))



def s_v(v):
    return ' '.join('{0}'.format(iv) for iv in v)

alpha = pi / 2.

with open('input.dat', 'w') as f:
    # ground
    f.write('1 -1 0 0  0 -.5 1 0 0 0 0 0 0 0 0 0\n')

    k = N

    for i in range(k-1,-1,-1):
        for x in range(0,i+1):
            for y in range(0,i+1):
                theta = 0
                a = 0
                b = 0
                c = 1
                n = sin(theta / 2) / norm((a, b, c))
                q1 = ((i+1 - 2 * x)/10. , (i+1 - 2 * y)/10. , ((k-i-1)*sqrt(2) + 1)/10.)

                r = btQuaternion(a * n, b * n, c * n, cos(theta / 2))
                r = mul(r, r1)
                r.normalize()

                o = (r.w(), r.x(), r.y(), r.z())
                v1 = (0, 0, 0, 0, 0, 0)

                f.write('0 0 1 {0} {1} {2}\n'.format(s_v(q1), s_v(o), s_v(v1)))
