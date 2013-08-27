#!/usr/bin/env python

from math import cos, sin, pi
from numpy.linalg import norm

from Siconos.Mechanics.ContactDetection.Bullet import btQuaternion

from Siconos.Mechanics.ContactDetection.Bullet.BulletWrap import __mul__ as mul


theta1 = pi / 2
a1 = 1
b1 = 0
c1 = 0
n1 = sin(theta1 / 2) / norm((a1, b1, c1))

r1 = btQuaternion(a1 * n1, b1 * n1, c1 * n1, cos(theta1 / 2))


def s_v(v):
    return ' '.join('{0}'.format(iv) for iv in v)

alpha = pi / 6

with open('input.dat', 'w') as f:
    # ground
    f.write('1 -1 0 0  0 -.5 1 0 0 0 0    0 0 0 0 0\n')

    # one object
    theta = 0
    a = 0
    b = 0
    c = 1
    n = sin(theta / 2) / norm((a, b, c))
    q = (0, 0, 30)

    r = btQuaternion(a * n, b * n, c * n, cos(theta / 2))
    r = mul(r, r1)
    r.normalize()

    o = (r.w(), r.x(), r.y(), r.z())

    v = (10, 0, 0, 5, 0, 0)

    f.write('0 0 1 {0} {1} {2}\n'.format(s_v(q), s_v(o), s_v(v)))
