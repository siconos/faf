#!/usr/bin/env python

import sys
sys.path.append('..')

from sympy import *
from funcodegen import funcodegen
from sympy.utilities.lambdify import lambdify
import numpy as np

import unittest

class Misc(unittest.TestCase):

    def test_1(self):
        x = Symbol('x', real=True)
        expr = x
        assert funcodegen('f', expr) == \
"""\
void f(
    double x,
    double *result)
{
    result[0] = x;
}
"""
    
    def test_2(self):
        x = Symbol('x', real=True)
        expr = sqrt(x)
        print '[{0}]'.format(funcodegen('f', expr))
        assert funcodegen('f', expr) == \
"""\
void f(
    double x,
    double *result)
{
    double x1 = 0.;
    x1 = sqrt(x);
    result[0] = x1;
}
"""

    def test_3(self):
        x = Symbol('x', real=True)
        expr = sqrt(x)*sqrt(x)*sqrt(x)
        assert funcodegen('f', expr) == \
"""\
void f(
    double x,
    double *result)
{
    double x1 = 0.;
    x1 = sqrt(x);
    result[0] = x*x1;
}
"""

    def test_4(self):
        x = Symbol('x', real=True)
        expr = 1/x
        assert funcodegen('f', expr) == \
"""\
void f(
    double x,
    double *result)
{
    double x1 = 0.;
    x1 = 1.0/x;
    result[0] = x1;
}
"""

    def test_5(self):
        x = Symbol('x', real=True)
        expr = sqrt(x*x*x)
        assert funcodegen('f', expr) == \
"""\
void f(
    double x,
    double *result)
{
    double x1 = 0.;
    double x2 = 0.;
    x1 = x*x*x;
    x2 = sqrt(x1);
    result[0] = x2;
}
"""

    def test_6(self):
        x = Symbol('x', real=True)
        y = Symbol('y', real=True)
        expr = x/sqrt(x*x+y*y)
        print funcodegen('f', expr)
        assert funcodegen('f', expr) == \
"""\
void f(
    double x,
    double y,
    double *result)
{
    double x1 = 0.;
    double x2 = 0.;
    double x3 = 0.;
    double x4 = 0.;
    x1 = x*x;
    x2 = y*y;
    x3 = sqrt(x1 + x2);
    x4 = 1.0/x3;
    result[0] = 1.0*x*x4;
}
"""

if __name__ == '__main__':
    unittest.main()
