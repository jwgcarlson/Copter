#!/usr/bin/python

from numpy import *
from pycopter import *

x = linspace(0, 4*pi, 20)
y = sin(x)

f1 = LinearSpline(x, y)
f2 = ShiftedLinearSpline(x, y)
f3 = CubicSpline(x, y)

x = linspace(0, 4*pi, 2000)
y = sin(x)
y1 = array([f1(xx) for xx in x])
y2 = array([f2(xx) for xx in x])
y3 = array([f3(xx) for xx in x])

chi1 = sum((y - y1)**2)
chi2 = sum((y - y2)**2)
chi3 = sum((y - y3)**2)
print "chi1 = %g, chi2 = %g, chi3 = %g" % (chi1,chi2,chi3)

try:
    import pylab
    pylab.plot(x,y, x,y1, x,y2, x,y3)
    pylab.show()
except:
    pass
