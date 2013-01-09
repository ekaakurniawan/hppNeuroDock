from LFSR import *

print "########################################################## Test LFSR ###"
lfsr = LFSR(lfsr = 1070, bit_len = 16)
print "### LFSR - Generate"
for x in xrange(1, 10):
    print lfsr.generate()
print "### LFSR - Random 0 to 1"
for x in xrange(1, 10):
    print lfsr.zero_to_one()
print "### LFSR - Random 0 to 2pi"
for x in xrange(1, 10):
    print lfsr.zero_to_2pi()
print "### LFSR - Random Sign"
for x in xrange(1, 10):
    print lfsr.sign()

from Axis3 import *
print ""
print "############################################################## Axis3 ###"
axis = Axis3()
print "### Axis3 - Init"
print axis
axis = Axis3(1.2, 2.3, 3.4)
print axis
print "### Axis3 - Normalize"
axis.normalize()
print axis

from Quaternion import *
print ""
print "######################################################### Quaternion ###"
q = Quaternion()
rng = LFSR(lfsr = 1070, bit_len = 16)

print "### Quaternion - Init"
print q
q = Quaternion(1.2, 2.3, 3.4, 4.5)
print q

print "### Quaternion - Copy"
q1 = q.copy()
q1.a = 100.002
q1.b = 200.003
q1.c = 300.004
q1.d = 400.005
print "q  = %s" % q
print "q1 = %s" % q1

print "### Quaternion - Uniform"
for x in xrange(1, 10):
    q.uniform(rng)
    angle, axis = q.get_angle_axis()
    print "%s, angle=%.4f, %s" % (q, angle, axis)

print "### Quaternion - Multiplication"
q.uniform(rng)
q1.uniform(rng)
print "(%s) x (%s)" % (q, q1)
print q * q1
q *= q1
print q

print "### Quaternion - Multiplication"
print "q = %s" % q
print "magnitude: abs(q)  = %s" % q.magnitude()
print "q1 = %s" % q1
print "magnitude: abs(q1) = %s" % q1.magnitude()
print "### Quaternion - Normalize"
print "q =            %s" % q
print "magnitude: abs(q)  = %s" % q.magnitude()
q.normalize()
print "normalized q = %s" % q
print "magnitude: abs(q)  = %s" % q.magnitude()

q.uniform(rng)
print "q =            %s" % q
print "magnitude: abs(q)  = %s" % q.magnitude()
q.normalize()
print "normalized q = %s" % q
print "magnitude: abs(q)  = %s" % q.magnitude()

q = Quaternion(a = 1000.0002, b = 2.03, c = 0.04, d = 40004.5)
print "q =            %s" % q
print "magnitude: abs(q)  = %s" % q.magnitude()
q.normalize()
print "normalized q = %s" % q
print "magnitude: abs(q)  = %s" % q.magnitude()

print "### Quaternion - Conjugate"
q.uniform(rng)
print "q             = %s" % q
q.conjugate()
print "q.conjugate() = %s" % q

print "### Quaternion - Identity"
q.identity()
print "q.identity  = %s" % q
print "q.magnitude = %s" % q.magnitude()

print "### Quaternion - Convert Quaternion to Angle-Axis"
q = Quaternion(a = 0.707, b = -0.240, c = -0.665, d = 0.000)
print q
print "Angle = %s, %s" % q.get_angle_axis()

print "### Quaternion - Convert Angle-Axis to Quaternion"
q.identity()
angle = 1.57
axis = Axis3(x = -0.340, y = -0.940, z = 0.000)
print "Angle = %s, %s" % (angle, axis)
q.set_angle_axis(angle, axis)
print q

