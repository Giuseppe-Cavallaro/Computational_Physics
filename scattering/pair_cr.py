
# PAIR CREATION

import sys
import math
import cmath
import numpy as np
import fourvector as fv

COUPLING = fv.COUPLING

# k1 - incoming photon
# k2 - incoming photon
# p1 - outgoing electron
# p2 - outgoing positron
def Creation(p1, p2, k1, k2, toPrint = False, coupling = COUPLING):
	# numeric factor due to averaging and coupling constant (power = 2 * vertices)
	factor = 0.25 * coupling**4
	res = 0.0
	for h1 in fv.Helicity:
		for h2 in fv.Helicity:
			for h3 in fv.Helicity:
				for h4 in fv.Helicity:
					ubar1 = fv.FourVector.InitVector(p1, fv.WhichVector.ubar,    h1)	# electron
					v2    = fv.FourVector.InitVector(p2, fv.WhichVector.v,       h2)	# positron
					ph1   = fv.FourVector.InitVector(k1, fv.WhichVector.epsilon, h3)	# photon
					ph2   = fv.FourVector.InitVector(k2, fv.WhichVector.epsilon, h4)	# photon

					numerator1 = (((ubar1 * ~ph1) * (~(p1 - k1) + p1.mass)) * ~ph2) * v2
					numerator2 = (((ubar1 * ~ph2) * (~(p1 - k2) + p1.mass)) * ~ph1) * v2

					temp_res = (numerator1 / ((p1 - k1).Squared() - p1.mass**2) +
							    numerator2 / ((p1 - k2).Squared() - p1.mass**2))

					if toPrint:
						print(f"h1 = {h1.value}, ubar1 = {ubar1}")
						print(f"   Dirac eq for ubar1 : {(ubar1 * ~p1)} = {ubar1 * p1.mass}")
						print(f"h2 = {h2.value}, v2    = {v2}")
						print(f"   Dirac eq for v2    : {(~p2 * v2)} = {v2 * (-1)*p2.mass}")
						print(f"temp res = {fv.format_num(temp_res)}")
						print()

					res += (temp_res * temp_res.conjugate())

	return factor * res

#-------------------------------------------------------------------------------------------------------------------------------------------------

argv = sys.argv
mass = 0.0
try:
	argv.remove("-p")
	toPrint = True
except Exception:
	toPrint = False
if len(argv) > 1:
	mass = float(argv[1])

kx = 5 / math.sqrt(2)
ky = 5 / math.sqrt(2)
kz = 0.1

px = 1.0
py = math.sqrt(kx**2 + ky**2 - mass**2 - px**2)	# due to momentum conservation
pz = kz

k1 = fv.FourVector([ kx,  ky, kz])			# photon
k2 = fv.FourVector([-kx, -ky, kz])			# photon
p1 = fv.FourVector([ px,  py, pz], mass)    # electron
p2 = fv.FourVector([-px, -py, pz], mass) 	# positron

print(f"PAIR CREATION RES = {fv.format_num(Creation(p1, p2, k1, k2, toPrint), precision = 14)}")
print(f"Mass = {mass}")
print()

s = (k1 + k2).Squared()
t = (k1 - p1).Squared()
u = (k1 - p2).Squared()

print(f"Momentum conservation: {p1 + p2} = {k1 + k2}")
print()
print(f"2*m^2 = {fv.format_num(2 * mass**2)}")
print(f"s+t+u = {fv.format_num(u+s+t)}")
print()
print(f"s     = {fv.format_num(s, precision = 14)}")
print(f"t     = {fv.format_num(t, precision = 14)}")
print(f"u     = {fv.format_num(u, precision = 14)}")
