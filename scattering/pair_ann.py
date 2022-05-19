
# PAIR ANNIHILATION

import sys
import math
import cmath
import numpy as np
import fourvector as fv

COUPLING = fv.COUPLING

# p1 - incoming electron
# p2 - incoming positron
# k1 - outgoing photon
# k2 - outgoing photon
def Annihilation(p1, p2, k1, k2, toPrint = False, coupling = COUPLING):
	# numeric factor due to averaging and coupling constant (power = 2 * vertices)
	factor = 0.25 * coupling**4
	res = 0.0
	for h1 in fv.Helicity:
		for h2 in fv.Helicity:
			for h3 in fv.Helicity:
				for h4 in fv.Helicity:
					u1    = fv.FourVector.InitVector(p1, fv.WhichVector.u,       h1)	# electron
					vbar2 = fv.FourVector.InitVector(p2, fv.WhichVector.vbar,    h2)	# positron
					ph1   = fv.FourVector.InitVector(k1, fv.WhichVector.epsilon, h3)	# photon
					ph2   = fv.FourVector.InitVector(k2, fv.WhichVector.epsilon, h4)	# photon

					numerator1 = (((vbar2 * ~ph2) * (~(p1 - k1) + p1.mass)) * ~ph1) * u1
					numerator2 = (((vbar2 * ~ph1) * (~(p1 - k2) + p1.mass)) * ~ph2) * u1

					temp_res = (numerator1 / ((p1 - k1).Squared() - p1.mass**2) +
							    numerator2 / ((p1 - k2).Squared() - p1.mass**2))

					if toPrint:
						print(f"h1 = {h1.value}, u1 = {u1}")
						print(f"   Dirac eq for u1      : {(~p1 * u1)} = {u1 * p1.mass}")
						print(f"h2 = {h2.value}, vbar2 = {vbar2}")
						print(f"   Dirac eq for vbar2   : {(vbar2 * ~p2)} = {vbar2 * (-1)*p2.mass}")
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

px = 1.0
py = 0.0
pz = 0.1

kx = 1 / math.sqrt(2)
ky = math.sqrt(mass**2 + px**2 + py**2 - kx**2) # due to momentum conservation
kz = pz

p1 = fv.FourVector([ px,  py, pz], mass)    # electron
p2 = fv.FourVector([-px, -py, pz], mass) 	# positron
k1 = fv.FourVector([ kx,  ky, kz])			# photon
k2 = fv.FourVector([-kx, -ky, kz])			# photon

print(f"PAIR ANNIHILATION RES = {fv.format_num(Annihilation(p1, p2, k1, k2, toPrint), precision = 15)}")
print(f"Mass = {mass}")
print()

s = (p1 + p2).Squared()
t = (p1 - k1).Squared()
u = (p1 - k2).Squared()

print(f"Momentum conservation: {p1 + p2} = {k1 + k2}")
print()
print(f"2*m^2 = {fv.format_num(2 * mass**2)}")
print(f"s+t+u = {fv.format_num(u+s+t)}")
print()
print(f"s     = {fv.format_num(s, precision = 15)}")
print(f"t     = {fv.format_num(t, precision = 15)}")
print(f"u     = {fv.format_num(u, precision = 15)}")
