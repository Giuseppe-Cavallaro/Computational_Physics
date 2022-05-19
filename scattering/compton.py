
# COMPTON

import sys
import math
import cmath
import numpy as np
import fourvector as fv

COUPLING = fv.COUPLING

# k1 - incoming electron
# k2 - incoming photon
# k3 - outgoing electron
# k4 - outgoing photon
def Compton(k1, k2, k3, k4, toPrint = False, coupling = COUPLING):
	# numeric factor due to averaging and coupling constant (power = 2 * vertices)
	factor = 0.25 * coupling**4
	res = 0.0
	for h1 in fv.Helicity:
		for h2 in fv.Helicity:
			for h3 in fv.Helicity:
				for h4 in fv.Helicity:
					ubar3 = fv.FourVector.InitVector(k3, fv.WhichVector.ubar,    h3)
					u1    = fv.FourVector.InitVector(k1, fv.WhichVector.u,       h1)
					ph2   = fv.FourVector.InitVector(k2, fv.WhichVector.epsilon, h2)
					ph4   = fv.FourVector.InitVector(k4, fv.WhichVector.epsilon, h4)

					numerator1 = (((ubar3 * ~ph4) * (~(k3 + k4) + k3.mass)) * ~ph2) * u1
					numerator2 = (((ubar3 * ~ph2) * (~(k3 - k2) + k3.mass)) * ~ph4) * u1

					temp_res = (numerator1 / ((k3 + k4).Squared() - k3.mass**2) +
							    numerator2 / ((k3 - k2).Squared() - k3.mass**2))

					if toPrint:
						print(f"h1 = {h1.value}, u1 = {u1}")
						print(f"   Dirac eq for u1   : {(~k1 * u1)} = {u1 * k1.mass}")
						print(f"h3 = {h3.value}, ubar3 = {ubar3}")
						print(f"   Dirac eq for ubar3: {(ubar3 * ~k3)} = {ubar3 * k3.mass}")
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
ky = math.sqrt(px**2 + py**2 - kx**2) # due to momentum conservation
kz = pz

k1 = fv.FourVector([ px,  py, pz], mass)    # electron
k2 = fv.FourVector([-px, -py, pz]) 			# photon
k3 = fv.FourVector([ kx,  ky, kz], mass)	# electron
k4 = fv.FourVector([-kx, -ky, kz])			# photon

print(f"COMPTON RES = {fv.format_num(Compton(k1, k2, k3, k4, toPrint), precision = 14)}")
print(f"Mass = {mass}")
print()

s = (k1 + k2).Squared()
t = (k1 - k3).Squared()
u = (k1 - k4).Squared()

print(f"Momentum conservation: {k1 + k2} = {k3 + k4}")
print()
print(f"2*m^2 = {fv.format_num(2 * mass**2)}")
print(f"s+t+u = {fv.format_num(u+s+t)}")
print()
print(f"s     = {fv.format_num(s, precision=14)}")
print(f"t     = {fv.format_num(t, precision=14)}")
print(f"u     = {fv.format_num(u, precision=14)}")
