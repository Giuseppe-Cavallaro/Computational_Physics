
# BHABHA

import sys
import math
import cmath
import numpy as np
import fourvector as fv

COUPLING = fv.COUPLING

# k1 - incoming electron
# k2 - incoming positron
# k3 - outgoing electron
# k4 - outgoing positron
def Bhabha(k1, k2, k3, k4, toPrint = False, coupling = COUPLING):
	# numeric factor due to averaging and coupling constant (power = 2 * vertices)
	factor = 0.25 * coupling**4
	res = 0.0
	for h1 in fv.Helicity:
		for h2 in fv.Helicity:
			for h3 in fv.Helicity:
				for h4 in fv.Helicity:
					u1    = fv.FourVector.InitVector(k1, fv.WhichVector.u,    h1) # electron
					vbar2 = fv.FourVector.InitVector(k2, fv.WhichVector.vbar, h2) # positron
					ubar3 = fv.FourVector.InitVector(k3, fv.WhichVector.ubar, h3) # electron
					v4    = fv.FourVector.InitVector(k4, fv.WhichVector.v,    h4) # positron

					numerator1   = fv.VertexMul(ubar3, v4) * fv.VertexMul(vbar2, u1, True)
					numerator2   = fv.VertexMul(ubar3, u1) * fv.VertexMul(vbar2, v4, True)

					# bhabha scattering has a relative minus between amplitudes
					temp_res = (numerator1 / (k1 + k2).Squared() - numerator2 / (k1 - k3).Squared())

					if toPrint:
						print(f"h1 = {h1.value}, u1 = {u1}")
						print(f"   Dirac eq for u1   : {(~k1 * u1)} = {u1 * k1.mass}")
						print(f"h2 = {h2.value}, vbar2 = {vbar2}")
						print(f"   Dirac eq for vbar2: {(vbar2 * ~k2)} = {vbar2 * (-1)*k2.mass}")
						print(f"h3 = {h3.value}, ubar3 = {ubar3}")
						print(f"   Dirac eq for ubar3: {(ubar3 * ~k3)} = {ubar3 * k3.mass}")
						print(f"h4 = {h4.value}, v4 = {v4}")
						print(f"   Dirac eq for v4   : {(~k4 * v4)} = {v4 * (-1)*k4.mass}")
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
k2 = fv.FourVector([-px, -py, pz], mass) 	# positron
k3 = fv.FourVector([ kx,  ky, kz], mass)	# electron
k4 = fv.FourVector([-kx, -ky, kz], mass)	# positron

print(f"BHABHA RES = {fv.format_num(Bhabha(k1, k2, k3, k4, toPrint), precision = 14)}")
print(f"Mass = {mass}")
print()

s = (k1 + k2).Squared()
t = (k1 - k3).Squared()
u = (k1 - k4).Squared()

print(f"Momentum conservation: {k1 + k2} = {k3 + k4}")
print()
print(f"4*m^2 = {fv.format_num(4 * mass**2)}")
print(f"s+t+u = {fv.format_num(u+s+t)}")
print()
print(f"s     = {fv.format_num(s, precision = 14)}")
print(f"t     = {fv.format_num(t, precision = 14)}")
print(f"u     = {fv.format_num(u, precision = 14)}")
