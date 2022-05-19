
# DOUBLE COMPTON

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
# k5 - outgoing photon
def DoubleCompton(k1, k2, k3, k4, k5, toPrint = False, coupling = COUPLING):
	# numeric factor due to averaging and coupling constant (power = 2 * vertices)
	factor = 0.25 * coupling**6
	res = 0.0
	for h1 in fv.Helicity:
		for h2 in fv.Helicity:
			for h3 in fv.Helicity:
				for h4 in fv.Helicity:
					for h5 in fv.Helicity:
						u1    = fv.FourVector.InitVector(k1, fv.WhichVector.u,       h1) # electron
						ph2   = fv.FourVector.InitVector(k2, fv.WhichVector.epsilon, h2) # photon
						ubar3 = fv.FourVector.InitVector(k3, fv.WhichVector.ubar,    h3) # electron
						ph4   = fv.FourVector.InitVector(k4, fv.WhichVector.epsilon, h4) # photon
						ph5   = fv.FourVector.InitVector(k5, fv.WhichVector.epsilon, h5) # photon

						numerator1 = (((((ubar3 * ~ph5) * (~(k3 + k5) + k3.mass)) * ~ph4) * (~(k1 + k2) + k1.mass)) * ~ph2) * u1
						numerator2 = (((((ubar3 * ~ph5) * (~(k3 + k5) + k3.mass)) * ~ph2) * (~(k1 - k4) + k1.mass)) * ~ph4) * u1
						numerator3 = (((((ubar3 * ~ph2) * (~(k3 - k2) + k3.mass)) * ~ph4) * (~(k1 - k5) + k1.mass)) * ~ph5) * u1

						numerator4 = (((((ubar3 * ~ph4) * (~(k3 + k4) + k3.mass)) * ~ph5) * (~(k1 + k2) + k1.mass)) * ~ph2) * u1
						numerator5 = (((((ubar3 * ~ph4) * (~(k3 + k4) + k3.mass)) * ~ph2) * (~(k1 - k5) + k1.mass)) * ~ph5) * u1
						numerator6 = (((((ubar3 * ~ph2) * (~(k3 - k2) + k3.mass)) * ~ph5) * (~(k1 - k4) + k1.mass)) * ~ph4) * u1

						temp_res = (numerator1 / ( ((k3 + k5).Squared() - k3.mass**2) * ((k1 + k2).Squared() - k1.mass**2) ) +
									numerator4 / ( ((k3 + k4).Squared() - k3.mass**2) * ((k1 + k2).Squared() - k1.mass**2) ) +
								    numerator2 / ( ((k3 + k5).Squared() - k3.mass**2) * ((k1 - k4).Squared() - k1.mass**2) ) +
								    numerator5 / ( ((k3 + k4).Squared() - k3.mass**2) * ((k1 - k5).Squared() - k1.mass**2) ) +
								    numerator3 / ( ((k3 - k2).Squared() - k3.mass**2) * ((k1 - k5).Squared() - k1.mass**2) ) +
								    numerator6 / ( ((k3 - k2).Squared() - k3.mass**2) * ((k1 - k4).Squared() - k1.mass**2) ) )

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

# incoming electron
k1x = 1.0
k1y = 0.0
k1z = 0.1
k1  = fv.FourVector([k1x, k1y, k1z], mass)

# incoming photon
k2x = -k1x
k2y = -k1y
k2z =  k1z
k2  = fv.FourVector([k2x, k2y, k2z])

# outgoing electron
k3x = -(k1x + k2x)
k3y = -(k1y + k2y)
k3z =  (k1z + k2z) / 3
k3  = fv.FourVector([k3x, k3y, k3z], mass)

# outgoing photon
k4x = 0.25
k4y = math.sqrt( ((k1[0]+k2[0]-k3[0])**2 / 4) - k3z**2 - k4x**2)
k4z = k3z
k4  = fv.FourVector([k4x, k4y, k4z])

# outgoing photon
k5x = -k4x
k5y = -k4y
k5z =  k3z
k5  = fv.FourVector([k5x, k5y, k5z])

print(f"DOUBLE COMPTON RES = {fv.format_num(DoubleCompton(k1, k2, k3, k4, k5, toPrint), precision = 14)}")
print(f"Mass = {mass}")
print()
print(f"Momentum conservation: {k1 + k2} = {k3 + k4 + k5}")
print()
print(f"2*m^2 = {fv.format_num(2 * mass**2)}")
print()
print(f"k1*k2 = {fv.format_num(k1*k2)}")
print(f"k1*k3 = {fv.format_num(k1*k3)}")
print(f"k1*k4 = {fv.format_num(k1*k4)}")
print(f"k2*k3 = {fv.format_num(k2*k3)}")
print(f"k2*k4 = {fv.format_num(k2*k4)}")
print(f"k3*k4 = {fv.format_num(k3*k4)}")
