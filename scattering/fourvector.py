import math
import cmath
import numpy as np
import warnings
import decimal
from enum import Enum

warnings.filterwarnings('error')

PRECISION = 3
FORMAT_SPACE = 25

COUPLING = 1

def format_num(n, space = 0, precision = PRECISION):
	if space == 0:
		space = ""
	r = float(n.real)
	i = float(n.imag)
	r_str = "{1:.{0}f}".format(precision, r).rstrip("0").rstrip(".")
	if r_str == "-0":
		r_str = "0"
	i_str = ""
	if i != 0.0:
		i_str = "{1:+.{0}f}".format(precision, i).rstrip("0").rstrip(".")
		if (i_str == "+0") or (i_str == "-0"):
			i_str = ""
		else:
			i_str = i_str[0] + "i" + i_str[1:]
	if r_str == "0" and i_str != "":
		if i_str[0] == "+":
			i_str = i_str[1:]
		res_str = "{1:>{0}}".format(space, i_str)
	else:
		res_str = "{1:>{0}}".format(space, r_str + i_str)
	return res_str

class Helicity(Enum):
	plus = 1
	minus = -1

class WhichVector(Enum):
	u = "u"
	ubar = "u_bar"
	v = "v"
	vbar = "v_bar"
	epsilon = "epsilon"

def check_u_ubar(V, mass = 0.0):
	space = FORMAT_SPACE
	print("CHECK U U_BAR, mass = {}".format(mass))
	v = FourVector(V, mass)
	u_plus     = FourVector.InitVector(v, WhichVector.u,    Helicity.plus)
	ubar_plus  = FourVector.InitVector(v, WhichVector.ubar, Helicity.plus)
	u_minus    = FourVector.InitVector(v, WhichVector.u,    Helicity.minus)
	ubar_minus = FourVector.InitVector(v, WhichVector.ubar, Helicity.minus)
	kslash = [[0.0, 0.0, 0.0, 0.0],
			  [0.0, 0.0, 0.0, 0.0],
			  [0.0, 0.0, 0.0, 0.0],
			  [0.0, 0.0, 0.0, 0.0]]
	for i in range(0, 4):
		for j in range(0, 4):
			kslash[i][j] = u_plus[i] * ubar_plus[j] + u_minus[i] * ubar_minus[j]

	for i in range(4):
		for j in range(4):
			print(format_num(kslash[i][j], space), end="")
		print()
	print()

	vtpz = v[0] + v[3]
	vtmz = v[0] - v[3]
	vp   = v[1] + 1j*v[2]
	vm   = v[1] - 1j*v[2]
	kslash2 = [[mass, 0.0, vtmz, -vm],
			   [0.0, mass, -vp, vtpz],
			   [vtpz, vm, mass, 0.0],
			   [vp, vtmz, 0.0, mass]]

	for i in range(4):
		for j in range(4):
			print(format_num(kslash2[i][j], space), end="")
		print()

def check_epstar_ep(k):
	space = FORMAT_SPACE
	print("CHECK EPSTAR_EP")
	epsilon_plus  = FourVector.InitVector(k, WhichVector.epsilon, Helicity.plus)
	epsilon_minus = FourVector.InitVector(k, WhichVector.epsilon, Helicity.minus)

	polsum = [[0.0, 0.0, 0.0, 0.0],
			  [0.0, 0.0, 0.0, 0.0],
			  [0.0, 0.0, 0.0, 0.0],
			  [0.0, 0.0, 0.0, 0.0]]

	for i in range(4):
		for j in range(4):
			polsum[i][j] = epsilon_plus[i].conjugate() * epsilon_plus[j] + epsilon_minus[i].conjugate() * epsilon_minus[j]

	for i in range(4):
		for j in range(4):
			print(format_num(polsum[i][j], space), end="")
		print()

	print()

	k_squared = k[0]**2 + k[1]**2 + k[2]**2
	for i in range(3):
		for j in range(3):
			polsum[i+1][j+1] = - k[i]*k[j] / k_squared
			if i == j:
				polsum[i+1][j+1] += 1

	for i in range(4):
		for j in range(4):
			print(format_num(polsum[i][j], space), end="")
		print()

class FourVector:
	def __init__(self, lst, mass = 0.0):
		try:
			self.v = lst.v
			self._mass = lst._mass
		except Exception:
			if len(lst) == 3:
				self.v = np.array([0.0] + lst)
				self.mass = mass
			else:
				self.v = np.array(lst)
				self._mass = None
		try:
			self.vectorType = lst.vectorType
			self.hel = lst.hel
		except Exception:
			self.vectorType = None
			self.hel = None

	def Squared(self):
		return self.v[0]*self.v[0].conjugate() - self.v[1]*self.v[1].conjugate() - self.v[2]*self.v[2].conjugate() - self.v[3]*self.v[3].conjugate()

	def threeSquared(self):
		return self.v[1]*self.v[1].conjugate() + self.v[2]*self.v[2].conjugate() + self.v[3]*self.v[3].conjugate()

	def mass(self, value):
		self._mass = value
		try:
			self.v[0] = math.sqrt(value**2 + self.threeSquared())
		except Warning:
			self.v[0] = cmath.sqrt(value**2 + self.threeSquared())
	mass = property(lambda x : x._mass, mass)

	@classmethod
	def InitVector(cls, input_vector, vectorType, hel, mass = 0.0):
		k = FourVector(input_vector, mass)
		if vectorType != WhichVector.epsilon:
			k.vectorType = vectorType
			k.hel = hel
			modk = math.sqrt(k.threeSquared())
			normchi = 1 / math.sqrt(2 * modk * (modk + k.v[3]))
			chip0 =   modk + k.v[3]
			chip1 =   k.v[1] + 1j * k.v[2]
			chim0 = - k.v[1] + 1j * k.v[2]
			chim1 = chip0

			normdm = math.sqrt(k.v[0] - modk) * normchi
			normdp = math.sqrt(k.v[0] + modk) * normchi

			k.v = np.array([0, 0, 0, 0], dtype=np.complex_)
			if k.isu:
				if hel == Helicity.plus:
					k.v[0] = normdm * chip0
					k.v[1] = normdm * chip1
					k.v[2] = normdp * chip0
					k.v[3] = normdp * chip1
				else:
					k.v[0] = normdp * chim0
					k.v[1] = normdp * chim1
					k.v[2] = normdm * chim0
					k.v[3] = normdm * chim1
			elif k.isv:
				if hel == Helicity.plus:
					k.v[0] = -normdp * chim0
					k.v[1] = -normdp * chim1
					k.v[2] =  normdm * chim0
					k.v[3] =  normdm * chim1
				else:
					k.v[0] =  normdm * chip0
					k.v[1] =  normdm * chip1
					k.v[2] = -normdp * chip0
					k.v[3] = -normdp * chip1

			if k.isbar:
				for i in range(4):
					k.v[i] = k.v[i].conjugate()
				tmp       = k.v[0]
				k.v[0] = k.v[2]
				k.v[2] = tmp
				tmp       = k.v[1]
				k.v[1] = k.v[3]
				k.v[3] = tmp

		elif vectorType == WhichVector.epsilon:
			k.mass = 0.0
			modk = math.sqrt(k.threeSquared())
			kt = math.sqrt(k.v[1]*k.v[1].conjugate() + k.v[2]*k.v[2].conjugate())
			kt_rt2 = math.sqrt(2.0) * kt
			modk_kt_rt2 = modk * kt_rt2

			ep1 = np.array([0.0, -k.v[2] / kt_rt2,            k.v[1] / kt_rt2,              0.0                ])
			ep2 = np.array([0.0, k.v[1]*k.v[3] / modk_kt_rt2, k.v[2]*k.v[3] / modk_kt_rt2, -kt**2 / modk_kt_rt2])

			if hel == Helicity.plus:
				k.v = -ep1 + 1j*ep2
			else:
				k.v = -ep1 - 1j*ep2
		return cls(k)

	@property
	def isu(self):
		if (self.vectorType == WhichVector.u) or (self.vectorType == WhichVector.ubar):
			return True
		else:
			return False

	@property
	def isv(self):
		if (self.vectorType == WhichVector.v) or (self.vectorType == WhichVector.vbar):
			return True
		else:
			return False

	@property
	def isbar(self):
		if (self.vectorType == WhichVector.ubar) or (self.vectorType == WhichVector.vbar):
			return True
		else:
			return False

	@property
	def iseps(self):
		if (self.vectorType == WhichVector.epsilon):
			return True
		else:
			return False

	#slash operator
	def __invert__(self):
		k = self.v
		matrix = [[0.0,          0.0,           k[0]-k[3],    -k[1]+1j*k[2]],
				  [0.0,          0.0,          -k[1]-1j*k[2],  k[0]+k[3]   ],
				  [k[0]+k[3],    k[1]-1j*k[2],  0.0,           0.0         ],
				  [k[1]+1j*k[2], k[0]-k[3],     0.0,           0.0         ]]
		return ComplexMatrix(matrix)

	def scalar(self, c):
		k = self.v * c
		return FourVector(k.tolist())

	def __add__(self, other):
		try:
			k = self.v + other.v
		except Exception:
			k = self.v + other
		return FourVector(k.tolist())

	def __radd__(self, other):
		k = self.v + other
		return FourVector(k.tolist())

	def __sub__(self, other):
		try:
			k = self.v - other.v
		except Exception:
			k = self.v - other
		return FourVector(k.tolist())

	def __rsub__(self, other):
		k = other - self.v
		return FourVector(k.tolist())

	def __mul__(self, other):
		try:
			k = self.v.dot(other)
		except Exception:
			try:
				k = self.v.dot(other.v)
			except Exception:
				try:
					k = self.v.dot(other.matrix)
				except Exception:
					k = self.v * other
		try:
			len(k)
			return FourVector(k.tolist())
		except Exception:
			return k

	def __truediv__(self, other):
		k = self.v / other
		return FourVector(k.tolist())

	def __rtruediv__(self, other):
		k = other / self.v
		return FourVector(k.tolist())

	def __getitem__(self, index):
		return self.v[index]

	def __str__(self):
		s_list = [format_num(elem) for elem in self.v.tolist()]
		s = ", ".join(s_list)
		return "(" + s + ")"

class ComplexMatrix:
	def __init__(self, matrix):
		self.matrix = np.array(matrix, dtype = np.complex_)

	def __mul__(self, other):
		try:
			self.v = self.v * other
			return self
		except Exception:
			v = self.matrix.dot(other.v)
			return FourVector(v.tolist())

	def __add__(self, other):
		for i in range(4):
			self.matrix[i][i] += other
		return ComplexMatrix(self.matrix.tolist())

	def __str__(self):
		sTot = []
		for row in self.matrix:
			s = ""
			for elem in row:
				s += format_num(elem, FORMAT_SPACE)
			sTot.append(s)
		return "\n".join(sTot)

# gamma matrices Weyl

gamma = []

gamma.append(ComplexMatrix([[  0,  0,  1,  0],
							[  0,  0,  0,  1],
							[  1,  0,  0,  0],
							[  0,  1,  0,  0]]))	#gamma 0

#controvariant (^mu)

gamma.append(ComplexMatrix([[  0,  0,  0,  1],
							[  0,  0,  1,  0],
							[  0, -1,  0,  0],
							[ -1,  0,  0,  0]]))	#gamma 1

gamma.append(ComplexMatrix([[  0,  0,  0,-1j],
							[  0,  0, 1j,  0],
							[  0, 1j,  0,  0],
							[-1j,  0,  0,  0]]))	#gamma 2

gamma.append(ComplexMatrix([[  0,  0,  1,  0],
							[  0,  0,  0, -1],
							[ -1,  0,  0,  0],
							[  0,  1,  0,  0]]))	#gamma 3

# covariant (_mu)

gamma.append(ComplexMatrix([[  0,  0, -1,  0],
							[  0,  0,  0,  1],
							[  1,  0,  0,  0],
							[  0, -1,  0,  0]]))	#gamma 3

gamma.append(ComplexMatrix([[  0,  0,  0, 1j],
							[  0,  0,-1j,  0],
							[  0,-1j,  0,  0],
							[ 1j,  0,  0,  0]]))	#gamma 2

gamma.append(ComplexMatrix([[  0,  0,  0, -1],
							[  0,  0, -1,  0],
							[  0,  1,  0,  0],
							[  1,  0,  0,  0]]))	#gamma 1

def VertexMul(s1, s2, controvariant = False):
	sign = 1
	if controvariant:
		sign = -1
	v = [0,0,0,0]
	for i in range(4):
		v[i] = (s1 * gamma[sign*i]) * s2
	return FourVector(v)
