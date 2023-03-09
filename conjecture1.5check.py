import mpmath
import numpy as np
import sympy as sp
import scipy as sci
from scipy.linalg import det
from numpy import *
from mpmath import *
from sympy.combinatorics.subsets import ksubsets
from sympy import sieve

#bernoulli number array
b = []
primes = list(sieve.primerange(2, 500))
text = open("conj15check.txt", 'w')

for i in range(0, 202):
	b.append(mpmath.bernfrac(i))

#helper function to compute powers mod p fast
def fast_power(a, n, p):
	if n == 0:
		return 1
	elif n == 1:
		return a % p
	else:
		temp = fast_power(a, n//2, p)
		if n % 2 == 0:
			return (temp * temp) % p
		else:
			return (temp * temp * a) % p

#helper function to compute inverses mod p fast
def inverse(a, p):
	if a % p == 0:
		raise ZeroDivisionError("denominator was 0 mod p")
	else:
		return fast_power(a, p - 2, p)

#helper function to compute bernoulli numbers mod p fast
def bern(n, p):
	if b[n][1] % p == 0:
		raise ZeroDivisionError("denominator was 0 mod p")
	else:
		return ((b[n][0] % p) * inverse(b[n][1] % p, p)) % p

#dimension of S_k
def cusp_dim(K):
	if K % 12 == 2:
		return K//12 - 1
	else:
		return K//12

# l_1 = 2, l_2 = 4, ...
def l_i(i, K):
	return 2 + 2 * (i - 1)

# k_i = K - l_i
def k_i(i, K):
	return K - l_i(i, K)

#helper function for the fourier coefficients mod p
def two_k_over_b(i, K, p):
	if bern(k_i(i, K), p) % p == 0:
		raise ZeroDivisionError("denominator was 0 mod p")
	else:
		return (2 * k_i(i, K) * inverse(bern(k_i(i, K), p), p)) % p

#helper function for the fourier coefficients mod p
def small_k_small_l(i, K, p):
	if bern(k_i(i, K), p) % p == 0 or l_i(i, K) % p == 0:
		raise ZeroDivisionError("denominator was 0 mod p")
	else:
		return ((k_i(i, K) * bern(l_i(i, K), p)) * inverse((l_i(i, K) * bern(l_i(i, K), p)), p)) % p 

#helper function for the fourier coefficients mod p
def large_K_small_l(i, K, p):
	if bern(K, p) % p == 0 or l_i(i, K) % p == 0:
		raise ZeroDivisionError("denominator was 0 mod p")
	else:
		return ((K * bern(l_i(i, K), p)) * inverse(l_i(i, K) * bern(K, p), p)) % p

#helper function for the fourier coefficients mod p
def sigma(n, k, p):
	sum = 0
	for i in range(1, n + 1):
		if n % i == 0:
			sum += fast_power(i, k, p)
			sum %= p
	return sum

#computes the fourier coefficients in equations 2.3 and 2.4 mod p
def four_coef(i, K, n, p):
	if i == 1:
		if (1 + small_k_small_l(i, K, p) - inverse(bern(k_i(i, K), p), p) - large_K_small_l(i, K, p)) % p == 0:
			raise ZeroDivisionError("denominator was 0 mod p")
		product = 0
		for m in range(1, n):
			product += (sigma(m, k_i(i, K) - 1, p) * sigma(n - m, l_i(i, K) - 1, p)) % p
			product %= p
		return ((small_k_small_l(i, K, p) * sigma(n, k_i(i, K) - 1, p) - (n * sigma(n, k_i(i, K) - 1, p) * inverse(bern(k_i(i, K), p), p))  + sigma(n, l_i(i, K) - 1, p) - two_k_over_b(i, K, p) * product
			- large_K_small_l(i, K, p) * sigma(n, K - 1, p)) * inverse(1 + small_k_small_l(i, K, p) - inverse(bern(k_i(i, K), p), p) - large_K_small_l(i, K, p), p)) % p
	else:
		if (1 + small_k_small_l(i, K, p) - large_K_small_l(i, K, p)) % p == 0:
			raise ZeroDivisionError("denominator was 0 mod p")
		product = 0
		for m in range(1, n):
			product += (sigma(m, k_i(i, K) - 1, p) * sigma(n - m, l_i(i, K) - 1, p)) % p
			product %= p
		return ((small_k_small_l(i, K, p) * sigma(n, k_i(i, K) - 1, p) + sigma(n, l_i(i, K) - 1, p) - two_k_over_b(i, K, p) * product
			- large_K_small_l(i, K, p) * sigma(n, K - 1, p)) * inverse(1 + small_k_small_l(i, K, p) - large_K_small_l(i, K, p), p)) % p

#designed for row >= col. returns the rank of the matrix modulo p
def row_reduce(matrix, p):
	row = len(matrix)
	col = len(matrix[0])
	if (col > row):
		return row_reduce([[matrix[r][c] for r in range(row)] for c in range(col)], p)
	for r in range(row):
		for c in range(col):
			matrix[r][c] %= p

	if col == 1:
		for r in range(row):
			if matrix[r][0] != 0:
				return 1
		return 0
	else:
		nonzero = -1
		for i in range(row):
			if matrix[i][0] != 0:
				nonzero = i
				break
		if nonzero != -1:
			if nonzero != 0:
				matrix[0] = [(matrix[0][j] - matrix[nonzero][j]) % p for j in range(col)]
			matrix[0] = [((matrix[0][i]) * inverse(matrix[0][0], p)) % p for i in range(col)]
			for r in range(1, row):
				matrix[r] = [ (matrix[r][c] - matrix[r][0] * matrix[0][c]) % p for c in range(col)]
			submatrix = [matrix[r][1:] for r in range(1, row)]
			return row_reduce(submatrix, p) + 1
		else:
			submatrix = [matrix[r][1:] for r in range(row)]
			return row_reduce(submatrix, p)

# checks for all choices of n k_i and l-i, the matrix corresponding to their normalized fourier coefficients has full rank.
for K in range(24, 202, 2):
	if K != 26:
		n = cusp_dim(K)
		length = 0

		for i in range(1, K):
			if k_i(i, K) < l_i(i, K):
				break
			length += 1

		for indices in list(ksubsets(range(1, length + 1), n)):
			flag = True
			for p in primes:
				if p > 1.2 * K:
					try:
						fc = []
						for i in indices:
							temp = []
							for j in range(1, n + 1):
								temp.append(four_coef(i, K, j, p))
							fc.append(temp)
						if row_reduce(np.array(fc), p) == n:
							flag = False
							break
					except ZeroDivisionError:
						pass
			if flag:
				text.write("K = " + str(K) + " counterexample, indices " + str(indices) + "\n")
				print("K = " + str(K) + " counterexample, indices " + str(indices) + "\n")
			else:
				text.write("K = " + str(K) + " done, indices " + str(indices) + "\n")
				print("K = " + str(K) + " done, indices " + str(indices))




