import mpmath
from mpmath import *

mp.dps = 2000
#bernoulli number array
b = []
#fourier coefficient array
fc = []
text = open("thm14check.txt", 'w')

for i in range(0, 400):
	b.append(mpmath.bernoulli(i))

# l_1 = 2, l_2 = 4, ...
def l_i(i, K):
	return 2 + 2 * (i - 1)

# k_i = K - l_i
def k_i(i, K):
	return K - l_i(i, K)

#dimension of S_k
def cusp_dim(K):
	if K % 12 == 2:
		return K//12 - 1
	else:
		return K//12

#helper function for the fourier coefficients
def two_k_over_b(i, K):
	return 2 * k_i(i, K) / b[k_i(i, K)]

#helper function for the fourier coefficients
def small_k_small_l(i, K):
	return (k_i(i, K) * b[l_i(i, K)])/(l_i(i, K) * b[k_i(i, K)])

#helper function for the fourier coefficients
def large_K_small_l(i, K):
	return (K * b[l_i(i, K)])/(l_i(i, K) * b[K])

#hardcoded to work for n prime
def sigma(n, k):
	return 1 + power(n, k)

#computes a_2(i), following equations 2.3 and 2.4
def four_coef_2(i, K):
	if i == 1:
		if 1 + small_k_small_l(i, K) - large_K_small_l(i, K) - 1/b[k_i(i, K)] == 0:
			return 0
		return (small_k_small_l(i, K) * sigma(2, k_i(i, K) - 1) - 2/b[k_i(i, K)] * sigma(2, k_i(i, K) - 1) + sigma(2, l_i(i, K) - 1) - two_k_over_b(i, K) 
			- large_K_small_l(i, K) * sigma(2, K - 1))/(1 + small_k_small_l(i, K) - large_K_small_l(i, K) - 1/b[k_i(i, K)])
	else:
		if 1 + small_k_small_l(i, K) - large_K_small_l(i, K)== 0:
			return 0
		return (small_k_small_l(i, K) * sigma(2, k_i(i, K) - 1) + sigma(2, l_i(i, K) - 1) - two_k_over_b(i, K) 
			- large_K_small_l(i, K) * sigma(2, K - 1))/(1 + small_k_small_l(i, K) - large_K_small_l(i, K))

#computes a_3(i), following equations 2.3 and 2.4
def four_coef_3(i, K):
	if i == 1:
		if 1 + small_k_small_l(i, K) - large_K_small_l(i, K) - 1/b[k_i(i, K)] == 0:
			return 0
		return (small_k_small_l(i, K) * sigma(2, k_i(i, K) - 1) - 3/b[k_i(i, K)] * sigma(3, k_i(i, K) - 1) + sigma(2, l_i(i, K) - 1) - two_k_over_b(i, K) 
			- large_K_small_l(i, K) * sigma(2, K - 1))/(1 + small_k_small_l(i, K) - large_K_small_l(i, K) - 1/b[k_i(i, K)])
	else:
		if 1 + small_k_small_l(i, K) - large_K_small_l(i, K)== 0:
			return 0
		return (small_k_small_l(i, K) * sigma(3, k_i(i, K) - 1) + sigma(3, l_i(i, K) - 1) - two_k_over_b(i, K) * (sigma(2, k_i(i, K) - 1) + sigma(2, l_i(i, K) - 1))
			- large_K_small_l(i, K) * sigma(3, K - 1))/(1 + small_k_small_l(i, K) - large_K_small_l(i, K))


# fc is an array such that fc[(K - 24)//2][i][j] is a_{j + 2}(i + 1) in the K case.
for K in range(24, 400, 2):
	temp = []
	for i in range(1, K//2):
		if l_i(i, K) > k_i(i, K):
			break
		else:
			temp.append((four_coef_2(i, K), four_coef_3(i, K)))
	fc.append(temp)


# computes the 3x3 determinant Delta_3 in section 3 for l_i, l_j, and l_k. Assumes i < j < k.
def det(i, j, k, K):
	return fc[K][i][1] * fc[K][j][0] - fc[K][i][1] * fc[K][k][0] - fc[K][j][1] * fc[K][i][0] + fc[K][j][1] * fc[K][k][0] + fc[K][k][1] * fc[K][i][0] -  fc[K][i][1] * fc[K][k][0] - fc[K][j][1] * fc[K][i][0] + fc[K][j][1] * fc[K][k][0] + fc[K][k][1] * fc[K][j][0]

# checks that the 3x3 determinant Delta_3 in section 3 is nonzero for K < 400; checks a 2x2 determinant if dim(S_k) = 2.
for K in range(24, 400, 2):
	# prints what weight K is being checked
	text.write(str(K) + "\n")
	print(str(K))
	if cusp_dim(K) == 2:
		for i in range(0, 1000):
			if i >= len(fc[(K - 24)//2]):
				break
			for j in range(i + 1, 1000):
				if j >= len(fc[(K - 24)//2]):
					break
				if abs(fc[(K - 24)//2][i][0] - fc[(K - 24)//2][j][0]) < power(10, -1):
					#outputs a counterexample (if one exists)
					text.write(str(i + 1) + ", " + str(j + 1) + ", " + str(K) + "\n")
	else:
		for i in range(0, 1000):
			if i >= len(fc[(K - 24)//2]):
				break
			for j in range(i + 1, 1000):
				if j >= len(fc[(K - 24)//2]):
					break
				for k in range(j + 1, 1000):
					if k >= len(fc[(K - 24)//2]):
						break
					if abs(det(i, j, k, (K - 24)//2)) < power(10, -3):
						#outputs a counterexample (if one exists)
						text.write(str(i + 1) + ", " + str(j + 1) + ", " + str(k + 1) + ", " + str(K) + "\n")




	
