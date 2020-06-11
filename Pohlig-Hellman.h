#ifndef POHLIG_HELLMAN_H
#define POHLIG_HELLMAN_H

#include <cassert>
#include <NTL\ZZ.h>
#include <NTL\matrix.h>
#include <vector>
#include <utility>
#include <NTL\tools.h>
#include "Bin-Method.h"

bool equiv_relation(const NTL::ZZ& a, const NTL::ZZ& b, const NTL::ZZ& m)
{
	assert(m >= 0);
	return (a - b) % m == 0;
}

std::vector<std::pair<NTL::ZZ, NTL::ZZ>> factorization(const NTL::ZZ& x)
{
	NTL::ZZ a = x;
	NTL::ZZ prime_div = (NTL::ZZ)2;
	NTL::ZZ deg = (NTL::ZZ)0;

	std::vector<std::pair<NTL::ZZ, NTL::ZZ>> res;
	while (a > 1)
	{
		while (a % prime_div == 0)
		{
			a /= prime_div;
			deg++;
		}
		if (deg != 0)
			res.emplace_back(std::make_pair(prime_div, deg));
		deg = 0;
		prime_div = NTL::NextPrime(prime_div + 1);
	}
	return res;
}

long max_in_vec(const std::vector<std::pair<NTL::ZZ, NTL::ZZ>>& v)
{
	long res = NTL::conv<long>(v[0].first);
	for (size_t i = 1; i < v.size(); i++)
		res = res < v[i].first ? NTL::conv<long>(v[i].first) : res;
	return res;
}

long max_deg_in_vec(const std::vector<std::pair<NTL::ZZ, NTL::ZZ>>& v)
{
	long res = NTL::conv<long>(v[0].second);
	for (size_t i = 1; i < v.size(); i++)
		res = res < v[i].second ? NTL::conv<long>(v[i].second) : res;
	return res;
}

NTL::Mat<NTL::ZZ> make_matrix(const std::vector<std::pair<NTL::ZZ, NTL::ZZ>>& v, const NTL::ZZ& a, const NTL::ZZ& p)
{
	NTL::Mat<NTL::ZZ> res;
	res.SetDims(v.size(), max_in_vec(v));
	for (long i = 0; i < v.size(); i++)
		for (long j = 0; j < v[i].first; j++)
			res[i][j] = bin_method(a, j * NTL::conv<long>(((p - 1) / v[i].first)));

	return res;
}

NTL::Mat<int> logs_of_b_to_a(const NTL::Mat<NTL::ZZ>& m, const std::vector<std::pair<NTL::ZZ, NTL::ZZ>>& v, const NTL::ZZ& a, const NTL::ZZ& b, const NTL::ZZ& p)
{
	NTL::Mat<int> res;
	res.SetDims(v.size(), 2);
	NTL::Mat<int> coefs;
	coefs.SetDims(v.size(), max_deg_in_vec(v));

	for (size_t i = 0; i < m.NumRows(); i++)
	{
		NTL::ZZ q_i = v[i].first;
		size_t x_0 = 0;
		while (!equiv_relation(m[i][x_0], bin_method(b, NTL::conv<long>((p - 1) / q_i)), p))
			x_0++;

		coefs[i][0] = x_0;
		NTL::ZZ sum = (NTL::ZZ)-1 * x_0;
		int k = 1;
		while (k < v[i].second)
		{
			for (size_t j = 0; j < v[i].first; j++)
				if (equiv_relation(m[i][j], bin_method(b * bin_method_mod(a, NTL::conv<long>(sum), p), NTL::conv<long>((p - 1) / bin_method(q_i, k + 1))), p))
				{
					coefs[i][k] = j;
					sum -= coefs[i][k] * NTL::power(q_i, k);
					k++;
					break;
				}
		}
	}

	for (size_t i = 0; i < coefs.NumRows(); i++)
	{
		res[i][0] = 0;
		res[i][1] = NTL::conv<int>(bin_method(v[i].first, NTL::conv<long>(v[i].second)));
		for (size_t j = 0; j < coefs.NumCols(); j++)
			if ((coefs[i][j] < v[i].first) && (coefs[i][j] > 0))
				res[i][0] += (coefs[i][j] * bin_method(v[i].first, NTL::conv<long>(j))) % res[i][1];
	}
	return res;
}

NTL::ZZ CRT(const NTL::Mat<int>& m)
{
	int M = 1;
	NTL::ZZ x = (NTL::ZZ)0;
	std::vector<int> M_i, M_i_1;
	int n = m.NumRows();
	M_i.resize(n);
	M_i_1.resize(n);
	for (size_t i = 0; i < n; i++)
		M *= m[i][1];

	for (size_t i = 0; i < n; i++)
	{
		M_i[i] = M / m[i][1];
		M_i_1[i] = NTL::conv<int>(NTL::PowerMod(M_i[i], -1, m[i][1]));
	}

	for (size_t i = 0; i < n; i++)
		x += (m[i][0] * M_i[i] * M_i_1[i]) % M;
	
	return x;
}

NTL::ZZ pohlig_hellman_alg(const NTL::ZZ& a, const NTL::ZZ& b, const NTL::ZZ& p)
{
	assert(NTL::GCD(b, p) == 1);
	std::vector<std::pair<NTL::ZZ, NTL::ZZ>> v = factorization(p - 1);
	NTL::Mat<NTL::ZZ> m = make_matrix(v, a, p);
	NTL::Mat<int> congruences = logs_of_b_to_a(m, v, a, b, p);
	NTL::ZZ res = CRT(congruences);
	return res < p ? res : res % p + 1;
}

#endif
