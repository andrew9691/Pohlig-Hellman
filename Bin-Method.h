#ifndef BIN_METHOD_H
#define BIN_METHOD_H

#include <cassert>
#include <NTL/ZZ.h>
#include <NTL/tools.h>
#include <math.h>

NTL::ZZ bin_method(const NTL::ZZ& x, const long& n)
{
	assert(n >= 0);
	if (n == 0)
		return (NTL::ZZ)1;
	NTL::ZZ N = NTL::conv<NTL::ZZ>(n);
	NTL::ZZ Y = (NTL::ZZ)1;
	NTL::ZZ Z = x;
	bool f;
	while (true)
	{
	m1:
		f = N % 2 == 0;
		N = std::floor(NTL::conv<double>(N / 2));
		if (f)
		{
		m2:
			Z *= Z;
			goto m1;
		}
		else
		{
			Y *= Z;
			if (N == 0)
				return Y;
			else
				goto m2;
		}
	}
}

NTL::ZZ bin_method_mod(const NTL::ZZ& x, const long& n, const NTL::ZZ& m)
{
	if (n < 0)
		return PowerMod(x, n, m);
	else
		return bin_method(x, n) % m;
}

#endif