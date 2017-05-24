//
//  Created by Seungcheon Baek (scienartbb@gmail.com)
//  Copyright CoCO SWING. All rights reserved.
//  ( http://cocoswing.com
//    http://github.com/scienartbb/ )
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// * Neither the name of tmatrix nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//         SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//


#if !defined(__SBMatrix_H__)
#define __SBMatrix_H__

#include <algorithm>
#include <sstream>
#include <complex>
#include <limits>
#include <float.h>
#include <string.h>
#include <assert.h>

#define _USE_MATH_DEFINES
#include <math.h>

#if defined(_WIN32)
#pragma warning(disable:4244)
#endif

#if !defined(min)
template<class T>
inline T min(const T& a, const T& b)
{
	return (a < b) ? a : b;
}
#endif

#if !defined(max)
template<class T>
inline T max(const T& a, const T& b)
{
	return (a > b) ? a : b;
}
#endif

template<class T>
inline T mod(const T& a, const T& b)
{
	return a % b;
}

template<>
inline float mod<float>(const float& a, const float& b)
{
	return fmod(a, b);
}

template<class T>
inline T abs(T a)
{
	return std::abs(a);
}

template<>
inline bool abs<bool>(bool a)
{
	return a;
}

template<>
inline unsigned char abs<unsigned char>(unsigned char a)
{
	return a;
}

template<>
inline unsigned short abs<unsigned short>(unsigned short a)
{
	return a;
}

template<>
inline unsigned int abs<unsigned int>(unsigned int a)
{
	return a;
}

template<class T>
inline T clamp(const T& a, const T& mn, const T& mx)
{
	if (a < mn)
		return mn;
	else if (a > mx)
		return mx;
	return a;
}

template<class T>
inline T fract(const T& a)
{
	return a - floor(a);
}

template<class T>
inline T fract(const T& a, const T& mn, const T& mx)
{
	const T d = mx - mn;
	if (a < 0)
	{
		return mod((a - mx), d) + mx;
	}
	return mod((a - mn), d) + mn;
}

template<class T>
inline T calcRateRepeat(const T& value, const T& by)
{
	return T(mod(value, by)) / by;
}

template<class T>
inline T calcRateReflect(const T& value, const T& by)
{
	const T f = T(mod(value, 2 * by)) / by;
	if (f >= 1)
		return 2 - f;
	return f;
}

template<class T>
inline T fractoffset(const T& a, const T& d, const T& o)
{
	// [a b c d e f g h i j k]
	// - (7,3,1) => h
	// - fractoffset(7,3,1) => g
	return T(a / d) * d + o;
}

template<class T>
inline T round(const T& a)
{
	const T fa = floor(a);
	if (a - fa < 0.5)
	{
		return fa;
	}
	return ceil(a);
}

template <class T>
inline T conj_(const T& a)
{
	return a;
}

template <class T>
inline std::complex<T> conj_(const std::complex<T>& a)
{
	return std::conj(a);
}

template <class T>
inline bool isZero(const T& a)
{
#if 1
	return abs(a) < std::numeric_limits<T>::epsilon();
#else
	return abs(a) < 0.0078125f;	// 2 ^ -7
#endif
}

template <class T>
inline bool isZero(const std::complex<T>& a)
{
	return isZero(a.real()) && isZero(a.imag());
}

template <class T>
inline bool isEqual(const T& a, const T& b)
{
	return isZero(a-b);
}

template <class T>
inline bool isIdentity(const T& a)
{
	return std::abs(a-1) < std::numeric_limits<T>::epsilon();
}

template <class T>
inline bool isIdentity(const std::complex<T>& a)
{
	return isIdentity(a.real()) && isZero(a.imag());
}

// to decrease floating point error
template <class T>
inline bool absless(const T& a, const T& b)
{
	return std::abs(a) < std::abs(b);
}

// to decrease floating point error
template <class T>
inline bool absless(const std::complex<T>& a, const std::complex<T>& b)
{
	return std::norm(a) < std::norm(b);
}

// a x + b = 0 �� �� x
template <class T>
inline bool roots(const T& a, const T& b, T* r)
{
	// a x + b = 0
	// x = - b / a
	if (isZero(a)) {
		return false;
	}
	*r = -b / a;
	return true;
};

//a x^2 + b x + c = 0 �� �� x
template <class T>
inline bool roots(const T& a, const T& b, const T& c, T* r1, T* r2)
{
	// a x^2 + b x + c = 0
	// -b + sqrt(b^2 - 4 a c) / 2 a
	if (isZero(a)) {
		roots(b, c, r1);
		*r2 = *r1;
		return true;
	}
	else {
		if (isZero(c)) {
			*r2 = 0;
			if (!roots(a,b, r1)) {
				return false;
			}
			return true;
		}
		else {
			const T d = b*b-4*a*c;
			if (d < 0) {
				return false;
			}
			const T sd = sqrt(d);
			const T ad = 0.5*a;
			*r1 = (-b + sd) * ad;
			*r2 = (-b - sd) * ad;
			return true;
		}
	}
	return false;
}

//a x^3 + b x^2 + c x + d = 0 �� �� x
template <class T>
inline bool roots(const T& a, const T& b, const T& c, const T& d, T* r1, T* r2, T* r3)
{
	// a x^2 + b x + c = 0
	if (isZero(a)) {
		if (!roots(b, c, d, r1, r2)) {
			return false;
		}
		*r3 = *r2;
		return true;
	}
	
	if (isZero(d)) {
		if (!roots(a, b, c, r2, r3)) {
			return false;
		}
		*r1 = 0;
		return true;
	}

	if (isZero(b) && isZero(c)) {
		*r1 = *r2 = *r3 = pow(-d/a, 1.0f/3);
		return true;
	}

	float s=0.5;
	s=(-(a*s*s*s+b*s*s+c*s+d)/(3*a*s*s+2*b*s+c))+s;
	s=(-(a*s*s*s+b*s*s+c*s+d)/(3*a*s*s+2*b*s+c))+s;
	s=(-(a*s*s*s+b*s*s+c*s+d)/(3*a*s*s+2*b*s+c))+s;
	s=(-(a*s*s*s+b*s*s+c*s+d)/(3*a*s*s+2*b*s+c))+s;
	s=(-(a*s*s*s+b*s*s+c*s+d)/(3*a*s*s+2*b*s+c))+s;
	s=(-(a*s*s*s+b*s*s+c*s+d)/(3*a*s*s+2*b*s+c))+s;
	s=(-(a*s*s*s+b*s*s+c*s+d)/(3*a*s*s+2*b*s+c))+s;
	s=(-(a*s*s*s+b*s*s+c*s+d)/(3*a*s*s+2*b*s+c))+s;
	s=(-(a*s*s*s+b*s*s+c*s+d)/(3*a*s*s+2*b*s+c))+s;
	s=(-(a*s*s*s+b*s*s+c*s+d)/(3*a*s*s+2*b*s+c))+s;
	*r1=s;
	const float b2 = b+s*a, c2 = -d/s;
	if (!roots(a, b2, c2, r2, r3)) {
		return false;
	}

	//// x^3 + e x^2 + f x + g = 0
	//// x = z-e/3
	//// z^3 + pz + q

	//const T e = b/a, f = c/a, g = d/a;
	//const T p=f-e*e/3, q=2*e*e*e/27-e*f/3+g;
	//const T t0=-q/2+sqrt(q*q+4*p*p*p/27)/2, u0 = q/2+sqrt(q*q+4*p*p*p/27)/2;
	////const T t1=-q/2-sqrt(q*q+4*p*p*p/27)/2, u1 = q/2-sqrt(q*q+4*p*p*p/27)/2;
	//T z1 = t0 < 0 ? pow(-t0, 1.0f/3) : pow(t0, 1.0f/3);
	//T z2 = u0 < 0 ? pow(-u0, 1.0f/3) : pow(u0, 1.0f/3);

	//const T x = (z1-z2)-e/3;
	//const T /*a2 = a, */b2 = b+x*a, c2 = -d/x;
	//if (!roots(a, b2, c2, r2, r3)) {
	//	return false;
	//}

	//*r1 = x;
	return true;
}

//a x^3 + b x^2 + c x + d = 0
template <class T>
inline float roots_cubic(const T& a, const T& b, const T& c, const T& d)
{
	float s=0.5;
	s=(-(a*s*s*s+b*s*s+c*s+d)/(3*a*s*s+2*b*s+c))+s;
	s=(-(a*s*s*s+b*s*s+c*s+d)/(3*a*s*s+2*b*s+c))+s;
	s=(-(a*s*s*s+b*s*s+c*s+d)/(3*a*s*s+2*b*s+c))+s;
	s=(-(a*s*s*s+b*s*s+c*s+d)/(3*a*s*s+2*b*s+c))+s;
	s=(-(a*s*s*s+b*s*s+c*s+d)/(3*a*s*s+2*b*s+c))+s;
	s=(-(a*s*s*s+b*s*s+c*s+d)/(3*a*s*s+2*b*s+c))+s;
	s=(-(a*s*s*s+b*s*s+c*s+d)/(3*a*s*s+2*b*s+c))+s;
	s=(-(a*s*s*s+b*s*s+c*s+d)/(3*a*s*s+2*b*s+c))+s;
	s=(-(a*s*s*s+b*s*s+c*s+d)/(3*a*s*s+2*b*s+c))+s;
	s=(-(a*s*s*s+b*s*s+c*s+d)/(3*a*s*s+2*b*s+c))+s;
	return s;
}

template <class T, unsigned int R, unsigned int C>
class SBMatrix
{
public:
	union {
		T m[R][C];
		//struct {
		//	T x, y, z, w;
		//};
		//struct {
		//	T r, g, b, a;
		//};
	};

	inline unsigned int rsize(void) const
	{
		return R;
	}

	inline unsigned int csize(void) const
	{
		return C;
	}

	friend std::ostream& operator <<(std::ostream& os, const SBMatrix<T, R, C>& v)
	{
		os << '[';
		unsigned int c, r;
		for (r = 0; r < R; ++r)
		{
			os.width(2);
			for (c = 0; c < C-1; ++c)
			{
				os << v.m[r][c] << " ";
			}
			os << v.m[r][c] << ';';
		}
		os << ']';
		return os;
	}

	friend std::istream& operator >>(std::istream& is, SBMatrix<T, R, C>& v)
	{
		std::stringbuf sb;
		is.get(sb);
		std::string s;
		s = sb.str();
		unsigned int r = 0, c = 0;
		const char* p = s.c_str();
		T t;
		while (*p)
		{
			switch (*p)
			{
			case ' ':
				++p;
				break;
			case ';':
				for (; c < C; ++c)
				{
					v.m[r][c] = 0;
				}
				++r;
				c = 0;
				++p;
				break;
			case '[':
				++p;
				break;
			case ']':
				if (r < R)
				{
					for (; c < C; ++c)
					{
						v.m[r][c] = 0;
					}
					++r;
				}
				++p;
				break;
			default:
				{
					std::istringstream iss(p, std::ios_base::in);
					iss >> t;
					if (c < C)
					{
						v.m[r][c++] = t;
					}
					p = strpbrk(p, " [];");
				}
				break;
			}
		}
		for (; r < R; ++r)
		{
			for (c = 0; c < C; ++c)
			{
				v.m[r][c] = 0;
			}
		}
		return is;
	}

private:

	SBMatrix<T, R, C>& operator=(const T& a)
	{
		for (unsigned int c = 0; c < C; ++c)
		{
			for (unsigned int r = 0; r < R; ++r)
			{
				m[r][c] = a;
			}
		}
		return *this;
	}

public:
	inline bool operator==(const SBMatrix<T, R, C>& a) const
	{
		for (unsigned int r = 0; r < R; ++r)
		{
			for (unsigned int c = 0; c < C; ++c)
			{
				if (!isEqual(a.m[r][c], m[r][c]))
				{
					return false;
				}
			}
		}
		return true;
	}

	inline bool operator!=(const SBMatrix<T, R, C>& a) const
	{
		for (unsigned int r = 0; r < R; ++r)
		{
			for (unsigned int c = 0; c < C; ++c)
			{
				if (!isEqual(a.m[r][c], m[r][c]))
				{
					return true;
				}
			}
		}
		return false;
	}

	inline bool operator<(const SBMatrix<T, R, C>& a) const
	{
		for (unsigned int r = 0; r < R; ++r)
		{
			for (unsigned int c = 0; c < C; ++c)
			{
				if (!isEqual(a.m[r][c], m[r][c]))
				{
					if (a.m[r][c] < m[r][c])
					{
						return true;
					}
					else {
						return false;
					}
				}
			}
		}
		return false;
	}

	inline bool operator>(const SBMatrix<T, R, C>& a) const
	{
		for (unsigned int r = 0; r < R; ++r)
		{
			for (unsigned int c = 0; c < C; ++c)
			{
				if (!isEqual(a.m[r][c], m[r][c]))
				{
					if (a.m[r][c] > m[r][c])
					{
						return true;
					}
					else {
						return false;
					}
				}
			}
		}
		return false;
	}

	inline bool operator<=(const SBMatrix<T, R, C>& a) const
	{
		return (*this < a) || (*this == a);
	}

	inline bool operator>=(const SBMatrix<T, R, C>& a) const
	{
		return (*this > a) || (*this == a);
	}

	inline SBMatrix<T, R, C> operator-(void) const
	{
		SBMatrix<T, R, C> res;
		for (unsigned int c = 0; c < C; ++c)
		{
			for (unsigned int r = 0; r < R; ++r)
			{
				res.m[r][c] = -m[r][c];
			}
		}
		return res;
	}
	inline SBMatrix<T, R, C>& operator=(const SBMatrix<T, R, C>& a)
	{
		for (unsigned int c = 0; c < C; ++c)
		{
			for (unsigned int r = 0; r < R; ++r)
			{
				m[r][c] = a.m[r][c];
			}
		}
		return *this;
	}

	inline SBMatrix<T, R, C> operator+(const SBMatrix<T, R, C>& a) const
	{
		SBMatrix<T, R, C> res;
		for (unsigned int c = 0; c < C; ++c)
		{
			for (unsigned int r = 0; r < R; ++r)
			{
				res.m[r][c] = m[r][c] + a.m[r][c];
			}
		}
		return res;
	}
	inline SBMatrix<T, R, C> operator+=(const SBMatrix<T, R, C>& a)
	{
		for (unsigned int c = 0; c < C; ++c)
		{
			for (unsigned int r = 0; r < R; ++r)
			{
				m[r][c] += a.m[r][c];
			}
		}
		return *this;
	}
	inline SBMatrix<T, R, C> operator+(const T& a) const
	{
		SBMatrix<T, R, C> res;
		for (unsigned int c = 0; c < C; ++c)
		{
			for (unsigned int r = 0; r < R; ++r)
			{
				res.m[r][c] = m[r][c] + a;
			}
		}
		return res;
	}
	inline SBMatrix<T, R, C> operator-(const SBMatrix<T, R, C>& a) const
	{
		SBMatrix<T, R, C> res;
		for (unsigned int c = 0; c < C; ++c)
		{
			for (unsigned int r = 0; r < R; ++r)
			{
				res.m[r][c] = m[r][c] - a.m[r][c];
			}
		}
		return res;
	}
	inline SBMatrix<T, R, C>& operator-=(const SBMatrix<T, R, C>& a)
	{
		for (unsigned int c = 0; c < C; ++c)
		{
			for (unsigned int r = 0; r < R; ++r)
			{
				m[r][c] -= a.m[r][c];
			}
		}
		return *this;
	}
	inline SBMatrix<T, R, C> operator-(const T& a) const
	{
		SBMatrix<T, R, C> res;
		for (unsigned int c = 0; c < C; ++c)
		{
			for (unsigned int r = 0; r < R; ++r)
			{
				res.m[r][c] = m[r][c] - a;
			}
		}
		return res;
	}
	inline SBMatrix<T, R, C> operator*(const T& a) const
	{
		SBMatrix<T, R, C> res;
		for (unsigned int c = 0; c < C; ++c)
		{
			for (unsigned int r = 0; r < R; ++r)
			{
				res.m[r][c] = m[r][c] * a;
			}
		}
		return res;
	}
	inline SBMatrix<T, R, C> operator*=(const T& a)
	{
		for (unsigned int c = 0; c < C; ++c)
		{
			for (unsigned int r = 0; r < R; ++r)
			{
				m[r][c] *= a;
			}
		}
		return *this;
	}

	friend SBMatrix<T, R, C> operator*(const T& a, const SBMatrix<T, R, C>& b)
	{
		return b * a;
	}

	template<unsigned int R1, unsigned int C1>
	inline SBMatrix<T, R1, C1> sub(const unsigned int r, const unsigned int c) const
	{
		SBMatrix<T, R1, C1> res;
		unsigned int r1, c1;
		for (c1 = 0; c1 < C1; ++c1)
		{
			for (r1 = 0; r1 < R1; ++r1)
			{
				res.m[r1][c1] = m[r+r1][c+c1];
			}
		}
		return res;
	}
	inline SBMatrix<T, 1, C> row(const unsigned int r) const
	{
		SBMatrix<T, 1, C> res;
		for (unsigned int c = 0; c < C; ++c)
		{
			res.m[0][c] = m[r][c];
		}
		return res;
	}

	inline SBMatrix<T, R, 1> col(const unsigned int c) const
	{
		SBMatrix<T, R, 1> res;
		for (unsigned int r = 0; r < R; ++r)
		{
			res.m[r][0] = m[r][c];
		}
		return res;
	}

	inline void row(const unsigned int r, const SBMatrix<T, 1, C>& a)
	{
		for (unsigned int c = 0; c < C; ++c)
		{
			m[r][c] = a.m[0][c];
		}
	}

	inline void col(const unsigned int c, const SBMatrix<T, R, 1>& a)
	{
		for (unsigned int r = 0; r < R; ++r)
		{
			m[r][c] = a.m[r][0];
		}
	}

	inline void swapRow(const unsigned int r1, const unsigned int r2)
	{
		SBMatrix<T, 1, C> t = row(r1);
		row(r1, row(r2));
		row(r2, t);
	}

	inline void swapCol(const unsigned int c1, const unsigned int c2)
	{
		const SBMatrix<T, R, 1> t = col(c1);
		col(c1, col(c2));
		col(c2, t);
	}

	template<unsigned int CA>
	inline SBMatrix<T, R, CA> operator*(const SBMatrix<T, C, CA>& a) const
	{
		SBMatrix<T, R, CA> res;
		if (R == 1 && CA == 1)
		{
			// [x1 y1 z1]*[x2;y2;z2]
			T f = 0;
			for (unsigned int c = 0; c < C; ++c)
				f += m[0][c] * a.m[c][0];
			res.m[0][0] = f;
		}
		else {
			for (unsigned int c = 0; c < CA; ++c)
			{
				for (unsigned int r = 0; r < R; ++r)
				{
					SBMatrix<T, 1, 1> v = row(r) * a.col(c);
					res.m[r][c] = v.m[0][0];
				}
			}
		}
		return res;
	}

	inline SBMatrix<T, R, C> operator/(const T& a) const
	{
		const T b = 1/a;
		SBMatrix<T, R, C> res;
		for (unsigned int c = 0; c < C; ++c)
		{
			for (unsigned int r = 0; r < R; ++r)
			{
				res.m[r][c] = m[r][c] * b;
			}
		}
		return res;
	}

	inline SBMatrix<T, R, C> operator/=(const T& a)
	{
		const T b = 1/a;
		for (unsigned int c = 0; c < C; ++c)
		{
			for (unsigned int r = 0; r < R; ++r)
			{
				m[r][c] *= b;
			}
		}
		return *this;
	}

	inline SBMatrix<T, R, C> diag(void) const
	{
		SBMatrix<T, R, C> res;
		for (unsigned int c = 0; c < C; ++c)
		{
			for (unsigned int r = 0; r < R; ++r)
			{
				if (r == c)
					res[r][c] = m[r][c];
				else
					res[r][c] = 0;
			}
		}
		return true;
	}
	inline SBMatrix<T, R, C> triu(void) const
	{
		SBMatrix<T, R, C> res;
		for (unsigned int c = 0; c < C; ++c)
		{
			for (unsigned int r = 0; r < R; ++r)
			{
				if (r <= c)
					res[r][c] = m[r][c];
				else
					res[r][c] = 0;
			}
		}
		return res;
	}
	inline SBMatrix<T, R, C> tril(void) const
	{
		SBMatrix<T, R, C> res;
		for (unsigned int c = 0; c < C; ++c)
		{
			for (unsigned int r = 0; r < R; ++r)
			{
				if (r >= c)
					res[r][c] = m[r][c];
				else
					res[r][c] = 0;
			}
		}
		return res;
	}
	inline SBMatrix<T, R, C> conj(void) const
	{
		SBMatrix<T, R, C> res;
		for (unsigned int c = 0; c < C; ++c)
		{
			for (unsigned int r = 0; r < R; ++r)
			{
				res.m[r][c] = conj_(m[r][c]);
			}
		}
		return res;
	}
	inline SBMatrix<T, C, R> trans(void) const
	{
		SBMatrix<T, C, R> res;
		for (unsigned int c = 0; c < C; ++c)
		{
			for (unsigned int r = 0; r < R; ++r)
			{
				res.m[c][r] = m[r][c];
			}
		}
		return res;
	}
	inline SBMatrix<T, C, R> hermitian(void) const
	{
		return conj().trans();
	}

	template<unsigned int R1, unsigned int C1>
	inline SBMatrix<T, R1, C1> except(const unsigned int re, const unsigned int ce) const
	{
		SBMatrix<T, R1, C1> res;
		unsigned int r2 = 0, c2 = 0;

		for (unsigned int c = 0; c < C1; ++c)
		{
			if (ce == c)
				c2 = 1;
			r2 = 0;
			for (unsigned int r = 0; r < R1; ++r)
			{
				if (re == r)
					r2 = 1;
				res.m[r][c] = m[r+r2][c+c2];
			}
		}
		return res;
	}

	inline T adj(const unsigned int r, const unsigned int c) const
	{
		if (R == 1 || C == 1)
			return 1;
		else if (R == 2 || C == 2)
		{
			return m[1-r][1-c];
		}
		SBMatrix<T, R-1, C-1> v;
		if (R-1 < 1 || C-1 < 1)
			return 0;
		const unsigned int rr = R - 1, cc = C - 1;
		unsigned int r2 = 0, c2 = 0;

		for (unsigned int c1 = 0; c1 < cc; ++c1)
		{
			if (c1 == c)
				c2 = 1;
			r2 = 0;
			for (unsigned int r1 = 0; r1 < rr; ++r1)
			{
				if (r1 == r)
					r2 = 1;
				v.m[r1][c1] = m[r1+r2][c1+c2];
			}
		}
		return v.det();
	}

	inline T cofactor(const unsigned int r, const unsigned int c) const
	{
		if ((r+c)&1)
			return -adj(r,c);
		return adj(r,c);
	}

	inline T det(void) const
	{
		if (R == 1 && C == 1)
			return m[0][0];
		T res = 0;
		for (unsigned int c = 0; c < C; ++c)
		{
			if (!isZero(m[0][c]))
				res += m[0][c] * cofactor(0, c);
		}
		return res;
	}

	inline SBMatrix<T, R, C> inv(void) const
	{
		SBMatrix<T, R, C> res;
		const T iD = det();
		if (isZero(iD))
		{
			for (unsigned int c = 0; c < C; ++c)
			{
				for (unsigned int r = 0; r < R; ++r)
				{
					res.m[r][c] = 0;
				}
			}
			return res;
		}
		for (unsigned int c = 0; c < C; ++c)
		{
			for (unsigned int r = 0; r < R; ++r)
			{
				res.m[r][c] = cofactor(r,c) / iD;
			}
		}
		return res.hermitian();
	}

	inline SBMatrix<T, R, R> projmat(void) const
	{
		// M inv(M' M) M'
		const SBMatrix<T, C, R> h = hermitian();
		const SBMatrix<T, C, C> d = h * *this;
		return *this * d.inv() * h;
	}

	inline SBMatrix<T, R, R> reflect(void) const
	{
		SBMatrix<T, R, R> I;
		I.identity();
		return projmat() * 2.0 - I;
	}

	inline T trace(void) const
	{
		T res = 0;
		const unsigned int dsize = R<C?R:C;
		for (unsigned int d = 0; d < dsize; ++d)
		{
			res += m[d][d];
		}
		return res;
	}

	inline unsigned int rank(void) const
	{
		if (C < R)
		{
			return trans().rank();
		}
		SBMatrix<T, R, C> temp;
		temp = *this;
		const unsigned int dsize = R<C?R:C;

		unsigned int res = 0;
		bool b = false;
		for (unsigned int d = 0; d < dsize; ++d)
		{
			b = isZero(temp.m[d][d]);
			if (b)
			{
				for (unsigned int r = d + 1; r < R; ++r)
				{
					if (!isZero(temp.m[r][d]))
					{
						temp.swapRow(r, d);
						b = false;
						break;
					}
				}
			}

			if (!b)
			{
				for (unsigned int r = d + 1; r < R; ++r)
				{
					if (!isZero(temp.m[r][d]))
					{
						// decreasing error
						if (absless(temp.m[r][d], temp.m[d][d]))
						{
							temp.row(r, temp.row(r) * (-temp.m[d][d] / temp.m[r][d]) + temp.row(d));
						}
						else {
							temp.row(r, temp.row(d) * (-temp.m[r][d] / temp.m[d][d]) + temp.row(r));
						}
					}
				}
				++ res;
			}
		}
		return res;
	}

	inline unsigned int nullity(void) const
	{
		// rank-nullity theorem
		// rank(A) + nullity(A) = n
		return C - rank();
	}

	inline unsigned int dim(void) const
	{
		if (C > 1)
			return C;
		return R;
	}

	inline bool pivot(const unsigned int d, unsigned int* pr, unsigned int* pc) const
	{
		if (isZero(m[d][d]))
		{
			for (unsigned int c = d+1; c < C; ++c)
			{
				if (!isZero(m[d][c]))
				{
					*pr = d;
					*pc = c;
					return true;
				}
			}
			for (unsigned int r = d + 1; r < R; ++r)
			{
				for (unsigned int c = d + 1; c < C; ++c)
				{
					if (!isZero(m[r][c]))
					{
						*pr = r;
						*pc = c;
						return true;
					}
				}
			}
			return false;
		}
		*pr = *pc = d;
		return true;
	}

	inline unsigned int findNonZeroRow(const unsigned int column) const
	{
		for (unsigned int r = column+1; r < R; ++r)
		{
			if (!isZero(m[r][column]))
			{
				return r;
			}
		}
		return -1;
	}

	inline void eliminateRowD(const unsigned int r, const unsigned int c)
	{
		for (unsigned int r1 = r + 1; r1 < R; ++r1)
		{
			if (!isZero(m[r1][c]))
			{
				// decreasing error
				if (absless(m[r1][c], m[r][c]))
				{
					const T coef = -m[r][c] / m[r1][c];
					row(r1, row(r) * coef + row(r1));
				}
				else {
					const T coef = -m[r1][c] / m[r][c];
					row(r1, row(r) * coef + row(r1));
				}
			}
		}
	}

	inline void eliminateRowU(const unsigned int r, const unsigned int c)
	{
		for (unsigned int r1 = 0; r1 < r; ++r1)
		{
			if (!isZero(m[r1][c]))
			{
				// decreasing error
				if (absless(m[r1][c], m[r][c]))
				{
					const T coef = -m[r][c] / m[r1][c];
					row(r1, row(r) * coef + row(r1));
				}
				else {
					const T coef = -m[r1][c] / m[r][c];
					row(r1, row(r) * coef + row(r1));
				}
			}
		}
	}

	// gaussian elimination
	inline SBMatrix<T, R, C> gausselim(void) const
	{
		SBMatrix<T, R, C> res;
		res = *this;

		unsigned int r = 0, c = 0;
		bool done = false;
		while (!done)
		{
			if (isZero(res.m[r][c]))
			{
				unsigned int pr;
				if (res.findNonZeroRow(c, &pr))
				{
					res.swapRow(r, pr);
				}
				else {
					++c;
					if (c >= C)
						done = true;
				}
			}
			else {
				res.eliminateRowD(r, c);
				if (!isIdentity(res.m[r][c]))
				{
					res.row(r, res.row(r) / res.m[r][c]);
				}
				++r;	++c;
				if (r >= R || c >= C)
					done = true;
			}
		}

		return res;
	}

	// gauss-jordan elimination
	inline SBMatrix<T, R, C> gjordelim(void) const
	{
		SBMatrix<T, R, C> res;
		res = *this;

		unsigned int r = 0, c = 0;
		bool done = false;
		while (!done)
		{
			if (isZero(res.m[r][c]))
			{
				unsigned int pr;
				if (res.findNonZeroRow(c, &pr))
				{
					res.swapRow(r, pr);
				}
				else {
					++c;
					if (c >= C)
						done = true;
				}
			}
			else {
				res.eliminateRowD(r, c);
				if (!isIdentity(res.m[r][c]))
				{
					res.row(r, res.row(r) / res.m[r][c]);
				}
				res.eliminateRowU(r, c);
				++r;	++c;
				if (r >= R || c >= C)
					done = true;
			}
		}

		return res;
	}

	inline SBMatrix<T, C, 1> rowSpace(unsigned int index) const
	{
		const SBMatrix<T, C, R> t = trans();
		const SBMatrix<T, C, R> b = t.gausselim();
		SBMatrix<T, C, 1> res;
		for (unsigned int r = 0; r < C; ++r)
		{
			res.m[r][0] = b.m[r][index];
		}
		return res;
	}

	inline SBMatrix<T, R, 1> colSpace(unsigned int index) const
	{
		const SBMatrix<T, R, C> b = gausselim();
		SBMatrix<T, R, 1> res;
		for (unsigned int r = 0; r < R; ++r)
		{
			res.m[r][0] = b.m[r][index];
		}
		return res;
	}

	inline bool isIdentity(void) const
	{
		for (unsigned int c = 0; c < C; ++c)
		{
			for (unsigned int r = 0; r < R; ++r)
			{
				if (r == c)
				{
					if (!isEqual(m[r][c], 1)) {
						return false;
					}
				}
				else
				{
					if (!isEqual(m[r][c], 0)) {
						return false;
					}
				}
			}
		}
		return true;
	}

	inline void identity(void)
	{
		for (unsigned int c = 0; c < C; ++c)
		{
			for (unsigned int r = 0; r < R; ++r)
			{
				if (r == c)
				{
					m[r][c] = 1;
				}
				else
				{
					m[r][c] = 0;
				}
			}
		}
	}

	inline bool isZeroes(void) const
	{
		for (unsigned int c = 0; c < C; ++c)
		{
			for (unsigned int r = 0; r < R; ++r)
			{
				if (::abs(m[r][c]) > std::numeric_limits<T>::epsilon())
					return false;
			}
		}
		return true;
	}

	inline void zeros(void)
	{
		for (unsigned int c = 0; c < C; ++c)
		{
			for (unsigned int r = 0; r < R; ++r)
			{
				m[r][c] = 0;
			}
		}
	}

	inline void zerosL(void)
	{
		for (unsigned int c = 0; c < C; ++c)
		{
			for (unsigned int r = c+1; r < R; ++r)
			{
				m[r][c] = 0;
			}
		}
	}

	inline void zerosU(void)
	{
		for (unsigned int c = 0; c < C; ++c)
		{
			unsigned int d = c < R ? c : R;
			for (unsigned int r = 0; r < d; ++r)
			{
				m[r][c] = 0;
			}
		}
	}

	inline void ones(void)
	{
		for (unsigned int c = 0; c < C; ++c)
		{
			for (unsigned int r = 0; r < R; ++r)
			{
				m[r][c] = 1;
			}
		}
	}

	inline void rand(void)
	{
		for (unsigned int c = 0; c < C; ++c)
		{
			for (unsigned int r = 0; r < R; ++r)
			{
				m[r][c] = T(::rand() % 1000) / T(1000);
			}
		}
	}

	inline void clamp(const T& mn, const T& mx)
	{
		unsigned int r, c;
		for (r = 0; r < R; ++r)
		{
			for (c = 0; c < C; ++c)
			{
				m[r][c] = clamp(m[r][c], mn, mx);
			}
		}
	}

	inline void pow(const T& b)
	{
		for (unsigned int c = 0; c < C; ++c)
		{
			for (unsigned int r = 0; r < R; ++r)
			{
				m[r][c] = ::pow(m[r][c], b);
			}
		}
	}

	inline void abs(void)
	{
		for (unsigned int c = 0; c < C; ++c)
		{
			for (unsigned int r = 0; r < R; ++r)
			{
				m[r][c] = ::abs(m[r][c]);
			}
		}
	}
};

template<class T>
class SBMatrix<T, 0, 0>
{
public:
	T m[1][1];

	T det() const
	{
		return 0;
	}
};

template<class T>
class SBMatrix<T, 1, 1>
{
public:
	union {
		T m[1][1];
		struct {
			T x;
		};
	};

	inline void zeros(void)
	{
		memset(m, 0, sizeof(m));
	}

	inline T det(void)
	{
		return m[0][0];
	}

	inline SBMatrix<T, 1, 1> operator-(void) const
	{
		SBMatrix<T, 1, 1> res = {-m[0][0]};
		return res;
	}

	inline SBMatrix<T, 1, 1> operator+(const SBMatrix<T, 1, 1>& a) const
	{
		SBMatrix<T, 1, 1> res;
		res.m[0][0] = m[0][0] + a.m[0][0];
		return res;
	}

	inline SBMatrix<T, 1, 1> operator-(const SBMatrix<T, 1, 1>& a) const
	{
		SBMatrix<T, 1, 1> res;
		res.m[0][0] = m[0][0] - a.m[0][0];
		return res;
	}

	inline SBMatrix<T, 1, 1> operator*(const T& a) const
	{
		SBMatrix<T, 1, 1> res;
		res.m[0][0] = m[0][0] * a;
		return res;
	}

	inline bool operator==(const SBMatrix<T, 1, 1>& a) const
	{
		if (!isEqual(a.m[0][0], m[0][0]))
		{
			return false;
		}
		return true;
	}

	inline bool operator!=(const SBMatrix<T, 1, 1>& a) const
	{
		if (!isEqual(a.m[0][0], m[0][0]))
		{
			return true;
		}
		return false;
	}

	inline SBMatrix<T, 1, 1> col(const unsigned int c) const
	{
		SBMatrix<T, 1, 1> res;
		res.m[0][0] = m[0][c];
		return res;
	}
};

template<class T>
class SBMatrix<T, 2, 1>
{
public:
	union {
		T m[2][1];
		struct {
			T x, y;
		};
	};

	inline void zeros(void)
	{
		memset(m, 0, sizeof(m));
	}

	inline SBMatrix<T, 2, 1> operator-(void) const
	{
		SBMatrix<T, 2, 1> res = {-m[0][0], -m[1][0]};
		return res;
	}

	inline SBMatrix<T, 2, 1> operator+(const SBMatrix<T, 2, 1>& a) const
	{
		SBMatrix<T, 2, 1> res;
		res.m[0][0] = m[0][0] + a.m[0][0];
		res.m[1][0] = m[1][0] + a.m[1][0];
		return res;
	}

	inline SBMatrix<T, 2, 1> operator-(const SBMatrix<T, 2, 1>& a) const
	{
		SBMatrix<T, 2, 1> res;
		res.m[0][0] = m[0][0] - a.m[0][0];
		res.m[1][0] = m[1][0] - a.m[1][0];
		return res;
	}

	inline SBMatrix<T, 2, 1> operator*(const T& a) const
	{
		SBMatrix<T, 2, 1> res;
		res.m[0][0] = m[0][0] * a;
		res.m[1][0] = m[1][0] * a;
		return res;
	}

	inline SBMatrix<T, 2, 1> operator/(const T& a) const
	{
#if 0
		SBMatrix<T, 2, 1> res;
		res.m[0][0] = m[0][0] / a;
		res.m[1][0] = m[1][0] / a;
		return res;
#else
		const T b = 1/a;
		return *this*b;
#endif
	}

	inline SBMatrix<T, 2, 1> operator*(const SBMatrix<T, 2, 1>& a) const
	{
		SBMatrix<T, 2, 1> res;
		res.m[0][0] = m[0][0] * a;
		res.m[1][0] = m[1][0] * a;
		return res;
	}

	inline SBMatrix<T, 2, 1> operator/=(const T& a)
	{
		const T b = 1/a;
		m[0][0] *= b;
		m[1][0] *= b;
		return *this;
	}

	inline SBMatrix<T, 2, 1> operator*=(const T& a)
	{
		m[0][0] *= a;
		m[1][0] *= a;
		return *this;
	}

	inline bool operator==(const SBMatrix<T, 2, 1>& a) const
	{
		if (!isEqual(a.m[0][0], m[0][0]))
		{
			return false;
		}
		if (!isEqual(a.m[1][0], m[1][0]))
		{
			return false;
		}
		return true;
	}

	inline bool operator!=(const SBMatrix<T, 2, 1>& a) const
	{
		if (!isEqual(a.m[0][0], m[0][0]))
		{
			return true;
		}
		if (!isEqual(a.m[1][0], m[1][0]))
		{
			return true;
		}
		return false;
	}

	inline SBMatrix<T, 2, 1> col(const unsigned int c) const
	{
		SBMatrix<T, 2, 1> res;
		res.m[0][0] = m[0][c];
		res.m[1][0] = m[1][c];
		return res;
	}
};

template<class T>
class SBMatrix<T, 3, 1>
{
public:
	union {
		T m[3][1];
		struct {
			T x, y, z;
		};
	};

	inline void zeros(void)
	{
		memset(m, 0, sizeof(m));
	}

	inline SBMatrix<T, 3, 1> operator-(void) const
	{
		const SBMatrix<T, 3, 1> res = {-m[0][0], -m[1][0], -m[2][0]};
		return res;
	}

	inline SBMatrix<T, 3, 1> operator+(const SBMatrix<T, 3, 1>& a) const
	{
		SBMatrix<T, 3, 1> res;
		res.m[0][0] = m[0][0] + a.m[0][0];
		res.m[1][0] = m[1][0] + a.m[1][0];
		res.m[2][0] = m[2][0] + a.m[2][0];
		return res;
	}

	inline SBMatrix<T, 3, 1> operator-(const SBMatrix<T, 3, 1>& a) const
	{
		SBMatrix<T, 3, 1> res;
		res.m[0][0] = m[0][0] - a.m[0][0];
		res.m[1][0] = m[1][0] - a.m[1][0];
		res.m[2][0] = m[2][0] - a.m[2][0];
		return res;
	}

	inline SBMatrix<T, 3, 1> operator*(const T& a) const
	{
		SBMatrix<T, 3, 1> res;
		res.m[0][0] = m[0][0] * a;
		res.m[1][0] = m[1][0] * a;
		res.m[2][0] = m[2][0] * a;
		return res;
	}

	inline SBMatrix<T, 3, 1> operator/(const T& a) const
	{
#if 0
		SBMatrix<T, 3, 1> res;
		res.m[0][0] = m[0][0] / a;
		res.m[1][0] = m[1][0] / a;
		res.m[2][0] = m[2][0] / a;
		return res;
#else
		const T b = 1/a;
		return *this*b;
#endif
	}

	inline SBMatrix<T, 3, 1> operator*=(const T& a)
	{
		m[0][0] *= a;
		m[1][0] *= a;
		m[2][0] *= a;
		return *this;
	}

	inline SBMatrix<T, 3, 1> operator/=(const T& a)
	{
		const T b = 1/a;
		m[0][0] *= b;
		m[1][0] *= b;
		m[2][0] *= b;
		return *this;
	}

	inline bool operator==(const SBMatrix<T, 3, 1>& a) const
	{
		if (!isEqual(a.m[0][0], m[0][0]))
		{
			return false;
		}
		if (!isEqual(a.m[1][0], m[1][0]))
		{
			return false;
		}
		if (!isEqual(a.m[2][0], m[2][0]))
		{
			return false;
		}
		return true;
	}

	inline bool operator!=(const SBMatrix<T, 3, 1>& a) const
	{
		if (!isEqual(a.m[0][0], m[0][0]))
		{
			return true;
		}
		if (!isEqual(a.m[1][0], m[1][0]))
		{
			return true;
		}
		if (!isEqual(a.m[2][0], m[2][0]))
		{
			return true;
		}
		return false;
	}

	inline SBMatrix<T, 3, 1> col(const unsigned int c) const
	{
		SBMatrix<T, 3, 1> res;
		res.m[0][0] = m[0][c];
		res.m[1][0] = m[1][c];
		res.m[2][0] = m[2][c];
		return res;
	}
};

template<class T>
class SBMatrix<T, 4, 1>
{
public:
	union {
		T m[4][1];
		struct {
			T x, y, z, w;
		};
	};

	inline void zeros(void)
	{
		memset(m, 0, sizeof(m));
	}

	inline SBMatrix<T, 4, 1> operator-(void) const
	{
		const SBMatrix<T, 4, 1> res = {-m[0][0], -m[1][0], -m[2][0], -m[3][0]};
		return res;
	}

	inline SBMatrix<T, 4, 1> operator+(const SBMatrix<T, 4, 1>& a) const
	{
		SBMatrix<T, 4, 1> res;
		res.m[0][0] = m[0][0] + a.m[0][0];
		res.m[1][0] = m[1][0] + a.m[1][0];
		res.m[2][0] = m[2][0] + a.m[2][0];
		res.m[3][0] = m[3][0] + a.m[3][0];
		return res;
	}

	inline SBMatrix<T, 4, 1> operator-(const SBMatrix<T, 4, 1>& a) const
	{
		SBMatrix<T, 4, 1> res;
		res.m[0][0] = m[0][0] - a.m[0][0];
		res.m[1][0] = m[1][0] - a.m[1][0];
		res.m[2][0] = m[2][0] - a.m[2][0];
		res.m[3][0] = m[3][0] - a.m[3][0];
		return res;
	}

	inline SBMatrix<T, 4, 1> operator*(const T& a) const
	{
		SBMatrix<T, 4, 1> res;
		res.m[0][0] = m[0][0] * a;
		res.m[1][0] = m[1][0] * a;
		res.m[2][0] = m[2][0] * a;
		res.m[3][0] = m[3][0] * a;
		return res;
	}

	inline SBMatrix<T, 4, 1> operator/(const T& a) const
	{
#if 0
		SBMatrix<T, 4, 1> res;
		res.m[0][0] = m[0][0] / a;
		res.m[1][0] = m[1][0] / a;
		res.m[2][0] = m[2][0] / a;
		res.m[3][0] = m[3][0] / a;
		return res;
#else
		const T b = 1/a;
		return *this*b;
#endif
	}

	inline SBMatrix<T, 4, 1> operator*=(const T& a)
	{
		m[0][0] *= a;
		m[1][0] *= a;
		m[2][0] *= a;
		m[3][0] *= a;
		return *this;
	}

	inline SBMatrix<T, 4, 1> operator/=(const T& a)
	{
		const T b = 1/a;
		m[0][0] *= b;
		m[1][0] *= b;
		m[2][0] *= b;
		m[3][0] *= b;
		return *this;
	}

	inline bool operator==(const SBMatrix<T, 4, 1>& a) const
	{
		if (!isEqual(a.m[0][0], m[0][0]))
		{
			return false;
		}
		if (!isEqual(a.m[1][0], m[1][0]))
		{
			return false;
		}
		if (!isEqual(a.m[2][0], m[2][0]))
		{
			return false;
		}
		if (!isEqual(a.m[3][0], m[3][0]))
		{
			return false;
		}
		return true;
	}

	inline bool operator!=(const SBMatrix<T, 4, 1>& a) const
	{
		if (!isEqual(a.m[0][0], m[0][0]))
		{
			return true;
		}
		if (!isEqual(a.m[1][0], m[1][0]))
		{
			return true;
		}
		if (!isEqual(a.m[2][0], m[2][0]))
		{
			return true;
		}
		if (!isEqual(a.m[3][0], m[3][0]))
		{
			return true;
		}
		return false;
	}

	inline SBMatrix<T, 4, 1> col(const unsigned int c) const
	{
		SBMatrix<T, 4, 1> res;
		res.m[0][0] = m[0][c];
		res.m[1][0] = m[1][c];
		res.m[2][0] = m[2][c];
		res.m[3][0] = m[3][c];
		return res;
	}
};


template <class T, unsigned int R, unsigned int C>
inline bool isZero(const SBMatrix<T, R, C>& a)
{
	for (unsigned int r = 0; r < R; ++r) {
		for (unsigned int c = 0; c < C; ++c) {
			if (!isZero(a.m[r][c])) {
				return false;
			}
		}
	}
	return true;
}


template<class T>
SBMatrix<T, 2, 1> mul(const SBMatrix<T, 3, 3>& a, const SBMatrix<T, 2, 1>& b)
{
	SBMatrix<T, 2, 1> r;
	r.m[0][0] = a.m[0][0] * b.m[0][0] + a.m[0][1] * b.m[1][0] + a.m[0][2];
	r.m[1][0] = a.m[1][0] * b.m[0][0] + a.m[1][1] * b.m[1][0] + a.m[1][2];
	const T w = a.m[2][0] * b.m[0][0] + a.m[2][1] * b.m[1][0] + a.m[2][2];
	return r * (1/w);
}

template<class T>
SBMatrix<T, 3, 1> mul(const SBMatrix<T, 4, 4>& a, const SBMatrix<T, 3, 1>& b)
{
	SBMatrix<T, 3, 1> r;
	r.m[0][0] = a.m[0][0] * b.m[0][0] + a.m[0][1] * b.m[1][0] + a.m[0][2] * b.m[2][0] + a.m[0][3];
	r.m[1][0] = a.m[1][0] * b.m[0][0] + a.m[1][1] * b.m[1][0] + a.m[1][2] * b.m[2][0] + a.m[1][3];
	r.m[2][0] = a.m[2][0] * b.m[0][0] + a.m[2][1] * b.m[1][0] + a.m[2][2] * b.m[2][0] + a.m[2][3];
	const T w = a.m[3][0] * b.m[0][0] + a.m[3][1] * b.m[1][0] + a.m[3][2] * b.m[2][0] + a.m[3][3];
	return r / w;
}

template<class T>
SBMatrix<T, 4, 1> mul(const SBMatrix<T, 4, 4>& a, const SBMatrix<T, 4, 1>& b)
{
	return a * b;
}

template<class T, unsigned int R>
inline T length2(const SBMatrix<T, R, 1>& a)
{
	T value = 0;
	for (unsigned int r = 0; r < R; ++r)
	{
		value += a.m[r][0] * a.m[r][0];
	}
	return value;
}

template<class T, unsigned int R>
inline T length(const SBMatrix<T, R, 1>& a)
{
	return sqrt(length2(a));
}

template<class T, unsigned int R>
inline SBMatrix<T, R, 1> normalize(const SBMatrix<T, R, 1>& a)
{
	const T t = length(a);
	assert(!isZero(t));
	const T l = 1 / t;
	return a * l;
}

template<class T, unsigned int R>
inline T dot(const SBMatrix<T, R, 1>& a, const SBMatrix<T, R, 1>& b)
{
	T value = 0;
	for (unsigned int r = 0; r < R; ++r)
	{
		value += a.m[r][0] * b.m[r][0];
	}
	return value;
}

template<class T>
inline T area(const SBMatrix<T, 2, 1>& a, const SBMatrix<T, 2, 1>& b)
{
	return a.m[0][0]*b.m[1][0]-a.m[1][0]*b.m[0][0];
}

template<class T>
inline SBMatrix<T, 3, 1> cross(const SBMatrix<T, 3, 1>& a, const SBMatrix<T, 3, 1>& b)
{
	SBMatrix<T, 3, 1> res;
	res.m[0][0] = a.m[1][0] * b.m[2][0] - a.m[2][0] * b.m[1][0];
	res.m[1][0] = a.m[2][0] * b.m[0][0] - a.m[0][0] * b.m[2][0];
	res.m[2][0] = a.m[0][0] * b.m[1][0] - a.m[1][0] * b.m[0][0];
	return res;
}

template<class T, unsigned int R>
inline SBMatrix<T, R, 1> proj(const SBMatrix<T, R, 1>& a, const SBMatrix<T, R, 1> &b)
{
	// proj a
	//     b
	return (dot(a, b) / dot(b, b)) * b;
}

template<class T, unsigned int C>
inline void permutation(const SBMatrix<unsigned int, 1, C>& zbased, SBMatrix<T, C, C>* p)
{
	p->zeros();
	p->isIdentity();
	for (unsigned int c = 0; c < C; ++c)
	{
		p->m[c][zbased.m[0][c]] = 1;
	}
}

template<class T, unsigned int R, unsigned int C, unsigned int C2>
SBMatrix<T, R, C+C2> augment(SBMatrix<T, R, C>& a, SBMatrix<T, R, C2>& b)
{
	SBMatrix<T, R, C+C2> res;
	const unsigned int c1 = C+C2-1;
	unsigned int r, c;
	for (r = 0; r < R; ++r)
	{
		for (c = 0; c < c1; ++c)
		{
			res.m[r][c] = a.m[r][c];
		}
		for (c = 0; c < C2; ++c)
		{
			res.m[r][C+c] = b.m[r][c];
		}
	}
	return res;
}

template<class T, unsigned int R>
inline T dist2(const SBMatrix<T, R, 1>& a, const SBMatrix<T, R, 1>& b)
{
	T s = 0;
	for (unsigned int r = 0; r < R; ++r)
	{
		T t = a.m[r][0]-b.m[r][0];
		s += t*t;
	}
	return s;
}

template<class T, unsigned int R>
inline T dist2(const SBMatrix<T, R, 1>& a, const SBMatrix<T, R, 1>& b, const SBMatrix<T, R, 1>& c)
{
	return dist2(a, b) + dist2(b,c);
}

template<class T, unsigned int R>
inline T dist2(const SBMatrix<T, R, 1>& a, const SBMatrix<T, R, 1>& b, const SBMatrix<T, R, 1>& c, const SBMatrix<T, R, 1>& d)
{
	return dist2(a, b) + dist2(b,c) + dist2(c,d);
}

template<class T, unsigned int R>
inline T dist(const SBMatrix<T, R, 1>& a, const SBMatrix<T, R, 1>& b)
{
	return sqrt(dist2(a, b));
}

template<class T, unsigned int R>
inline T dist(const SBMatrix<T, R, 1>& a, const SBMatrix<T, R, 1>& b, const SBMatrix<T, R, 1>& c)
{
	return dist(a, b) + dist(b,c);
}

template<class T, unsigned int R>
inline T dist(const SBMatrix<T, R, 1>& a, const SBMatrix<T, R, 1>& b, const SBMatrix<T, R, 1>& c, const SBMatrix<T, R, 1>& d)
{
	return dist(a, b) + dist(b,c) + dist(c,d);
}

template<class T, unsigned int R>
inline bool lu(const SBMatrix<T, R, R>& a, SBMatrix<T, R, R>* l, SBMatrix<T, R, R>* u)
{
	// lu, m == n
	// l lower triangular, u upper triangular
	*u = a;
	l->zerosU();
	for (unsigned int d = 0; d < R; ++d)
	{
		unsigned int pr, pc;
		if (u->pivot(d, &pr, &pc))
		{
			if (pr != d)
				return false;
			for (unsigned int r = d + 1; r < R; ++r)
			{
				l->m[r][d] = u->m[r][d] / u->m[d][d];
				for (unsigned int c = d; c < R; ++c)
				{
					u->m[r][c] += -l->m[r][d] * u->m[d][c];
				}
			}
		}
		else {
			return false;
		}
	}
	return true;
}

template<class T, unsigned int R>
inline bool qr(const SBMatrix<T, R, R>& a, SBMatrix<T, R, R>* q, SBMatrix<T, R, R>* r)
{
	// a=q*r
	struct V {
		SBMatrix<T, R, 1> a, u, e;
	};
	V v[R];
	r->zerosL();
	for (unsigned int x = 0; x < R; ++x)
	{
		v[x].u = v[x].a = a.col(x);
		for (int y = x - 1; y >= 0; --y)
		{
			r->m[y][x] = dot(v[x].a, v[y].e);
			v[x].u -= r->m[y][x] * v[y].e;
		}
		// sqrt(V'*V)
		r->m[x][x] = length(v[x].u);
		v[x].e = v[x].u / r->m[x][x];

		q->col(x, v[x].e);
	}
	return true;
}

template<class T>
inline SBMatrix<T, 2, 2> rotate2(const T& rad)
{
	SBMatrix<T, 2, 2> res;
	const T c = cos(rad), s = sin(rad);
	res.m[0][0] = c;	res.m[0][1] = -s;
	res.m[1][0] = s;	res.m[1][1] = c;
	return res;
}

template<class T>
inline SBMatrix<T, 3, 3> rotate3(const T& rad)
{
	SBMatrix<T, 3, 3> res;
	const T c = cos(rad), s = sin(rad);
	res.m[0][0] = c;	res.m[0][1] = -s;	res.m[0][2] = 0;
	res.m[1][0] = s;	res.m[1][1] = c;	res.m[1][2] = 0;
	res.m[2][0] = 0;	res.m[2][1] = 0;	res.m[2][2] = 1;
	return res;
}

template<class T>
inline SBMatrix<T, 4, 4> rotateX(const T& rad)
{
	SBMatrix<T, 4, 4> res;
	const T c = cos(rad), s = sin(rad);
	res.m[0][0] = 1;	res.m[0][1] = 0;	res.m[0][2] = 0;	res.m[0][3] = 0;
	res.m[1][0] = 0;	res.m[1][1] = c;	res.m[1][2] = -s;	res.m[1][3] = 0;
	res.m[2][0] = 0;	res.m[2][1] = s;	res.m[2][2] = c;	res.m[2][3] = 0;
	res.m[3][0] = 0;	res.m[3][1] = 0;	res.m[3][2] = 0;	res.m[3][3] = 1;
	return res;
}

template<class T>
inline SBMatrix<T, 4, 4> rotateY(const T& rad)
{
	SBMatrix<T, 4, 4> res;
	const T c = cos(rad), s = sin(rad);
	res.m[0][0] = c;	res.m[0][1] = 0;	res.m[0][2] = s;	res.m[0][3] = 0;
	res.m[1][0] = 0;	res.m[1][1] = 1;	res.m[1][2] = 0;	res.m[1][3] = 0;
	res.m[2][0] = -s;	res.m[2][1] = 0;	res.m[2][2] = c;	res.m[2][3] = 0;
	res.m[3][0] = 0;	res.m[3][1] = 0;	res.m[3][2] = 0;	res.m[3][3] = 1;
	return res;
}

template<class T>
inline SBMatrix<T, 4, 4> rotateZ(const T& rad)
{
	SBMatrix<T, 4, 4> res;
	const T c = cos(rad), s = sin(rad);
	res.m[0][0] = c;	res.m[0][1] = -s;	res.m[0][2] = 0;	res.m[0][3] = 0;
	res.m[1][0] = s;	res.m[1][1] = c;	res.m[1][2] = 0;	res.m[1][3] = 0;
	res.m[2][0] = 0;	res.m[2][1] = 0;	res.m[2][2] = 1;	res.m[2][3] = 0;
	res.m[3][0] = 0;	res.m[3][1] = 0;	res.m[3][2] = 0;	res.m[3][3] = 1;
	return res;
}

template<class T>
inline SBMatrix<T, 4, 4> rotateXDegree(const T& degree)
{
	return rotateX(degree / 180.0f * M_PI);
}

template<class T>
inline SBMatrix<T, 4, 4> rotateYDegree(const T& degree)
{
	return rotateY(degree / 180.0f * M_PI);
}

template<class T>
inline SBMatrix<T, 4, 4> rotateZDegree(const T& degree)
{
	return rotateZ(degree / 180.0f * M_PI);
}

template<class T>
inline SBMatrix<T, 4, 4> rotate(const T& x, const T& y, const T& z, const T& rad)
{
	SBMatrix<T, 4, 4> res;
	const T s=sin(rad), c=cos(rad), c1 = 1-c;

	res.m[0][0] = c + x * x * c1;
	res.m[0][1] = x * y * c1 - z * s;
	res.m[0][2] = x * z * c1 + y * s;
	res.m[0][3] = 0;
	res.m[1][0] = x * y * c1 + z * s;
	res.m[1][1] = c + y * y * c1;
	res.m[1][2] = y * z * c1 - x * s;
	res.m[1][3] = 0;
	res.m[2][0] = x * z * c1 - y * s;
	res.m[2][1] = y * z * c1 + x * s;
	res.m[2][2] = c + z * z * c1;
	res.m[2][3] = 0;
	res.m[3][0] = 0;
	res.m[3][1] = 0;
	res.m[3][2] = 0;
	res.m[3][3] = 1;

	return res;
}

template<class T>
inline SBMatrix<T, 4, 4> rotateDegree(const T& x, const T& y, const T& z, const T& degree)
{
	return rotate(x, y, z, degree / 180.0f * M_PI);
}

template<class T>
inline SBMatrix<T, 3, 3> translate(const T& x, const T& y)
{
	SBMatrix<T, 3, 3> res;
	res.m[0][0] = 1;	res.m[0][1] = 0;	res.m[0][2] = x;
	res.m[1][0] = 0;	res.m[1][1] = 1;	res.m[1][2] = y;
	res.m[2][0] = 0;	res.m[2][1] = 0;	res.m[2][2] = 1;
	return res;
}

template<class T>
inline SBMatrix<T, 4, 4> translate(const T& x, const T& y, const T& z)
{
	SBMatrix<T, 4, 4> res;
	res.m[0][0] = 1;	res.m[0][1] = 0;	res.m[0][2] = 0;	res.m[0][3] = x;
	res.m[1][0] = 0;	res.m[1][1] = 1;	res.m[1][2] = 0;	res.m[1][3] = y;
	res.m[2][0] = 0;	res.m[2][1] = 0;	res.m[2][2] = 1;	res.m[2][3] = z;
	res.m[3][0] = 0;	res.m[3][1] = 0;	res.m[3][2] = 0;	res.m[3][3] = 1;
	return res;
}

template<class T>
inline SBMatrix<T, 2, 2> scale2(const T& x, const T& y)
{
	SBMatrix<T, 2, 2> res;
	res.m[0][0] = x;	res.m[0][1] = 0;
	res.m[1][0] = 0;	res.m[1][1] = y;
	return res;
}

template<class T>
inline SBMatrix<T, 3, 3> scale(const T& x, const T& y)
{
	SBMatrix<T, 3, 3> res;
	res.m[0][0] = x;	res.m[0][1] = 0;	res.m[0][2] = 0;
	res.m[1][0] = 0;	res.m[1][1] = y;	res.m[1][2] = 0;
	res.m[2][0] = 0;	res.m[2][1] = 0;	res.m[2][2] = 1;
	return res;
}

template<class T>
inline SBMatrix<T, 4, 4> scale(const T& x, const T& y, const T& z)
{
	SBMatrix<T, 4, 4> res;
	res.m[0][0] = x;	res.m[0][1] = 0;	res.m[0][2] = 0;	res.m[0][3] = 0;
	res.m[1][0] = 0;	res.m[1][1] = y;	res.m[1][2] = 0;	res.m[1][3] = 0;
	res.m[2][0] = 0;	res.m[2][1] = 0;	res.m[2][2] = z;	res.m[2][3] = 0;
	res.m[3][0] = 0;	res.m[3][1] = 0;	res.m[3][2] = 0;	res.m[3][3] = 1;
	return res;
}

template<class T>
inline SBMatrix<T, 4, 4> skew(const SBMatrix<T, 3, 1>& a1, const SBMatrix<T, 3, 1>& a2, const T& rad)
{
	const SBMatrix<T, 3, 1> v1 = a1 / length(a1);
	const SBMatrix<T, 3, 1> v2 = a2 / length(a2);
	const T d = dot(v1, v2), axisangle = acos(d);
	if (rad >= axisangle || rad <= (axisangle - M_PI))
	{
		SBMatrix<T, 4, 4> res;
		res.identity();
		return res;
	}
	SBMatrix<T, 3, 1> c;
	c = cross(v1, v2);
	c /= length(c);

	SBMatrix<T, 3, 1> v1o = cross(v2, c);
	SBMatrix<T, 4, 4> r;
	r.m[0][0] = c.m[0][0];	r.m[0][1] = v1o.m[0][0];	r.m[0][2] = v2.m[0][0];	r.m[0][3] = 0;
	r.m[1][0] = c.m[1][0];	r.m[1][1] = v1o.m[1][0];	r.m[1][2] = v2.m[1][0];	r.m[1][3] = 0;
	r.m[2][0] = c.m[2][0];	r.m[2][1] = v1o.m[2][0];	r.m[2][2] = v2.m[2][0];	r.m[2][3] = 0;
	r.m[3][0] = 0;	r.m[3][1] = 0;	r.m[3][2] = 0;	r.m[3][3] = 1;

	T perp = sqrt(1-d*d);
	T s = tan(rad+acos(perp)) * perp - d;

	SBMatrix<T, 4, 4> k;
	k.m[0][0] = 1;	k.m[0][1] = 0;	k.m[0][2] = 0;	k.m[0][3] = 0;
	k.m[1][0] = 0;	k.m[1][1] = 1;	k.m[1][2] = s;	k.m[1][3] = 0;
	k.m[2][0] = 0;	k.m[2][1] = 0;	k.m[2][2] = 1;	k.m[2][3] = 0;
	k.m[3][0] = 0;	k.m[3][1] = 0;	k.m[3][2] = 0;	k.m[3][3] = 1;
	return r.trans() * k * r;
}

template<class T>
inline SBMatrix<T, 4, 4> skewDegree(const SBMatrix<T, 3, 1>& a1, const SBMatrix<T, 3, 1>& a2, const T& degree)
{
	return skew(a1, a2, degree / 180.0f * M_PI);
}

template<class T>
inline SBMatrix<T, 4, 4> lookat(const SBMatrix<T, 3, 1>& eye, const SBMatrix<T, 3, 1>& at, const SBMatrix<T, 3, 1>& up)
{
	SBMatrix<T, 4, 4> res;

	SBMatrix<T, 3, 1> x, y, z;
	z = normalize(eye - at);
	x = normalize(cross(up, z));
	y = normalize(cross(z, x));

	res.m[0][0] = x.m[0][0];	res.m[0][1] = y.m[0][0];	res.m[0][2] = z.m[0][0];	res.m[0][3] = -dot(x, eye);
	res.m[1][0] = x.m[1][0];	res.m[1][1] = y.m[1][0];	res.m[1][2] = z.m[1][0];	res.m[1][3] = -dot(y, eye);
	res.m[2][0] = x.m[2][0];	res.m[2][1] = y.m[2][0];	res.m[2][2] = z.m[2][0];	res.m[2][3] = -dot(z, eye);
	res.m[3][0] = 0;	res.m[3][1] = 0;	res.m[3][2] = 0;	res.m[3][3] = 1;
	return res;
}

template<class T>
inline SBMatrix<T, 4, 4> nproj(const T& w, const T& h, const T& n, const T& f)
{
	SBMatrix<T, 4, 4> res;
	// x: [-w/2,w/2] -> [-1,1], y: [-h/2,h/2] -> [-1,1], z -> 0
	res.m[0][0] = 2/w;	res.m[0][1] = 0;	res.m[0][2] = 0;	res.m[0][3] = 0;
	res.m[1][0] = 0;	res.m[1][1] = 2/h;	res.m[1][2] = 0;	res.m[1][3] = 0;
	res.m[2][0] = 0;	res.m[2][1] = 0;	res.m[2][2] = 0;	res.m[2][3] = 0;
	res.m[3][0] = 0;	res.m[3][1] = 0;	res.m[3][2] = 0;	res.m[3][3] = 1;
	return res;
}

template<class T>
inline SBMatrix<T, 4, 4> oproj(const T& w, const T& h, const T& n, const T& f)
{
	return scale(2 / w, 2 / h, -1 / (f - n));
}

template<class T>
inline SBMatrix<T, 4, 4> pproj(const T& fovY, const T& aspect_xpery, const T& n, const T& f)
{
	SBMatrix<T, 4, 4> res;
	const T h = 2 * tan(fovY / 2) * n, w = aspect_xpery * h;

	res.m[0][0] = 2*n/w;	res.m[0][1] = 0;	res.m[0][2] = 0;	res.m[0][3] = 0;
	res.m[1][0] = 0;	res.m[1][1] = 2*n/h;	res.m[1][2] = 0;	res.m[1][3] = 0;
	res.m[2][0] = 0;	res.m[2][1] = 0;	res.m[2][2] = -(f+n)/(f-n);	res.m[2][3] = -(2*f*n)/(f-n);
	res.m[3][0] = 0;	res.m[3][1] = 0;	res.m[3][2] = -1;	res.m[3][3] = 0;
	return res;
}

template<class T>
inline SBMatrix<T, 4, 4> pprojDegree(const T& fovYdegree, const T& aspect_xpery, const T& n, const T& f)
{
	return pproj(fovYdegree / 180.0f * M_PI, aspect_xpery, n, f);
}

template<class T>
inline SBMatrix<T, 4, 4> euler(const T& head, const T& pitch, const T& roll)
{
	SBMatrix<T, 4, 4> res;
	const T hc = cos(head), hs = sin(head);
	const T pc = cos(pitch), ps = sin(pitch);
	const T rc = cos(roll), rs = sin(roll);
	res.m[0][0] = rc * hc - rs * ps * hs;
	res.m[0][1] = -rs * pc;
	res.m[0][2] = rc * hs + rs * ps * hc;
	res.m[0][3] = 0;
	res.m[1][0] = rs * hc + rc * ps * hs;
	res.m[1][1] = rc * pc;
	res.m[1][2] = rs * hs - rc * ps * hc;
	res.m[1][3] = 0;
	res.m[2][0] = -pc * hs;
	res.m[2][1] = ps;
	res.m[2][2] = pc * hc;
	res.m[2][3] = 0;
	res.m[3][0] = 0;
	res.m[3][1] = 0;
	res.m[3][2] = 0;
	res.m[3][3] = 1;
	return res;
}

template<class T>
inline SBMatrix<T, 4, 4> affineTransform2D(const T& sx, const T& sy, const T& rz, const T& tx, const T& ty)
{
	//[ sx*cos(rz), -sy*sin(rz), 0, tx]
	//[ sx*sin(rz),  sy*cos(rz), 0, ty]
	//[          0,           0, 1,  0]
	//[          0,           0, 0,  1]
	const float zc = cos(rz), zs = sin(rz);
	SBMatrix<T, 4, 4> res = {
		sx*zc, -sy*zs, 0, tx
		, sx*zs,  sy*zc, 0, ty
		, 0, 0, 1, 0
		, 0, 0, 0, 1};
	return res;
}

template<class T>
inline SBMatrix<T, 4, 4> affineTransform(const T& sx, const T& sy, const T& sz, const T& rx, const T& ry, const T& rz, const T& tx, const T& ty, const T& tz)
{
	//[ sx*cos(ry)*cos(rz), -sy*(cos(rx)*sin(rz) - cos(rz)*sin(rx)*sin(ry)),  sz*(sin(rx)*sin(rz) + cos(rx)*cos(rz)*sin(ry)), tx]
	//[ sx*cos(ry)*sin(rz),  sy*(cos(rx)*cos(rz) + sin(rx)*sin(ry)*sin(rz)), -sz*(cos(rz)*sin(rx) - cos(rx)*sin(ry)*sin(rz)), ty]
	//[        -sx*sin(ry),                              sy*cos(ry)*sin(rx),                              sz*cos(rx)*cos(ry), tz]
	//[                  0,                                               0,                                               0,  1]
	const float xc = cos(rx), yc = cos(ry), zc = cos(rz), 
		xs = sin(rx), ys = sin(ry), zs = sin(rz);
	SBMatrix<T, 4, 4> res = {
		sx*yc*zc, -sy*(xc*zs - zc*xs*ys),  sz*(xs*zs + xc*zc*ys), tx
		, sx*yc*zs,  sy*(xc*zc + xs*ys*zs), -sz*(zc*xs - xc*ys*zs), ty
		, -sx*ys, sy*yc*xs, sz*xc*yc, tz
		, 0, 0, 0, 1};
	return res;
}

template<class T>
inline T lerp(const T& a1, const T& a2, const float t)
{
	return a1 + (a2 - a1) * t;
}

template<class T>
inline SBMatrix<T, 2, 1> deCasteljau(const SBMatrix<T, 2, 1>& p0, const SBMatrix<T, 2, 1>& p1) {
	return (p0+p1)*0.5f;
}

template<class T>
inline SBMatrix<T, 2, 1> quadBezier(const SBMatrix<T, 2, 1>& p0, const SBMatrix<T, 2, 1>& p1, const SBMatrix<T, 2, 1>& p2, const float s)
{
	const float u = 1 - s;
	return p0*(u*u)+p1*(2*s*u)+p2*(s*s);
}

template<class T>
inline SBMatrix<T, 2, 1> quadBezierDerivative(const SBMatrix<T, 2, 1>& p0, const SBMatrix<T, 2, 1>& p1, const SBMatrix<T, 2, 1>& p2, const float s)
{
	const SBMatrix<T, 2, 1> d0=p1-p0, d1=p2-p1;
	const bool b0 = isZero(d0), b1 = isZero(d1);
	if (b0 && b1) {
		SBMatrix<T, 2, 1> v;
		v.zeros();
		assert(false);
		return v;
	}
	else if (b0) {
		return d1;
	}
	else if (b1) {
		return d0;
	}

	assert(!b0);
	assert(!b1);
	SBMatrix<T, 2, 1> v = (d0*(1-s)+d1*s)*2;
	return normalize(v);
}

template<class T>
inline SBMatrix<T, 2, 1> cubicBezier(const SBMatrix<T, 2, 1>& p0, const SBMatrix<T, 2, 1>& p1, const SBMatrix<T, 2, 1>& p2, const SBMatrix<T, 2, 1>& p3, const float s)
{
	const float u = 1 - s;
	return p0*(u*u*u)+p1*(3*s*u*u)+p2*(3*s*s*u)+p3*(s*s*s);
}

template<class T>
inline SBMatrix<T, 2, 1> cubicBezierDerivative(const SBMatrix<T, 2, 1>& p0, const SBMatrix<T, 2, 1>& p1, const SBMatrix<T, 2, 1>& p2, const SBMatrix<T, 2, 1>& p3, const float s)
{
	const SBMatrix<T, 2, 1> d0=p1-p0,d1=p2-p1,d2=p3-p2;
	const bool b0 = isZero(d0), b1 = isZero(d1), b2 = isZero(d2);

	if (b0 && b1 && b2) {
		SBMatrix<T, 2, 1> v;
		v.zeros();
		assert(false);
		return v;
	}
	else if (b0 && b1) {
		return d2;
	}
	else if (b0 && b2) {
		return d1;
	}
	else if (b1 && b2) {
		return d0;
	}
	else if (b0) {
		return quadBezierDerivative(p1, p2, p3, s);
	}
	else if (b1) {
		return quadBezierDerivative(p0, p2, p3, s);
	}
	else if (b2) {
		return quadBezierDerivative(p0, p1, p2, s);
	}
	
	assert(!b0);
	assert(!b1);
	assert(!b2);
	const float u = 1 - s;
	SBMatrix<T, 2, 1> v = (d0*(u*u)+d1*(2*s*u)+d2*(s*s))*3;
	return normalize(v);
}

template<class T>
inline SBMatrix<T, 2, 1> arc(const T& a0, const T& a1, const SBMatrix<T, 3, 3>& m, const float s)
{
	const T a = s*(a1-a0)+a0;
	SBMatrix<T, 3, 1> v0 = {cos(a),sin(a),1.0};
	v0 = m * v0;
	SBMatrix<T, 2, 1> v1 = {v0.x, v0.y};
	return v1;
}

template<class T>
inline SBMatrix<T, 2, 1> arc(const T& rh, const T& rv, const T& rad, const T& cx, const T& cy, const float t)
{
	// [cr -sr;sr cr][rh 0;0 rv][ct;st]+[cx;cy]
	// [cr*rh -sr*rv;sr*rh cr*rv][ct;st]+[cx;cy]
	// [cr*rh*ct-sr*rv*st+cx;sr*rh*ct+cr*rv*st+cy]
	const T cr = cos(rad), sr = sin(rad);
	const T ct = cos(t), st = sin(t);
	const T rhct=rh*ct, rvst=rv*st;
	//SBMatrix<T, 2, 1> r = {cr*rh*ct-sr*rv*st+cx, sr*rh*ct+cr*rv*st+cy};
	SBMatrix<T, 2, 1> r = {cr*rhct-sr*rvst+cx, sr*rhct+cr*rvst+cy};
	return r;
}

template<class T>
inline SBMatrix<T, 2, 1> arcDerivative(const T& rh, const T& rv, const T& rad, const bool ccw, const float t)
{
	const SBMatrix<T, 2, 1> u0 = {cos(t), sin(t)};
	const SBMatrix<T, 2, 1> u1 = ccw ? perpCCW(u0) : perpCW(u0);
	const SBMatrix<T, 2, 1> u = rotate2(rad)*scale2(rh, rv)*u1;
	return normalize(u);
}

template<class T>
inline SBMatrix<T, 2, 1> circle(const T& radius, const T& cx, const T& cy, const float t)
{
	const SBMatrix<T, 2, 1> r = {
		cx+radius*cos(t), 
		cy+radius*sin(t)
	};
	return r;
}

template<class T>
inline SBMatrix<T, 2, 1> circleDerivative(const T& radius, const bool ccw, const float t)
{
	const SBMatrix<T, 2, 1> u0 = {cos(t), sin(t)};

	return ccw ? perpCCW(u0) : perpCW(u0);
}

template<class T>
inline T hermite(const T& p0, const T& p1, const T& t0, const T& t1, const float s)
{
	const float s2 = s * s, s3 = s2 * s;
	return (2*s3-3*s2+1)*p0+(s3-2*s2+s)*t0+(-2*s3+3*s2)*p1+(s3-s2)*t1;
}

template<class T>
inline void assign(const SBMatrix<T, 3, 1>& a, SBMatrix<T, 4, 1>* p)
{
	p->m[0][0] = a.m[0][0];	p->m[1][0] = a.m[1][0];	p->m[2][0] = a.m[2][0];	p->m[3][0] = 1;
}

template<class T>
inline void assign(const SBMatrix<T, 4, 1>& a, SBMatrix<T, 3, 1>* p)
{
	p->m[0][0] = a.m[0][0] / a.m[3][0];	p->m[1][0] = a.m[1][0] / a.m[3][0];	p->m[2][0] = a.m[2][0] / a.m[3][0];
}

template<class T>
inline SBMatrix<T, 2, 1> perpCW(const SBMatrix<T, 2, 1>& a)
{
	const SBMatrix<T, 2, 1> v = {a.m[1][0], -a.m[0][0]};
	return v;
}

template<class T>
inline SBMatrix<T, 2, 1> perpCCW(const SBMatrix<T, 2, 1>& a)
{
	const SBMatrix<T, 2, 1> v = {-a.m[1][0], a.m[0][0]};
	return v;
}

template<class T>
inline SBMatrix<T, 3, 3> ellipseMatrix(const T rh, const T rv, const T rad, const T cx, const T cy)
{
	const T c = cos(rad), s = sin(rad);
	const SBMatrix<T, 3, 3> m = {
		rh*c, -rv*s, cx,
		rh*s, rv*c, cy,
		0, 0, 1,
	};
	return m;
}

template<class T>
inline SBMatrix<T, 3, 3> ellipseInvMatrix(const T rh, const T rv, const T rad, const T cx, const T cy)
{
	const T c = cos(rad), s = sin(rad);
	const SBMatrix<T, 3, 3> m = {
		c/rh, s/rh, cx*-c/rh-cy*s/rh,
		-s/rv, c/rv, cx*s/rv-cy*c/rv,
		0, 0, 1,
	};
	return m;
}

template<class T>
inline T calcContinuousAngleCCW(const T t0, const T t1)
{
	// % CCW -> runs on positive direction
	// % t0 has to be smaller than t1

	if (t1 < t0)
		return t1+(T)(M_PI*2);

	return t1;
}

template<class T>
inline T calcContinuousAngleCW(const T t0, const T t1)
{
	// % CW -> runs on negative direction
	// % t0 has to be greater than t1

	if (t0 < t1)
		return t1-(T)(M_PI*2);

	return t1;
}

template<class T>
inline bool findUnitCircles(const SBMatrix<T, 2, 1>& v0, const SBMatrix<T, 2, 1>& v1, SBMatrix<T, 2, 1>* c0, SBMatrix<T, 2, 1>* c1)
{
	const SBMatrix<T, 2, 1> d = v0-v1;
	const SBMatrix<T, 2, 1> m = (v0+v1)*0.5;
	T dq = dot(d,d);
	if (!isZero(dq))
	{
		if (abs(dq-4)<0.0001f) {
			dq = 4;
		}
		const T rq = 1/dq-0.25f;
		if (rq >= 0)
		{
			const SBMatrix<T, 2, 1> r = d*sqrt(rq);
			*c0 = m+perpCW(r);
			*c1 = m+perpCCW(r);
			return true;
		}
	}
	return false;
}

template<class T>
inline bool findEllipses(const T rh, const T rv, const T rad, const SBMatrix<T, 2, 1>& v0, const SBMatrix<T, 2, 1>& v1, SBMatrix<T, 2, 1>* c0, SBMatrix<T, 2, 1>* c1)
{
	const T c = cos(rad), s= sin(rad);

	const SBMatrix<T, 2, 2> irs = {
		c/rh, s/rv, -s/rh, c/rv
	};
	const SBMatrix<T, 2, 1> t0=irs*v0;
	const SBMatrix<T, 2, 1> t1=irs*v1;

	SBMatrix<T, 2, 1> s0, s1;
	if (!findUnitCircles(t0, t1, &s0, &s1))
	{
		return false;
	}
	const SBMatrix<T, 2, 2> rs = {
		c*rh, -s*rv, s*rh, c*rv
	};
	*c0=rs*s0;
	*c1=rs*s1;
	return true;
}

template<class T>
inline void calcEllipticArcPosition(const SBMatrix<T, 3, 3>& emat, const T theta, T* px, T* py) {
	const SBMatrix<T, 3, 1> v = {cos(theta), sin(theta), 1};
	const SBMatrix<T, 3, 1> r = emat * v;
	*px = r.m[0][0];	*py = r.m[1][0];
}

template<class T>
inline void calcEllipticArcDerivativeCCW(const SBMatrix<T, 3, 3>& emat, const T theta, T* px, T* py) {
	const SBMatrix<T, 2, 1> e = {cos(theta), sin(theta)};
	const SBMatrix<T, 2, 1> d = perpCCW(e);
	*px = emat.m[0][0]*d.m[0][0]+emat.m[0][1]*d.m[1][0];
	*py = emat.m[1][0]*d.m[0][0]+emat.m[1][1]*d.m[1][0];
}

template<class T>
inline void calcEllipticArcDerivativeCW(const SBMatrix<T, 3, 3>& emat, const T theta, T* px, T* py) {
	const SBMatrix<T, 2, 1> e = {cos(theta), sin(theta)};
	const SBMatrix<T, 2, 1> d = perpCW(e);
	*px = emat.m[0][0]*d.m[0][0]+emat.m[0][1]*d.m[1][0];
	*py = emat.m[1][0]*d.m[0][0]+emat.m[1][1]*d.m[1][0];
}

template<class T>
inline bool calcArgOfQuadAt(const T& p0, const T& p1, const T& p2, const T& at, T* r) {
	const T a = p0-2*p1+p2, b = 2*(p1-p0), c = p0-at;
	T t1, t2;
	roots(a, b, c, &t1, &t2);

	if (t1 >= 0 && t1 <= 1) {
		*r = t1;
		return true;
	}
	if (t2 >= 0 && t2 <= 1) {
		*r = t2;
		return true;
	}
	return false;
}

template<class T>
inline bool calcArgOfQuadAt(const SBMatrix<T, 2, 1>& p0, const SBMatrix<T, 2, 1>& p1, const SBMatrix<T, 2, 1>& p2, const SBMatrix<T, 2, 1>& p3, const SBMatrix<T, 2, 1>& at, T* r) {
	if (calcArgOfCubicAt(p0.x, p1.x, p2.x, p3.x, at.x, r)) {
		return true;
	}
	return calcArgOfCubicAt(p0.y, p1.y, p2.y, p3.y, at.y, r);
}

template<class T>
inline bool calcArgOfCubicAt(const T& p0, const T& p1, const T& p2, const T& p3, const T& at, T* r) {
	const T a = -p0+3*p1-3*p2+p3, b = 3*(p0-2*p1+p2), c = 3*(p1-p0), d = p0-at;
	T t1, t2, t3;
	if (!roots(a, b, c, d, &t1, &t2, &t3)) {
		return false;
	}
	if (t1 >= 0 && t1 <= 1) {
		*r = t1;
		return true;
	}
	if (t2 >= 0 && t2 <= 1) {
		*r = t2;
		return true;
	}
	if (t3 >= 0 && t3 <= 1) {
		*r = t3;
		return true;
	}
	return false;
}

template<class T>
inline bool calcArgOfCubicAt(const SBMatrix<T, 2, 1>& p0, const SBMatrix<T, 2, 1>& p1, const SBMatrix<T, 2, 1>& p2, const SBMatrix<T, 2, 1>& p3, const SBMatrix<T, 2, 1>& at, T* r) {
	if (calcArgOfCubicAt(p0.x, p1.x, p2.x, p3.x, at.x, r)) {
		return true;
	}
	return calcArgOfCubicAt(p0.y, p1.y, p2.y, p3.y, at.y, r);
}

template<class T>
inline T calcArgOfArcAt(const SBMatrix<T, 3, 3>& imat, const T t0, const T t1, const SBMatrix<T, 2, 1>& at) {
	SBMatrix<T, 2, 1> iat = mul(imat, at);
	const T t = atan2(iat.y, iat.x);
	//t=t0+s(t1-t0)
	// s = (t-t0)/(t1-t0)
	const T s = (t-t0)/(t1-t0);
	return s;
}

template<class T>
T calcAASignCCW(const SBMatrix<T, 2, 1>& a, const SBMatrix<T, 2, 1>& b)
{
	return area(a, b)<0.0?-1.0:1.0;
}

template<class T>
T calcAASignCW(const SBMatrix<T, 2, 1>& a, const SBMatrix<T, 2, 1>& b)
{
	return area(a, b)<0.0?1.0:-1.0;
}

template<class T>
inline SBMatrix<T, 2, 1> calcAADirCCW(const SBMatrix<T, 2, 1>& d0, const SBMatrix<T, 2, 1>& d1)
{
	const SBMatrix<T, 2, 1> d0n = {-d0.x, -d0.y};
	const T d = dot(d0n, d1);
	if (d < -0.9921875) {	// d0, d1
		return perpCW(d0);
	}
	else if (d > 0.9921875) {
		return d0;
	}
	const T s = calcAASignCCW(d0n, d1);
	return normalize(d0n+d1)*(s*sqrt(2.0/(1.0-d)));
}

template<class T>
inline SBMatrix<T, 2, 1> calcAADirCW(const SBMatrix<T, 2, 1>& d0, const SBMatrix<T, 2, 1>& d1)
{
	const SBMatrix<T, 2, 1> d0n = {-d0.x, -d0.y};
	const T d = dot(d0n, d1);
	if (d < -0.9921875) {	// d0, d1
		return perpCCW(d0);
	}
	else if (d > 0.9921875) {
		return d0;
	}
	const T s = calcAASignCW(d0n, d1);
	return normalize(d0n+d1)*(s*sqrt(2.0/(1.0-d)));
}

template<class T>
inline void getAngle(const T bx, const T by, T* result)
{
	*result = fmod(2*M_PI + (by > 0.0 ? 1.0 : -1.0) * acos( bx / sqrt(bx * bx + by * by) ), 2*M_PI);
}

template<class T>	
inline void calArcBoundBox(const SBMatrix<T, 2, 1> p0, const T rh, const T rv, const SBMatrix<T, 2, 1> center, const T phi, const bool largeArc, const bool sweep, const SBMatrix<T, 2, 1> p1, T &xmin, T &ymin, T &xmax, T &ymax)
{
	T rx = rh;
	T ry = rv;
	T cx = center.x;
	T cy = center.y;
	T x1 = p0.x;
	T y1 = p0.y;
	T x2 = p1.x;
	T y2 = p1.y;

	if (rx <= 0.0 || ry <= 0.0)
	{
		xmin = (x1 < x2 ? x1 : x2);
		xmax = (x1 > x2 ? x1 : x2);
		ymin = (y1 < y2 ? y1 : y2);
		ymax = (y1 > y2 ? y1 : y2);
		return;
	}

	T txmin, txmax, tymin, tymax;

	if (phi == 0 || phi == M_PI)
	{
		xmin = cx - rx;
		getAngle(-rx, (T)0, &txmin);
		xmax = cx + rx;
		getAngle(rx, (T)0, &txmax);
		ymin = cy - ry;
		getAngle((T)0, -ry, &tymin);
		ymax = cy + ry;
		getAngle((T)0, ry, &tymax);
	} else if (phi == M_PI / 2.0 || phi == 3.0*M_PI/2.0)
	{
		xmin = cx - ry;
		getAngle(-ry, (T)0, &txmin);
		xmax = cx + ry;
		getAngle(ry, (T)0, &txmax);
		ymin = cy - rx;
		getAngle((T)0, -rx, &tymin);
		ymax = cy + rx;
		getAngle((T)0, rx, &tymax);
	}  else {
		txmin = -atan(ry*tan(phi)/rx);
		txmax = M_PI - atan (ry*tan(phi)/rx);
		xmin = cx + rx*cos(txmin)*cos(phi) - ry*sin(txmin)*sin(phi);
		xmax = cx + rx*cos(txmax)*cos(phi) - ry*sin(txmax)*sin(phi);
		if (xmin > xmax)
		{
			std::swap(xmin,xmax);
			std::swap(txmin,txmax);
		}
		T tmpY = cy + rx*cos(txmin)*sin(phi) + ry*sin(txmin)*cos(phi);
		getAngle(xmin - cx, tmpY - cy, &txmin);
		tmpY = cy + rx*cos(txmax)*sin(phi) + ry*sin(txmax)*cos(phi);
		getAngle(xmax - cx, tmpY - cy, &txmax);


		tymin = atan(ry/(tan(phi)*rx));
		tymax = atan(ry/(tan(phi)*rx))+M_PI;
		ymin = cy + rx*cos(tymin)*sin(phi) + ry*sin(tymin)*cos(phi);
		ymax = cy + rx*cos(tymax)*sin(phi) + ry*sin(tymax)*cos(phi);
		if (ymin > ymax)
		{
			std::swap(ymin,ymax);
			std::swap(tymin,tymax);
		}
		T tmpX = cx + rx*cos(tymin)*cos(phi) - ry*sin(tymin)*sin(phi);
		getAngle(tmpX - cx, ymin - cy, &tymin);
		tmpX = cx + rx*cos(tymax)*cos(phi) - ry*sin(tymax)*sin(phi);
		getAngle(tmpX - cx, ymax - cy, &tymax);
	}

	T angle1, angle2;
	getAngle(x1 - cx, y1 - cy, &angle1);
	getAngle(x2 - cx, y2 - cy, &angle2);

	if (!sweep)
		std::swap(angle1, angle2);

	bool otherArc = false;
	if (angle1 > angle2)
	{
		std::swap(angle1, angle2);
		otherArc = true;
	}

	if ((!otherArc && (angle1 > txmin || angle2 < txmin)) || (otherArc && !(angle1 > txmin || angle2 < txmin)))
		xmin = x1 < x2 ? x1 : x2;
	if ((!otherArc && (angle1 > txmax || angle2 < txmax)) || (otherArc && !(angle1 > txmax || angle2 < txmax)))
		xmax = x1 > x2 ? x1 : x2;
	if ((!otherArc && (angle1 > tymin || angle2 < tymin)) || (otherArc && !(angle1 > tymin || angle2 < tymin)))
		ymin = y1 < y2 ? y1 : y2;
	if ((!otherArc && (angle1 > tymax || angle2 < tymax)) || (otherArc && !(angle1 > tymax || angle2 < tymax)))
		ymax = y1 > y2 ? y1 : y2;
}

template<class T>
inline T calcDet(const SBMatrix<T, 2, 1>& v0, const SBMatrix<T, 2, 1>& v1)
{
	return v0.x * v1.y - v1.x * v0.y;
}

template<class T>
inline T calcDet(const SBMatrix<T, 2, 1>& v0, const SBMatrix<T, 2, 1>& v1, const SBMatrix<T, 2, 1>& v2)
{
	const SBMatrix<T, 2, 1> u0 = v1-v0;
	const SBMatrix<T, 2, 1> u1 = v2-v1;
	return calcDet(u0, u1);
}

template<class T>
inline T calcCircleIntegral(const T r, const T x) {
	const T x2 = x*x;
	const T r2 = r*r;
	const T srx = sqrt(r2-x2);
	return 0.5f*(x*srx+r2*atan(x/srx));
}

inline unsigned int calcPowerOfTwo(const unsigned int a)
{
	static const unsigned int ptwo[] = {
		0x1<<0x00,0x1<<0x01,0x1<<0x02,0x1<<0x03,
		0x1<<0x04,0x1<<0x05,0x1<<0x06,0x1<<0x07,
		0x1<<0x08,0x1<<0x09,0x1<<0x0a,0x1<<0x0b,
		0x1<<0x0c,0x1<<0x0d,0x1<<0x0e,0x1<<0x0f,
		0x1<<0x10,0x1<<0x11,0x1<<0x12,0x1<<0x13,
		0x1<<0x14,0x1<<0x15,0x1<<0x16,0x1<<0x17,
		0x1<<0x18,0x1<<0x19,0x1<<0x1a,0x1<<0x1b,
		0x1<<0x1c,0x1<<0x1d,0x1<<0x1e,//0x1<<0x1f,
	};

	for (unsigned int x = 0; x < 31; ++x) {	// integer -> 32 bits
		if (a <= ptwo[x]) {
			return ptwo[x];
		}
	}

	return 0;
}

inline bool isPowerOfTwo(const unsigned int a)
{
	static const unsigned int ptwo[] = {
		0x1<<0x00,0x1<<0x01,0x1<<0x02,0x1<<0x03,
		0x1<<0x04,0x1<<0x05,0x1<<0x06,0x1<<0x07,
		0x1<<0x08,0x1<<0x09,0x1<<0x0a,0x1<<0x0b,
		0x1<<0x0c,0x1<<0x0d,0x1<<0x0e,0x1<<0x0f,
		0x1<<0x10,0x1<<0x11,0x1<<0x12,0x1<<0x13,
		0x1<<0x14,0x1<<0x15,0x1<<0x16,0x1<<0x17,
		0x1<<0x18,0x1<<0x19,0x1<<0x1a,0x1<<0x1b,
		0x1<<0x1c,0x1<<0x1d,0x1<<0x1e,//0x1<<0x1f,
	};

	for (unsigned int x = 0; x < 31; ++x) {	// integer -> 32 bits
		if (a == ptwo[x]) {
			return true;
		}
	}

	return false;
}

#endif
