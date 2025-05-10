/**
 *  Copyright (C) 2011 by Morten S. Mikkelsen
 *
 *  This software is provided 'as-is', without any express or implied
 *  warranty.  In no event will the authors be held liable for any damages
 *  arising from the use of this software.
 *
 *  Permission is granted to anyone to use this software for any purpose,
 *  including commercial applications, and to alter it and redistribute it
 *  freely, subject to the following restrictions:
 *
 *  1. The origin of this software must not be misrepresented; you must not
 *     claim that you wrote the original software. If you use this software
 *     in a product, an acknowledgment in the product documentation would be
 *     appreciated but is not required.
 *  2. Altered source versions must be plainly marked as such, and must not be
 *     misrepresented as being the original software.
 *  3. This notice may not be removed or altered from any source distribution.
 */

#ifndef __MINIUTILS_H__
#define __MINIUTILS_H__

#include <math.h>

typedef unsigned char uint8;
typedef char int8;
typedef int int32;
typedef long long int int64;
typedef unsigned long long int uint64;
typedef unsigned int uint32;


void QuickSort(uint32* pSortBuffer, int iLeft, int iRight);
void QuickSort(uint32* pSortBuffer, uint32 * pDataBuffer, int iLeft, int iRight);

inline const int Abs(const int iX) { return iX<0?(-iX):iX; }

// vector 3 stuff
struct Vec3
{
	Vec3( const Vec3 &v ) : x(v.x), y(v.y), z(v.z) {}
	Vec3( float fX, float fY, float fZ ) : x(fX), y(fY), z(fZ) {}
	Vec3() : x(0.0f), y(0.0f), z(0.0f) {}


	float x, y, z;
};

inline const float operator *( const Vec3 &v1, const Vec3 &v2)
{
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

inline const Vec3 operator *(const float fS, const Vec3 &v)
{
	Vec3 vRes;

	vRes.x = fS * v.x;
	vRes.y = fS * v.y;
	vRes.z = fS * v.z;

	return vRes;
}

inline const Vec3 operator *(const Vec3 &v, const float fS)
{
	return fS * v;
}


inline const Vec3 Cross(const Vec3 &v1, const Vec3 &v2)
{
	Vec3 v;
	
	v.x = v1.y*v2.z - v2.y*v1.z;
	v.y = v1.z*v2.x - v2.z*v1.x;
	v.z = v1.x*v2.y - v2.x*v1.y;
	
	return v;
}

inline const Vec3 operator +( const Vec3 &v1, const Vec3 &v2 )
{
	Vec3 vRes;

	vRes.x = v1.x + v2.x;
	vRes.y = v1.y + v2.y;
	vRes.z = v1.z + v2.z;

	return vRes;
}

inline Vec3& operator +=(Vec3 &v1, const Vec3 &v2)
{
	v1.x += v2.x;
	v1.y += v2.y;
	v1.z += v2.z;

	return v1;
}


inline const Vec3 operator -( const Vec3 &v1, const Vec3 &v2 )
{
	Vec3 vRes;

	vRes.x = v1.x - v2.x;
	vRes.y = v1.y - v2.y;
	vRes.z = v1.z - v2.z;

	return vRes;
}


inline Vec3& operator *=(Vec3 &v, const float fS )
{
	v.x *= fS;
	v.y *= fS;
	v.z *= fS;

	return v;
}



inline float LengthSquared( const Vec3 &v )
{
	return v.x*v.x + v.y*v.y + v.z*v.z;
}

inline float Length( const Vec3 &v )
{
	return (float) sqrt(LengthSquared(v));	// this doesn't need to be sqrtf
}

inline const Vec3		Normalize( const Vec3 &v )
{
	return (1 / Length(v)) * v;
}

inline bool	operator ==( const Vec3 &v1, const Vec3 &v2 )
{
	return (v1.x == v2.x) && (v1.y == v2.y) && (v1.z == v2.z);
}









// matrix 3x3 stuff
struct Mat33
{
	Mat33() 
	{
		for(int j=0; j<3; j++)
		{
			for(int i=0; i<3; i++)
				m_fMat[i + j*3] = i!=j ? 0.0f : 1.0f;
		}
	}
	Mat33( const Mat33 &m );


	float m_fMat[3*3];
};

inline void	SetRow(Mat33 * pM, int iRow, const Vec3 &v )
{
	pM->m_fMat[iRow+0*3] = v.x;
	pM->m_fMat[iRow+1*3] = v.y;
	pM->m_fMat[iRow+2*3] = v.z;
}

inline void	SetColumn(Mat33 * pM, int iColumn, const Vec3 &v )
{
	pM->m_fMat[0+iColumn*3] = v.x;
	pM->m_fMat[1+iColumn*3] = v.y;
	pM->m_fMat[2+iColumn*3] = v.z;
}

inline const Vec3	operator *( const Mat33 &m, const Vec3 &v )
{
	Vec3 r;

	r.x = m.m_fMat[0+0*3]*v.x + m.m_fMat[0+1*3]*v.y + m.m_fMat[0+2*3]*v.z;
	r.y = m.m_fMat[1+0*3]*v.x + m.m_fMat[1+1*3]*v.y + m.m_fMat[1+2*3]*v.z;
	r.z = m.m_fMat[2+0*3]*v.x + m.m_fMat[2+1*3]*v.y + m.m_fMat[2+2*3]*v.z;

	return r;
}

void	LoadIdentity(Mat33 * pM);
const Mat33	Transpose( const Mat33 &m );
float		Determinant( const Mat33 &m );


#endif