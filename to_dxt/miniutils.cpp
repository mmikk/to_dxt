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

#include "miniutils.h"
#include <assert.h>





static unsigned int uSeed=39871946;

void QuickSort(uint32* pSortBuffer, int iLeft, int iRight)
{
	// Random
	unsigned int t=uSeed&31;
	t=(uSeed<<t)|(uSeed>>(32-t));
	uSeed=uSeed+t+3;
	// Random end

	int iL=iLeft, iR=iRight;
	int n = (iR-iL)+1;
	assert(n>=0);
	int index = (int) (uSeed%n);

	uint32 uMid=pSortBuffer[index + iL];
	uint32 uTmp;

	do
	{
		while(pSortBuffer[iL] < uMid)
			++iL;
		while(pSortBuffer[iR] > uMid)
			--iR;

		if(iL <= iR)
		{
			uTmp = pSortBuffer[iL];
			pSortBuffer[iL] = pSortBuffer[iR];
			pSortBuffer[iR] = uTmp;
			++iL; --iR;
		}
	}
	while(iL <= iR);

	if(iLeft < iR)
		QuickSort(pSortBuffer, iLeft, iR);
	if(iL < iRight)
		QuickSort(pSortBuffer, iL, iRight);
}

void QuickSort(uint32* pSortBuffer, uint32 * pDataBuffer, int iLeft, int iRight)
{
	// Random
	unsigned int t=uSeed&31;
	t=(uSeed<<t)|(uSeed>>(32-t));
	uSeed=uSeed+t+3;
	// Random end

	int iL=iLeft, iR=iRight;
	int n = (iR-iL)+1;
	assert(n>=0);
	int index = (int) (uSeed%n);

	uint32 uMid=pSortBuffer[index + iL];
	uint32 uTmp;

	do
	{
		while(pSortBuffer[iL] < uMid)
			++iL;
		while(pSortBuffer[iR] > uMid)
			--iR;

		if(iL <= iR)
		{
			uTmp = pSortBuffer[iL];
			pSortBuffer[iL] = pSortBuffer[iR];
			pSortBuffer[iR] = uTmp;

			uTmp = pDataBuffer[iL];
			pDataBuffer[iL] = pDataBuffer[iR];
			pDataBuffer[iR] = uTmp;

			++iL; --iR;
		}
	}
	while(iL <= iR);

	if(iLeft < iR)
		QuickSort(pSortBuffer, iLeft, iR);
	if(iL < iRight)
		QuickSort(pSortBuffer, iL, iRight);
}


void	LoadIdentity(Mat33 * pM)
{
	pM->m_fMat[0+0*3] = 1;
	pM->m_fMat[1+0*3] = 0;
	pM->m_fMat[2+0*3] = 0;
	pM->m_fMat[0+1*3] = 0;
	pM->m_fMat[1+1*3] = 1;
	pM->m_fMat[2+1*3] = 0;
	pM->m_fMat[0+2*3] = 0;
	pM->m_fMat[1+2*3] = 0;
	pM->m_fMat[2+2*3] = 1;
}

const Mat33	Transpose( const Mat33 &m )
{
	Mat33 r;

	r.m_fMat[0+0*3] = m.m_fMat[0+0*3];
	r.m_fMat[1+0*3] = m.m_fMat[0+1*3];
	r.m_fMat[2+0*3] = m.m_fMat[0+2*3];

	r.m_fMat[0+1*3] = m.m_fMat[1+0*3];
	r.m_fMat[1+1*3] = m.m_fMat[1+1*3];
	r.m_fMat[2+1*3] = m.m_fMat[1+2*3];

	r.m_fMat[0+2*3] = m.m_fMat[2+0*3];
	r.m_fMat[1+2*3] = m.m_fMat[2+1*3];
	r.m_fMat[2+2*3] = m.m_fMat[2+2*3];

	return r;
}

float		Determinant( const Mat33 &m )
{
	return m.m_fMat[0+0*3]*(m.m_fMat[1+1*3]*m.m_fMat[2+2*3] - m.m_fMat[1+2*3]*m.m_fMat[2+1*3]) -
		   m.m_fMat[0+1*3]*(m.m_fMat[1+0*3]*m.m_fMat[2+2*3] - m.m_fMat[1+2*3]*m.m_fMat[2+0*3]) +
		   m.m_fMat[0+2*3]*(m.m_fMat[1+0*3]*m.m_fMat[2+1*3] - m.m_fMat[1+1*3]*m.m_fMat[2+0*3]);
}

Mat33::Mat33( const Mat33 &m )
{
	m_fMat[0+0*3] = m.m_fMat[0+0*3];
	m_fMat[1+0*3] = m.m_fMat[1+0*3];
	m_fMat[2+0*3] = m.m_fMat[2+0*3];
	m_fMat[0+1*3] = m.m_fMat[0+1*3];
	m_fMat[1+1*3] = m.m_fMat[1+1*3];
	m_fMat[2+1*3] = m.m_fMat[2+1*3];
	m_fMat[0+2*3] = m.m_fMat[0+2*3];
	m_fMat[1+2*3] = m.m_fMat[1+2*3];
	m_fMat[2+2*3] = m.m_fMat[2+2*3];
}