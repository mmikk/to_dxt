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

#include "encoder_dxt.h"
#include <assert.h>
#include <string.h>

// this scalar 2*3 is used to maintain
// the interpolated dxt palette in integer form.
// this will allow us to achieve order
// independant accumulation of error
#define TWO_X_THREE		6


void FindBestTrans(const Vec3 vVerts[], const int nr_verts, Mat33 * pRot, Vec3 * pPos);

float fmaxf(const float a, const float b) { return a > b ? a : b; }

int Min(const int32 A, const int32 B) { return A<B?A:B; }
int Max(const int32 A, const int32 B) { return A>B?A:B; }
int iAbs(const int32 iX) { return iX < 0 ? (-iX) : iX; }

struct SDXT_RGB_Pal;

int ColorAs888(const Vec3 &v0);
void ComputeColorsAs565(const Vec3 &v0, const Vec3 &v1, int * pColor0, int * pColor1);
void ComputeColorAs565(const Vec3 &v0, int * pColor0);
void GetColors(SDXT_RGB_Pal * pRGB_Pal, const int32 c0, const int32 c1);
void CompareResult(const uint32 block[], const int iBlockColors, const SDXT_RGB_Pal &rgb_pal, const int32 color0, const int32 color1, int * piCurDist, int * piBestC0, int * piBestC1);

// keeping unsigned data in signed integers (top bit reserved), top bit must remain zero.
bool IsValidAdd32bit(const int32 iA, const int32 iB) { return ((int64) iA) < (0x7fffffff-iB); }

// these need to work for floats too. So no ANDing, ORing or whatever here.
//#define GPU_8_TO_5( X_8BIT ) ( (X_8BIT*0x1f + 0x7f) / 0xff )
//#define GPU_5_TO_8( X_5BIT ) ( (X_5BIT*0xff + 0xf) / 0x1f )
//#define GPU_8_TO_6( X_8BIT ) ( (X_8BIT*0x3f + 0x7f) / 0xff )
//#define GPU_6_TO_8( X_6BIT ) ( (X_6BIT*0xff + 0x1f) / 0x3f )

// according to the d3d 10 specification (which is not public).
int d3d_spec_5_to_8(int iValue) { return ((iValue<<3) | ((iValue>>2)&0x7))&0xff; }
int d3d_spec_6_to_8(int iValue) { return ((iValue<<2) | ((iValue>>4)&0x3))&0xff; }

int d3d_spec_8_to_5(int iValue)
{
	static bool bTableNotInitialized = true;
	static uint8 uBest_5bit_rep[256];
	if(bTableNotInitialized)
	{
		for(int i=0; i<256; i++)
		{
			int iBestValue = 0;
			int iMinError = iAbs(i-d3d_spec_5_to_8(0));
			for(int j=1; j<0x20; j++)
			{
				const int iCurValue = j;
				const int iCurError = iAbs(i-d3d_spec_5_to_8(j));
				if(iCurError < iMinError)
				{
					iMinError = iCurError;
					iBestValue = iCurValue;
				}
			}

			uBest_5bit_rep[i] = iBestValue;
		}

		bTableNotInitialized = false;
	}

	// result
	return uBest_5bit_rep[iValue&0xff];
}

int d3d_spec_8_to_6(int iValue)
{
	static bool bTableNotInitialized = true;
	static uint8 uBest_6bit_rep[256];
	if(bTableNotInitialized)
	{
		for(int i=0; i<256; i++)
		{
			int iBestValue = 0;
			int iMinError = iAbs(i-d3d_spec_6_to_8(0));
			for(int j=1; j<0x40; j++)
			{
				const int iCurValue = j;
				const int iCurError = iAbs(i-d3d_spec_6_to_8(j));
				if(iCurError < iMinError)
				{
					iMinError = iCurError;
					iBestValue = iCurValue;
				}
			}

			uBest_6bit_rep[i] = iBestValue;
		}

		bTableNotInitialized = false;
	}

	// result
	return uBest_6bit_rep[iValue&0xff];
}



int GetC565_from_C888(const int color);
int GetC888_from_C565(const int color);

int Process3Fast(const uint32 uLocalCols[], const Vec3 vSortedVerts[], const int iNumCols, const bool bIsFourOnLine, const Vec3 &vFinPos, const Mat33 &mFinRot, int * pC0, int * pC1);

void CorrectColors(int * pC0, int * pC1, const bool bIsFourOnLine)
{
	if( ((pC0[0] > pC1[0]) ^ bIsFourOnLine) )
	{	const int iTmp = pC0[0]; pC0[0] = pC1[0]; pC1[0] = iTmp; }

	if( (pC0[0]==pC1[0]) && bIsFourOnLine )
	{
		if( pC0[0]==0 ) ++pC0[0];
		else if( ((pC1[0]>>0)&0x1f) != 0 ) pC1[0] -= (1<<0);
		else if( ((pC1[0]>>5)&0x3f) != 0 ) pC1[0] -= (1<<5);
		else if( ((pC1[0]>>11)&0x1f) != 0 ) pC1[0] -= (1<<11);
	}
	assert(((pC0[0] > pC1[0])&&bIsFourOnLine) || ((pC0[0] <= pC1[0])&&(!bIsFourOnLine)));
}


void FindBest5BitPair(const uint8 uValue, int * pA, int * pB, const bool bIsFourOnLine)
{
	static bool bTablesInitialized = false;
	static uint8 uBestI_5_IsFour[256];		// 4 x 0.25 kB
	static uint8 uBestJ_5_IsFour[256];
	static uint8 uBestI_5_IsNotFour[256];
	static uint8 uBestJ_5_IsNotFour[256];

	if(!bTablesInitialized)
	{
		for(int c=0; c<256; c++)
		{
			int iBestI_4, iBestJ_4, iSmallestError_4=0;
			int iBestI_3, iBestJ_3, iSmallestError_3=0;
			int iValue = c;
			for(int i=0; i<0x20; i++)
				for(int j=0; j<0x20; j++)
				{
					const int iNom_4 = d3d_spec_5_to_8(i)*2 + d3d_spec_5_to_8(j);
					const int iNom_3 = d3d_spec_5_to_8(i) + d3d_spec_5_to_8(j);
					const int iCurError_4 = iAbs( 3*iValue - iNom_4 );
					const int iCurError_3 = iAbs( 2*iValue - iNom_3 );

					if((i==0 && j==0) || iCurError_4<iSmallestError_4) { iBestI_4=i; iBestJ_4=j; iSmallestError_4 = iCurError_4; }
					if((i==0 && j==0) || iCurError_3<iSmallestError_3) { iBestI_3=i; iBestJ_3=j; iSmallestError_3 = iCurError_3; }
				}

			// fill up tables
			uBestI_5_IsFour[c] = iBestI_4; uBestJ_5_IsFour[c] = iBestJ_4;
			uBestI_5_IsNotFour[c] = iBestI_3; uBestJ_5_IsNotFour[c] = iBestJ_3;
		}

		bTablesInitialized = true;
	}

	if(bIsFourOnLine)
	{ pA[0] = uBestI_5_IsFour[uValue]; pB[0] = uBestJ_5_IsFour[uValue]; }
	else
	{ pA[0] = uBestI_5_IsNotFour[uValue]; pB[0] = uBestJ_5_IsNotFour[uValue]; }
}

void FindBest6BitPair(const uint8 uValue, int * pA, int * pB, const bool bIsFourOnLine)
{
	static bool bTablesInitialized = false;
	static uint8 uBestI_6_IsFour[256];		// 4 x 0.25 kB
	static uint8 uBestJ_6_IsFour[256];
	static uint8 uBestI_6_IsNotFour[256];
	static uint8 uBestJ_6_IsNotFour[256];

	if(!bTablesInitialized)
	{
		for(int c=0; c<256; c++)
		{
			int iBestI_4, iBestJ_4, iSmallestError_4=0;
			int iBestI_3, iBestJ_3, iSmallestError_3=0;
			int iValue = c;
			for(int i=0; i<0x40; i++)
				for(int j=0; j<0x40; j++)
				{
					const int iNom_4 = d3d_spec_6_to_8(i)*2 + d3d_spec_6_to_8(j);
					const int iNom_3 = d3d_spec_6_to_8(i) + d3d_spec_6_to_8(j);
					const int iCurError_4 = iAbs( 3*iValue - iNom_4 );
					const int iCurError_3 = iAbs( 2*iValue - iNom_3 );

					if((i==0 && j==0) || iCurError_4<iSmallestError_4) { iBestI_4=i; iBestJ_4=j; iSmallestError_4 = iCurError_4; }
					if((i==0 && j==0) || iCurError_3<iSmallestError_3) { iBestI_3=i; iBestJ_3=j; iSmallestError_3 = iCurError_3; }
				}

				// fill up tables
				uBestI_6_IsFour[c] = iBestI_4; uBestJ_6_IsFour[c] = iBestJ_4;
				uBestI_6_IsNotFour[c] = iBestI_3; uBestJ_6_IsNotFour[c] = iBestJ_3;
		}

		bTablesInitialized = true;
	}

	if(bIsFourOnLine)
	{ pA[0] = uBestI_6_IsFour[uValue]; pB[0] = uBestJ_6_IsFour[uValue]; }
	else
	{ pA[0] = uBestI_6_IsNotFour[uValue]; pB[0] = uBestJ_6_IsNotFour[uValue]; }
}

struct SDXT_RGB_Pal
{
	int iColors[4];
	int iPalMul;	// 2/3
};

int ProcessOneColor(const int32 iColor, int * pC0, int * pC1, const bool bIsFourOnLine)
{
	int iValidInputs = 2;
	int iInputs[6];

	int c0_565 = GetC565_from_C888(iColor);
	int c1_565 = c0_565;

	CorrectColors(&c0_565, &c1_565, bIsFourOnLine);

	iInputs[0] = c0_565;
	iInputs[1] = c1_565;

	int Ar, Br, Ag, Bg, Ab, Bb;
	FindBest5BitPair((iColor>>0)&0xff, &Ar, &Br, bIsFourOnLine);
	FindBest6BitPair((iColor>>8)&0xff, &Ag, &Bg, bIsFourOnLine);
	FindBest5BitPair((iColor>>16)&0xff, &Ab, &Bb, bIsFourOnLine);

	c0_565 = (Ar<<11) | (Ag<<5) | (Ab<<0);
	c1_565 = (Br<<11) | (Bg<<5) | (Bb<<0);

	CorrectColors(&c0_565, &c1_565, bIsFourOnLine);

	iInputs[2] = c0_565;
	iInputs[3] = c1_565;

	int iCurDist = -1;
	for(int i=0; i<iValidInputs; i++)
	{
		SDXT_RGB_Pal rgb_pal;
		const int c0_565 = iInputs[2*i+0];
		const int c1_565 = iInputs[2*i+1];
		GetColors(&rgb_pal, c0_565, c1_565);
		CompareResult((const uint32 *) &iColor, 1, rgb_pal, c0_565, c1_565, &iCurDist, pC0, pC1);
	}

	return iCurDist;
}


int ProceedMethod3(const uint32 uCols[], const int iNumCols, int * pC0, int * pC1, const bool bIsFourOnLine, const bool bApplyFastMethod)
{
	assert(iNumCols<=16);

	if(iNumCols == 0)
	{	pC0[0] = 0; pC1[0] = 0; return 0;	}

	int i=1;
	bool bFoundSecondColor = false;
	while(i<iNumCols && (!bFoundSecondColor))
		bFoundSecondColor = (uCols[0]&0x00ffffff)!=(uCols[i++]&0x00ffffff);

	if(!bFoundSecondColor)
		return iNumCols * ProcessOneColor((uCols[0]&0x00ffffff), pC0, pC1, bIsFourOnLine);	// mul error sum by iNumCols

	Vec3 vVerts[16];
	for(int k=0; k<iNumCols; k++)
		vVerts[k] = Vec3( (float) ((uCols[k]>>0)&0xff), (float) ((uCols[k]>>8)&0xff), (float) ((uCols[k]>>16)&0xff) );


	Mat33 rot;
	Vec3 pos;
	FindBestTrans(vVerts, iNumCols, &rot, &pos);


	// local space to block is pos + rot * vVertLocal
	// block to local is (vVert-pos) * rot

	// make the best fit line the X-axis in local space.
	const Vec3 vFinPos = pos;
	const Mat33 mFinRot = rot;

	const Mat33 mInvFinRot = Transpose(mFinRot);

	for(int k=0; k<iNumCols; k++)
		vVerts[k] = mInvFinRot * (vVerts[k]-vFinPos);

	// Sort based on vVerts[k].x (vVerts, uLocalCols)
	float fMin = 0;
	for(i=0; i<iNumCols; i++)
	{
		if(i==0) fMin = vVerts[i].x;
		else if(vVerts[i].x < fMin) fMin = vVerts[i].x;
	}
	uint32 uSortValues[16];	// make all floats unsigned for integer sorting
	uint32 uIndices[16];
	for(i=0; i<iNumCols; i++)
	{
		uIndices[i] = i;
		((float *) uSortValues)[i] = fmaxf(vVerts[i].x - fMin, 0);
		assert((uSortValues[i]&0x80000000) == 0);
	}
	if(iNumCols>1)
		QuickSort(uSortValues, uIndices, 0, iNumCols-1);

	Vec3 vSortedVerts[16];
	uint32 uLocalCols[16];
	for(i=0; i<iNumCols; i++)
	{
		uLocalCols[i] = uCols[uIndices[i]];
		vSortedVerts[i] = vVerts[uIndices[i]];
	}

	// early out.
	if(bApplyFastMethod) return Process3Fast(uLocalCols, vSortedVerts, iNumCols, bIsFourOnLine, vFinPos, mFinRot, pC0, pC1);

	// preprocess a lookup-table
	struct SSumsTable
	{
		double dSumX, dSumY, dSumZ;
		int iTwoSumR, iTwoSumG, iTwoSumB;
		int iSqSum;

		SSumsTable() {}

		SSumsTable(int i)
		{
			dSumX = 0; dSumY = 0; dSumZ = 0;
			iTwoSumR = 0; iTwoSumG = 0; iTwoSumB = 0; iSqSum = 0;
		}
	};

	SSumsTable sSumsTable[17][17];	// the total size of this tables is 11560 bytes (roughly 11.3kB)
									// TODO: verify lookups are generally hits and not misses in the cache.

	//memset(sSumsTable, 0, sizeof(sSumsTable));

	for(int x=0; x<=iNumCols; x++)
	{
		sSumsTable[x][x] = SSumsTable(0);	// Clear diagonal.

		for(int y=(x+1); y<=iNumCols; y++)
		{
			const float fVals[] = {vSortedVerts[y-1].x, vSortedVerts[y-1].y, vSortedVerts[y-1].z};
			const int iVals[] = {(int) ((uLocalCols[y-1]>>0)&0xff), (int) ((uLocalCols[y-1]>>8)&0xff), (int) ((uLocalCols[y-1]>>16)&0xff)};
			const int iLenSQ = iVals[0]*iVals[0]+iVals[1]*iVals[1]+iVals[2]*iVals[2];

			SSumsTable &cur_entry = sSumsTable[x][y];
			const SSumsTable &prev_entry = sSumsTable[x][y-1];

			cur_entry.dSumX = fVals[0] + prev_entry.dSumX;
			cur_entry.dSumY = fVals[1] + prev_entry.dSumY;
			cur_entry.dSumZ = fVals[2] + prev_entry.dSumZ;

			cur_entry.iTwoSumR = 2*iVals[0] + prev_entry.iTwoSumR;
			cur_entry.iTwoSumG = 2*iVals[1] + prev_entry.iTwoSumG;
			cur_entry.iTwoSumB = 2*iVals[2] + prev_entry.iTwoSumB;

			cur_entry.iSqSum = prev_entry.iSqSum + iLenSQ;
		}
	}

	// initialize error by center color
	// use SUMi_0_to_n (R[i] - iCenR)^2 = ( SUMi_0_to_n R[i]^2 ) + n*iCenR^2 -2*iCenR*( SUMi_0_to_n R[i])
	// and of course use the same property for green and blue.
	// 888 to 565
	int iRcen = (sSumsTable[0][iNumCols].iTwoSumR / (2*iNumCols))&0xff;
	int iGcen = (sSumsTable[0][iNumCols].iTwoSumG / (2*iNumCols))&0xff;
	int iBcen = (sSumsTable[0][iNumCols].iTwoSumB / (2*iNumCols))&0xff;
	const int single_color = GetC565_from_C888( ((iBcen<<16)|(iGcen<<8)|(iRcen<<0)) );

	iRcen = d3d_spec_5_to_8( d3d_spec_8_to_5( iRcen ) );
	iGcen = d3d_spec_6_to_8( d3d_spec_8_to_6( iGcen ) );
	iBcen = d3d_spec_5_to_8( d3d_spec_8_to_5( iBcen ) );
	int iSumOfSquaresError = (TWO_X_THREE*TWO_X_THREE)*( sSumsTable[0][iNumCols].iSqSum + iNumCols*(iRcen*iRcen+iGcen*iGcen+iBcen*iBcen) - (iRcen*sSumsTable[0][iNumCols].iTwoSumR + iGcen*sSumsTable[0][iNumCols].iTwoSumG + iBcen*sSumsTable[0][iNumCols].iTwoSumB) );
	
	*pC0 = single_color;
	*pC1 = single_color;

	// find the best match!
	const int tot_num = iNumCols;
	for(int m=0; m<=tot_num; m++)
		for(int n=0;n<=(tot_num-m); n++)
		{
			const int iPend = tot_num-m-n;
			const int iPstart = bIsFourOnLine ? 0 : iPend;

			for(int p=iPstart; p<=iPend; p++)
			{
				const int q = bIsFourOnLine ? (tot_num-m-n-p) : 0;

				int A = bIsFourOnLine ? (9*m + 4*n + p) : (4*m+n);
				int B = bIsFourOnLine ? (2*n + 2*p) : n;
				int C = bIsFourOnLine ? (9*q + 4*p + n) : (4*p+n);

				const int iDeterminant = A*C - B*B;

				if(iDeterminant != 0)
				{
					const int iGroup0_index = 0;	// and m elems from then on
					const int iGroup1_index = m;	// and n elems from then on
					const int iGroup2_index = m+n;	// and p elems from then on
					const int iGroup3_index = m+n+p;// and q elems from then on
					const int iGroup4_index = m+n+p+q;	// end of group3 (and not a real group)


					const SSumsTable &table0 = sSumsTable[iGroup0_index][iGroup1_index];
					const SSumsTable &table1 = sSumsTable[iGroup1_index][iGroup2_index];
					const SSumsTable &table2 = sSumsTable[iGroup2_index][iGroup3_index];
					const SSumsTable &table3 = sSumsTable[iGroup3_index][iGroup4_index];

					const double dSum1X = table0.dSumX;		// i_to_m
					const double dSum2X = table1.dSumX;		// j_to_n
					const double dSum3X = table2.dSumX;		// k_to_p
					const double dSum4X = table3.dSumX;		// l_to_q

					const double C1x = bIsFourOnLine ? (9*dSum1X  +  6*dSum2X  +  3*dSum3X) : (4*dSum1X + 2*dSum2X);
					const double C2x = bIsFourOnLine ? (9*dSum4X  +  6*dSum3X  +  3*dSum2X) : (4*dSum3X + 2*dSum2X);

					const double dSum1Y = table0.dSumY;
					const double dSum2Y = table1.dSumY;
					const double dSum3Y = table2.dSumY;
					const double dSum4Y = table3.dSumY;

					const double C1y = bIsFourOnLine ? (9*dSum1Y  +  6*dSum2Y  +  3*dSum3Y) : (4*dSum1Y + 2*dSum2Y);
					const double C2y = bIsFourOnLine ? (9*dSum4Y  +  6*dSum3Y  +  3*dSum2Y) : (4*dSum3Y + 2*dSum2Y);

					const double dSum1Z = table0.dSumZ;
					const double dSum2Z = table1.dSumZ;
					const double dSum3Z = table2.dSumZ;
					const double dSum4Z = table3.dSumZ;

					const double C1z = bIsFourOnLine ? (9*dSum1Z  +  6*dSum2Z  +  3*dSum3Z) : (4*dSum1Z + 2*dSum2Z);
					const double C2z = bIsFourOnLine ? (9*dSum4Z  +  6*dSum3Z  +  3*dSum2Z) : (4*dSum3Z + 2*dSum2Z);

					// 4 points on line
					// (9*m + 4*n + p)*c0.x + (2*n + 2*p)*c1.x  = 9*SUMi_to_m  (  C[i].x )  +  6*SUMj_to_n  (  C[j].x )  +  3*SUMk_to_p (  C[k].x )
					// (2*p+2*n)*c0.x +  (9*q + 4*p + n)*c1.x   = 9*SUMl_to_q  (  C[l].x )  +  6*SUMk_to_p  (  C[k].x )  +  3*SUMj_to_n (  C[j].x )

					// 3 points on line
					// (4*m+n)*c0.x + n*c1.x = 4 * (SUMi_to_m C[i].x) + 2 * (SUMj_to_n C[j].x)
					// n*c0.x + (4*p+n)*c1.x = 4 * (SUMk_to_p C[k].x) + 2 * (SUMj_to_n C[j].x)

					// solve equation
					// | A B | |x|   |C1|
					// | B C | |y| = |C2|

					// |x|			         | C -B | |C1|
					// |y| = (1/(A*C-B*B)) * |-B  A | |C2|

					const float x0 = (float) ( (C*C1x - B*C2x)/iDeterminant );
					const float x1 = (float) ( ((-B)*C1x + A*C2x)/iDeterminant );

					const float y0 = (float) ( (C*C1y - B*C2y)/iDeterminant );
					const float y1 = (float) ( ((-B)*C1y + A*C2y)/iDeterminant );

					const float z0 = (float) ( (C*C1z - B*C2z)/iDeterminant );
					const float z1 = (float) ( ((-B)*C1z + A*C2z)/iDeterminant );



					// evaluate error
					const Vec3 vC0 = vFinPos + mFinRot*Vec3(x0,y0,z0);
					const Vec3 vC1 = vFinPos + mFinRot*Vec3(x1,y1,z1);

					int color0, color1;
					ComputeColorsAs565(vC0, vC1, &color0, &color1);
					bool bDoFlip = false;
					if( ((color0 > color1) ^ bIsFourOnLine) )
					{	const int iTmp = color0; color0 = color1; color1 = iTmp; bDoFlip=true; }

					// special case: color0 == color1 and bIsFourOnLine
					// this will have to be allowed for now to compute a correct error value.
					// the code will correct it before returning from this function.
					assert(((color0 >= color1)&&bIsFourOnLine) || ((color0 <= color1)&&(!bIsFourOnLine)));

					SDXT_RGB_Pal rgb_pal;
					const int iPalEntries = bIsFourOnLine ? 4 : 3;
					GetColors(&rgb_pal, color0, color1);
					if(color0==color1 && bIsFourOnLine)	// special case, remove black
					{	rgb_pal.iColors[3] = rgb_pal.iColors[0]; assert(rgb_pal.iColors[0]==rgb_pal.iColors[1]); assert(rgb_pal.iColors[1]==rgb_pal.iColors[2]);	}

					if(bDoFlip)
					{
						for(int i=0; i<(iPalEntries/2); i++)	// swap entries to match group ordering
						{ const int iTmp = rgb_pal.iColors[i]; rgb_pal.iColors[i] = rgb_pal.iColors[iPalEntries-i-1]; rgb_pal.iColors[iPalEntries-i-1] = iTmp; }
					}

					const int offs[] = {iGroup0_index, iGroup1_index, iGroup2_index, iGroup3_index, iGroup4_index};


					int iSumSqLen = 0;
					for(int index=0; index<iPalEntries; index++)
					{
						const int iR = rgb_pal.iPalMul*((rgb_pal.iColors[index]>>0)&0x3ff);
						const int iG = rgb_pal.iPalMul*((rgb_pal.iColors[index]>>10)&0x3ff);
						const int iB = rgb_pal.iPalMul*((rgb_pal.iColors[index]>>20)&0x3ff);

						const int iLeft = offs[index];
						const int iRight = offs[index+1];

						const int iNrElems = iRight-iLeft;
						const SSumsTable &table = sSumsTable[iLeft][iRight];

						const int iSumSq = table.iSqSum;
						const int iTwoSumR = table.iTwoSumR;
						const int iTwoSumG = table.iTwoSumG;
						const int iTwoSumB = table.iTwoSumB;

						// use SUMi_0_to_n (R[i] - iCenR)^2 = ( SUMi_0_to_n R[i]^2 ) + n*iCenR^2 -2*iCenR*( SUMi_0_to_n R[i])
						// and of course use the same property for green and blue.
						const int iLenSQGroup = (TWO_X_THREE*TWO_X_THREE)*iSumSq + iNrElems*(iR*iR+iG*iG+iB*iB) - TWO_X_THREE*(iR*iTwoSumR + iG*iTwoSumG + iB*iTwoSumB);
						iSumSqLen += iLenSQGroup;



						/*for(int i=offs[index]; i<offs[index+1]; i++)
						{
						int iR_in = (uLocalCols[i]>>0)&0xff;
						int iG_in = (uLocalCols[i]>>8)&0xff;
						int iB_in = (uLocalCols[i]>>16)&0xff;

						const int iDr = iR_in - iR;
						const int iDg = iG_in - iG;
						const int iDb = iB_in - iB;
						const int iLenSQ = iDr*iDr + iDg*iDg + iDb*iDb;
						assert( IsValidAdd32bit(iSumSqLen, iLenSQ) );
						iSumSqLen += iLenSQ;
						}*/
					}

					if(iSumSqLen<iSumOfSquaresError)
					{
						iSumOfSquaresError = iSumSqLen;
						*pC0 = color0;
						*pC1 = color1;
					}
				}
			}
		}

	// correct special case, force pC0[0] > pC1[0]
	// this only makes a difference when using alpha
	CorrectColors(pC0, pC1, bIsFourOnLine);

	// return result
	return iSumOfSquaresError;
}

// returns best index
const int GetBestFit(int * piBestDist, const uint32 col_in, const SDXT_RGB_Pal &rgb_pal);

#define DXT1_NOALPHA	0
#define DXT1_ALPHA		1
#define DXT3_ALPHA		2

int CurMode(const int mode=-1)
{
	static int g_iMode = 0;
	if(mode>=0) g_iMode = mode;
	return g_iMode;
}

void ConvertBlockHeuristic3(const uint32 block_in[], const int mode,  const int refinement_value_in, int32 * block_out, uint32 uResultingPixels[])
{
	CurMode(mode);
	const bool bApplyFastMethod = refinement_value_in <= 1;

	int iC0, iC1;
	int sum_dist = -1;
	int * pBestC0 = &iC0; int * pBestC1 = &iC1; int * piCurDist = &sum_dist;

	// sort out unbelievably annoying special cases here
	uint32 block[16];
	memset(block, 0, sizeof(block));
	int iNumActualBlockColors = 0;
	bool bMustForceAlphaDXT1 = false;	// force case color0 <= color1
	for(int k=0; k<16; k++)
	{
		if( CurMode()==DXT1_ALPHA && ((block_in[k]&0xff000000)==0) )
			bMustForceAlphaDXT1 = true;

		if( CurMode()!=DXT1_ALPHA || ((block_in[k]&0xff000000)!=0) )
			block[iNumActualBlockColors++] = (block_in[k]&0x00ffffff);
	}

	if(iNumActualBlockColors > 1)
		QuickSort(block, 0, iNumActualBlockColors-1);

	int bestC0 = -1;
	int bestC1 = -1;

	// insert input colors into uRemainingCols
	// sorted by distance from zero (or by luminosity)
	uint32 uRemainingCols[16];
	int32 iSquaredLengths[16];
	for(int k=0; k<iNumActualBlockColors; k++)
	{
		uRemainingCols[k] = block[k];
		int32 iR = (block[k]>>0)&0xff;
		int32 iG = (block[k]>>8)&0xff;
		int32 iB = (block[k]>>16)&0xff;
		iSquaredLengths[k] = iR*iR+iG*iG+iB*iB;
		assert((iSquaredLengths[k]&0x80000000)==0);
	}

	if(iNumActualBlockColors > 1)
		QuickSort((uint32 *) iSquaredLengths, uRemainingCols, 0, iNumActualBlockColors-1);

	int32 iSumErrorZeroEntries = 0;
	int iCurBestErrorSumSq = -1;
	for(int c=0; c<=iNumActualBlockColors; c++)
	{
		const int iNumRemainingCols = iNumActualBlockColors - c;

		// do the "4 colors on a line" method the first time
		if(c==0 && (!bMustForceAlphaDXT1)) iCurBestErrorSumSq = ProceedMethod3(uRemainingCols, iNumRemainingCols, &bestC0, &bestC1, true, bApplyFastMethod);

		if( (CurMode()==DXT1_NOALPHA) || (CurMode()==DXT1_ALPHA && c==0) )
		{
			// now in reversed order, "3 colors on a line and 4th color fixed to zero"
			int C0, C1;		// skip the first c entries
			const int iCurErrorSumSq = ProceedMethod3(uRemainingCols+c, iNumRemainingCols, &C0, &C1, false, bApplyFastMethod);

			assert( (c==0) || IsValidAdd32bit(iSumErrorZeroEntries, iSquaredLengths[c-1]) );
			if(c>0) iSumErrorZeroEntries += iSquaredLengths[c-1];


			if( (iCurBestErrorSumSq==-1) || ((iCurErrorSumSq+(TWO_X_THREE*TWO_X_THREE)*iSumErrorZeroEntries) < iCurBestErrorSumSq) )
			{
				iCurBestErrorSumSq = iCurErrorSumSq+(TWO_X_THREE*TWO_X_THREE)*iSumErrorZeroEntries;
				bestC0 = C0; bestC1 = C1;
			}
		}
	}


	// refine the result
	const int iDim = refinement_value_in==0 ? 0 : 2;
	const int R0 = (bestC0>>11)&0x1f;	const int R1 = (bestC1>>11)&0x1f;
	const int G0 = (bestC0>>5)&0x3f;	const int G1 = (bestC1>>5)&0x3f;
	const int B0 = (bestC0>>0)&0x1f;	const int B1 = (bestC1>>0)&0x1f;
	const int r0_min = Max(R0-iDim, 0); const int r0_max = Min(R0+iDim, 0x1f);
	const int g0_min = Max(G0-iDim, 0); const int g0_max = Min(G0+iDim, 0x3f);
	const int b0_min = Max(B0-iDim, 0); const int b0_max = Min(B0+iDim, 0x1f);
	const int r1_min = Max(R1-iDim, 0); const int r1_max = Min(R1+iDim, 0x1f);
	const int g1_min = Max(G1-iDim, 0); const int g1_max = Min(G1+iDim, 0x3f);
	const int b1_min = Max(B1-iDim, 0); const int b1_max = Min(B1+iDim, 0x1f);

	int iCurDist = -1;
	if(refinement_value_in==0 || refinement_value_in==3)
	{
		for(int r0=r0_min; r0<=r0_max; r0++)
			for(int g0=g0_min; g0<=g0_max; g0++)
				for(int b0=b0_min; b0<=b0_max; b0++)
					for(int r1=r1_min; r1<=r1_max; r1++)
						for(int g1=g1_min; g1<=g1_max; g1++)
							for(int b1=b1_min; b1<=b1_max; b1++)
							{
								int C0 = ((r0&0x1f)<<11) | ((g0&0x3f)<<5) | ((b0&0x1f)<<0);
								int C1 = ((r1&0x1f)<<11) | ((g1&0x3f)<<5) | ((b1&0x1f)<<0);
								if( (bMustForceAlphaDXT1 && (C0>C1)) ||
									(CurMode()==DXT3_ALPHA && (C0<=C1)) )
								{ const int iTmp = C0; C0 = C1; C1 = iTmp; }	// swap

								SDXT_RGB_Pal rgb_pal;
								GetColors(&rgb_pal, C0, C1);
								CompareResult(block, iNumActualBlockColors, rgb_pal, C0, C1, &iCurDist, &bestC0, &bestC1);	
							}
	}
	else
	{
		const int cB = ((R1&0x1f)<<11) | ((G1&0x3f)<<5) | ((B1&0x1f)<<0);
		for(int r0=r0_min; r0<=r0_max; r0++)
			for(int g0=g0_min; g0<=g0_max; g0++)
				for(int b0=b0_min; b0<=b0_max; b0++)
				{
					int C0 = ((r0&0x1f)<<11) | ((g0&0x3f)<<5) | ((b0&0x1f)<<0);
					int C1 = cB;
					if( (bMustForceAlphaDXT1 && (C0>C1)) ||
						(CurMode()==DXT3_ALPHA && (C0<=C1)) )
					{ const int iTmp = C0; C0 = C1; C1 = iTmp; }	// swap

					SDXT_RGB_Pal rgb_pal;
					GetColors(&rgb_pal, C0, C1);
					CompareResult(block, iNumActualBlockColors, rgb_pal, C0, C1, &iCurDist, &bestC0, &bestC1);	
				}

		const int cA = ((R0&0x1f)<<11) | ((G0&0x3f)<<5) | ((B0&0x1f)<<0);
		for(int r1=r1_min; r1<=r1_max; r1++)
			for(int g1=g1_min; g1<=g1_max; g1++)
				for(int b1=b1_min; b1<=b1_max; b1++)
				{
					int C0 = cA;
					int C1 = ((r1&0x1f)<<11) | ((g1&0x3f)<<5) | ((b1&0x1f)<<0);
					if( (bMustForceAlphaDXT1 && (C0>C1)) ||
						(CurMode()==DXT3_ALPHA && (C0<=C1)) )
					{ const int iTmp = C0; C0 = C1; C1 = iTmp; }	// swap

					SDXT_RGB_Pal rgb_pal;
					GetColors(&rgb_pal, C0, C1);
					CompareResult(block, iNumActualBlockColors, rgb_pal, C0, C1, &iCurDist, &bestC0, &bestC1);	
				}
	}
	


	// overwrite result if better.
	if(((*piCurDist)==-1) || ((*piCurDist)>iCurDist))
	{
		*pBestC0 = bestC0;
		*pBestC1 = bestC1;
		*piCurDist = iCurDist;
	}

	//
	/*const int iBestOrgC0 = *pBestC0;
	const int iBestOrgC1 = *pBestC1;
	const int iBestOrgDist = *piCurDist;

	*pBestC0 = -1;
	*pBestC1 = -1;
	*piCurDist = -1;
	for(int C1=0; C1<0x10000; C1++)
		for(int C0=0; C0<0x10000; C0++)
		{
			int iColors[4];
			GetColors(iColors, C0, C1);
			CompareResult(block, iNumActualBlockColors, iColors, C0, C1, piCurDist, pBestC0, pBestC1);	
		}

	static int iCounter = 0;
	printf("block_nr: %d, org_dist: %d, new_dist: %d\n", iCounter, iBestOrgDist, piCurDist[0]);
	printf("org_col_pair: (%d, %d, %d), (%d, %d, %d)\n", (iBestOrgC0>>11)&0x1f, (iBestOrgC0>>5)&0x3f, (iBestOrgC0>>0)&0x1f, (iBestOrgC1>>11)&0x1f, (iBestOrgC1>>5)&0x3f, (iBestOrgC1>>0)&0x1f);
	printf("new_col_pair: (%d, %d, %d), (%d, %d, %d)\n", (pBestC0[0]>>11)&0x1f, (pBestC0[0]>>5)&0x3f, (pBestC0[0]>>0)&0x1f, (pBestC1[0]>>11)&0x1f, (pBestC1[0]>>5)&0x3f, (pBestC1[0]>>0)&0x1f);
	++iCounter;*/
	//



	// output dxt block
	// block_out
	// uResultingPixels
	SDXT_RGB_Pal rgb_pal;
	GetColors(&rgb_pal, pBestC0[0], pBestC1[0]);

	// output dxt codes in a little/big endian friendly way.
	((uint8 *) block_out)[0] = (uint8) ((pBestC0[0]>>0)&0xff);	// lo c0
	((uint8 *) block_out)[1] = (uint8) ((pBestC0[0]>>8)&0xff);	// hi c0
	((uint8 *) block_out)[2] = (uint8) ((pBestC1[0]>>0)&0xff);	// lo c1
	((uint8 *) block_out)[3] = (uint8) ((pBestC1[0]>>8)&0xff);	// hi c1

	for(int i=0; i<4; i++)
	{
		uint8 uIndices = 0;
		for(int j=0; j<4; j++)
		{
			int iBestDist;
			
			int cur_index; // index == 3, is reserved for alpha in dxt1_alpha mode
			if( CurMode()==DXT1_ALPHA && ((block_in[i*4+j]&0xff000000)==0) )
			{
				cur_index = 3;
				uResultingPixels[i*4+j] = 0;	// rgb forced to black
			}
			else
			{
				cur_index = GetBestFit(&iBestDist, block_in[i*4+j], rgb_pal);

				int r0 = ( rgb_pal.iPalMul * ((rgb_pal.iColors[cur_index]>>0)&0x3ff) ) / 3;
				int g0 = ( rgb_pal.iPalMul * ((rgb_pal.iColors[cur_index]>>10)&0x3ff) ) / 3;
				int b0 = ( rgb_pal.iPalMul * ((rgb_pal.iColors[cur_index]>>20)&0x3ff) ) / 3;
				r0 = Min((r0>>1)+(r0&0x1), 0xff);	// round to nearest
				g0 = Min((g0>>1)+(g0&0x1), 0xff);
				b0 = Min((b0>>1)+(b0&0x1), 0xff);
				const int color888 = (b0<<16) | (g0<<8) | (r0<<0);

				uResultingPixels[i*4+j] = (color888&0x00ffffff) | (block_in[i*4+j]&0xff000000);

				assert( (!bMustForceAlphaDXT1) || (pBestC0[0]<=pBestC1[0]) );
				assert( (cur_index<3) || (CurMode()!=DXT1_ALPHA) || (pBestC0[0]>pBestC1[0]) );
			}
			assert(cur_index>=0 && cur_index<4);

			if((pBestC0[0]>pBestC1[0]) || (CurMode()!=DXT1_NOALPHA && CurMode()!=DXT1_ALPHA))
			{
				if(cur_index == 3) cur_index = 1;
				else if(cur_index == 1) cur_index = 2;
				else if(cur_index == 2) cur_index = 3;
			}
			else
			{
				if(cur_index == 1) cur_index = 2;
				else if(cur_index == 2) cur_index = 1;
			}

			/*
			// 00 = color_0, 01 = color_1, 10 = color_2, 11 = color_3
			// These 2-bit codes correspond to the 2-bit fields 
			// stored in the 64-bit block.
			color_2 = (2 * color_0 + color_1 + 1) / 3;
			color_3 = (color_0 + 2 * color_1 + 1) / 3;

			// Three-color block: derive the other color.
			// 00 = color_0,  01 = color_1,  10 = color_2,  
			// 11 = transparent.
			// These 2-bit codes correspond to the 2-bit fields 
			// stored in the 64-bit block. 
			color_2 = (color_0 + color_1) / 2;    
			color_3 = transparent;    
			*/


			uIndices |= (cur_index<<(2*j));
		}

		((uint8 *) block_out)[4+i] = uIndices;
	}
}




//////////// utils /////////////

const bool IsVecValidBasisVector(const Vec3 &v)
{
	const unsigned int * uVec = (const unsigned int *) &v;
	// if vector is valid
	if(		(uVec[0]!=0xffc00000 && uVec[1]!=0xffc00000 && uVec[2]!=0xffc00000) && (Length(v)>0.0001f) 	)
		return true;
	else
		return false;
}

int solveCubic(const double eqn[], double res[]);

bool CompareRows(const double row1[], const double row2[], const int iEntries)
{
	bool bMatchNotDone = true;
	int i=0;
	while(i<iEntries && bMatchNotDone)
		if(	fabs(row1[i]) != fabs(row2[i]) )
			bMatchNotDone = false;
		else
			++i;

	return (i==iEntries) || (fabs(row1[i]) > fabs(row2[i]));
}

void FindBestTrans(const Vec3 vVerts[], const int nr_verts, Mat33 * pRot, Vec3 * pPos)
{
	double mat_data[9];
	int k;
	for(k=0; k<9; k++ )
		mat_data[k] = 0;

	double cx = 0, cy = 0, cz = 0;
	for(k=0; k<nr_verts; k++)
	{
		cx += vVerts[k].x;
		cy += vVerts[k].y;
		cz += vVerts[k].z;
	}
	cx /= nr_verts; cy /= nr_verts; cz /= nr_verts;
	const Vec3 vCen = Vec3((const float) cx, (const float) cy, (const float) cz);
	*pPos = vCen;
	for(k=0; k<nr_verts; k++)
	{
		const Vec3 vV = vVerts[k] - vCen;
		const double sxx = vV.x*vV.x;
		const double syy = vV.y*vV.y;
		const double szz = vV.z*vV.z;
		const double sxy = vV.x*vV.y;
		const double sxz = vV.x*vV.z;
		const double syz = vV.y*vV.z;

		mat_data[0] += sxx; mat_data[1] += sxy; mat_data[2] += sxz;
		mat_data[3] += sxy; mat_data[4] += syy; mat_data[5] += syz;
		mat_data[6] += sxz; mat_data[7] += syz; mat_data[8] += szz;
	}

	// apply a uniform scale to keep future
	// calculations within a reasonable range.
	double SUM = 0;
	for(k=0; k<9; k++)
		SUM += mat_data[k];
	for(k=0; k<9; k++)
		mat_data[k] /= (SUM/16);

	// find eigenvalues by solving
	// X^3 + A*X^2 + B*X + C

	// negated trace
	const double A = -(mat_data[0]+mat_data[4]+mat_data[8]);

	// not sure if there's anything famous about this
	const double B = -(mat_data[3]*mat_data[1] + mat_data[6]*mat_data[2] + mat_data[7]*mat_data[5]
	-mat_data[0]*mat_data[4] -mat_data[8]*mat_data[0] -mat_data[8]*mat_data[4]  );

	// negated determinant
	const double C = -( mat_data[3]*mat_data[7]*mat_data[2] + mat_data[6]*mat_data[1]*mat_data[5] + 
		mat_data[8]*mat_data[0]*mat_data[4] - mat_data[0]*mat_data[7]*mat_data[5]
		- mat_data[8]*mat_data[3]*mat_data[1] - mat_data[4]*mat_data[6]*mat_data[2] );




	const double coeffs[] = {C, B, A, 1};
	//const float x0 = 3, x1 = 3, x2 = 5;
	//const double coeffs[] = {-x0*x1*x2, x2*x0+x2*x1+x0*x1, -(x0+x1+x2), 1};

	// find eigenvalues
	double res[3];
	int iNrRoots = solveCubic( coeffs, res);
	//if(C==0) { iNrRoots=1; res[0]=0; }


	// sort eigenvalues
	assert(iNrRoots>0 && iNrRoots<4);
	if(iNrRoots > 1)
	{
		if(fabs(res[0])<fabs(res[1])) { const double fTmp = res[0]; res[0] = res[1]; res[1] = fTmp; }
		if(iNrRoots > 2)
		{
			if(fabs(res[1])<fabs(res[2])) { const double fTmp = res[1]; res[1] = res[2]; res[2] = fTmp; }
			if(fabs(res[0])<fabs(res[1])) { const double fTmp = res[0]; res[0] = res[1]; res[1] = fTmp; }
		}
	}

	for(int i=1; i<iNrRoots; i++)
		assert(fabs(res[i-1])>=fabs(res[i]));


	int iNumEigenVectors = 0;
	Vec3 vResultingBasis[3];

	// Apply Gaussian Elimination
	for(k=0; k<iNrRoots; k++)
	{
		// make a copy of the input matrix.
		double mat_work[9];
		for(int i=0; i<9; i++)
			mat_work[i] = mat_data[i];

		// subtract the current eigenvalue
		// from the diagonal of the copy
		mat_work[0*3+0] = mat_data[0*3+0] - res[k];
		mat_work[1*3+1] = mat_data[1*3+1] - res[k];
		mat_work[2*3+2] = mat_data[2*3+2] - res[k];

		// place row with numerically largest first component at the top.
		const bool bThirdRowBetter = CompareRows(mat_work+2*3, mat_work+1*3, 3);
		int iReplacementRow = bThirdRowBetter ? 2 : 1;
		if( CompareRows(mat_work+iReplacementRow*3, mat_work+0*3, 3) )		// must swap rows
			for(int q=0; q<3; q++)
			{ const double dTmp = mat_work[0*3+q]; mat_work[0*3+q] = mat_work[iReplacementRow*3+q]; mat_work[iReplacementRow*3+q] = dTmp; }

		// scale first row by non zero component
		bool bFoundComponent = false;
		int c=0;
		while(c<3 && (!bFoundComponent))
			if(mat_work[0*3+c]!=0) bFoundComponent = true;
			else ++c;
		const double dDivValue = bFoundComponent ? mat_work[0*3+c] : 0;
		if(bFoundComponent)
		{
			mat_work[0*3+2] /= dDivValue; mat_work[0*3+1] /= dDivValue; mat_work[0*3+0] /= dDivValue; mat_work[0*3+c] = 1;

			// subtract first row by appropriate multiplier from row two and three.
			for(int r=1; r<3; r++)
				if(mat_work[r*3+c]!=0)
				{
					for(int c1=(c+1); c1<3; c1++)
					{	mat_work[r*3+c1] -= (mat_work[0*3+c1] * mat_work[r*3+c]); }
						
					mat_work[r*3+c] = 0;
				}
		}

		// now focus on the 2x2 submatrix we get by removing first row and column.
		// As second row, choose between second and third row the one which has the
		// numerically largest second component.
		if( CompareRows(mat_work+(2*3+1), mat_work+(1*3+1), 2) ) // swap rows if necessary
			for(int q=1; q<3; q++)
			{ const double dTmp = mat_work[1*3+q]; mat_work[1*3+q] = mat_work[2*3+q]; mat_work[2*3+q] = dTmp; }

		// scale second row by diagonal component
		if(mat_work[1*3+1]!=0)
		{
			// scale
			mat_work[1*3+2] /= mat_work[1*3+1]; mat_work[1*3+1]=1;

			// subtract from third row, second row by an appropriate multiplier
			mat_work[2*3+2] -= (mat_work[1*3+2] * mat_work[2*3+1]);
			mat_work[2*3+1] = 0;
		}
		else if(mat_work[1*3+2]!=0)
		{
			mat_work[1*3+2]=1; mat_work[2*3+2] = 0;
		}

		// Cross first and second row to obtain a solution
		double nx = mat_work[0*3+1]*mat_work[1*3+2] - mat_work[0*3+2]*mat_work[1*3+1];
		double ny = mat_work[0*3+2]*mat_work[1*3+0] - mat_work[0*3+0]*mat_work[1*3+2];
		double nz = mat_work[0*3+0]*mat_work[1*3+1] - mat_work[0*3+1]*mat_work[1*3+0];

		const double dLen = sqrt(nx*nx+ny*ny+nz*nz);
		if(dLen>0)
		{
			vResultingBasis[iNumEigenVectors] = Vec3((float) (nx/dLen), (float) (ny/dLen), (float) (nz/dLen));

			// Apply the Gram-Schmidt Orthonormalization Process
			const int iCurI = iNumEigenVectors;
			for(int i=0; i<iCurI; i++)
				vResultingBasis[iCurI] = vResultingBasis[iCurI] - vResultingBasis[i]*(vResultingBasis[i]*vResultingBasis[iCurI]);

			if(IsVecValidBasisVector(vResultingBasis[iCurI]))
			{
				vResultingBasis[iCurI] = Normalize( vResultingBasis[iCurI] );
				++iNumEigenVectors;
			}
		}
	}


	// extract primary eigenvector
	Vec3 vXAxis, vYAxis, vZAxis;
	if(iNumEigenVectors>0) vXAxis = vResultingBasis[0];
	else vXAxis = Vec3(1,0,0);

	// validate second eigenvector
	if(iNumEigenVectors>1)
		vYAxis = vResultingBasis[1];
	else
	{
		if( fabs(vXAxis.x)<fabs(vXAxis.y) && fabs(vXAxis.x)<fabs(vXAxis.z)  )
			vYAxis = Vec3(1,0,0);
		else if( fabs(vXAxis.y)<fabs(vXAxis.z) )
			vYAxis = Vec3(0,1,0);
		else
			vYAxis = Vec3(0,0,1);
	}

	vYAxis = Normalize( vYAxis - vXAxis*(vXAxis*vYAxis) );
	vZAxis = Cross(vXAxis, vYAxis);		


	SetColumn(pRot, 0, vXAxis);
	SetColumn(pRot, 1, vYAxis);
	SetColumn(pRot, 2, vZAxis);

	assert(Determinant(*pRot)>0);
}

#ifndef M_PI
#define M_PI	3.1415926535897932384626433832795
#endif

int solveCubic(const double eqn[], double res[])
{
	double a, b, c, q, r, Q, R;
	double c3, Q3, R2, CR2, CQ3;

	// If the cubic coefficient is zero, we have a quadratic equation.
	c3 = eqn[3];
	if (c3 == 0)
	{
		assert(false);
		//return QuadCurve2D.solveQuadratic(eqn, res);
	}

	a = eqn[2] / c3;
	b = eqn[1] / c3;
	c = eqn[0] / c3;

	// We now need to solve x^3 + ax^2 + bx + c = 0.
	q = a * a - 3 * b;
	r = 2 * a * a * a - 9 * a * b + 27 * c;

	Q = q / 9;
	R = r / 54;

	Q3 = Q * Q * Q;
	R2 = R * R;

	CR2 = 729 * r * r;
	CQ3 = 2916 * q * q * q;

	if (R == 0 && Q == 0)
	{
		// The GNU Scientific Library would return three identical
		// solutions in this case.
		res[0] = -a/3;
		return 1;
	}

	if (CR2 == CQ3) 
	{
		/* this test is actually R2 == Q3, written in a form suitable
		for exact computation with integers */

		/* Due to finite precision some double roots may be missed, and
		considered to be a pair of complex roots z = x +/- epsilon i
		close to the real axis. */

		double sqrtQ = sqrt(Q);

		if (R > 0)
		{
			res[0] = -2 * sqrtQ - a/3;
			res[1] = sqrtQ - a/3;
		}
		else
		{
			res[0] = -sqrtQ - a/3;
			res[1] = 2 * sqrtQ - a/3;
		}
		return 2;
	}

	if (CR2 < CQ3) /* equivalent to R2 < Q3 */
	{
		double sqrtQ = sqrt(Q);
		double sqrtQ3 = sqrtQ * sqrtQ * sqrtQ;
		double myfrac = R / sqrtQ3;
		double clamped_frac = (myfrac<(-1)) ? (-1) : ( (myfrac>1) ? 1 : myfrac );
		double theta = acos(clamped_frac);
		double norm = -2 * sqrtQ;
		res[0] = norm * cos(theta / 3) - a / 3;
		res[1] = norm * cos((theta + 2.0 * M_PI) / 3) - a/3;
		res[2] = norm * cos((theta - 2.0 * M_PI) / 3) - a/3;

		// The GNU Scientific Library sorts the results. We don't.
		return 3;
	}

	double sgnR = (R >= 0 ? 1 : -1);
	double A = -sgnR * pow( (R<0?(-R):R) + sqrt(R2 - Q3), 1.0/3.0);
	double B = Q / A ;
	res[0] = A + B - a/3;
	return 1;
}

int ColorAs888(const Vec3 &v0)
{
	const int iR0 = Min(v0.x<0?0:((int) (v0.x+0.5)), 0xff);
	const int iG0 = Min(v0.y<0?0:((int) (v0.y+0.5)), 0xff);
	const int iB0 = Min(v0.z<0?0:((int) (v0.z+0.5)), 0xff);
	return (iR0<<0) | (iG0<<8) | (iB0<<16);
}

void ComputeColorsAs565(const Vec3 &v0, const Vec3 &v1, int * pColor0, int * pColor1)
{
	// snap C1 and C2 to closest 5:6:5 color
	const int iColor0 = ColorAs888(v0);
	const int iColor1 = ColorAs888(v1);

	*pColor0 = GetC565_from_C888(iColor0);
	*pColor1 = GetC565_from_C888(iColor1);
}

void ComputeColorAs565(const Vec3 &v0, int * pColor0)
{
	const int iColor0 = ColorAs888(v0);
	*pColor0 = GetC565_from_C888(iColor0);
}

/*int GetC565_from_C888(const int color)
{
	return (((color>>19)&0x1f)<<0) | (((color>>10)&0x3f)<<5) | (((color>>3)&0x1f)<<11);
}*/
int GetC565_from_C888(const int color)
{
	const int iRed = d3d_spec_8_to_5( ((color>>0)&0xff) );
	const int iGreen = d3d_spec_8_to_6( ((color>>8)&0xff) );
	const int iBlue = d3d_spec_8_to_5( ((color>>16)&0xff) );
	return (iRed<<11) | (iGreen<<5) | (iBlue<<0);
	//return (((color>>19)&0x1f)<<11) | (((color>>10)&0x3f)<<5) | (((color>>3)&0x1f)<<0);
}
int GetC888_from_C565(const int color)
{
	const int iRed = d3d_spec_5_to_8( ((color>>11)&0x1f) );
	const int iGreen = d3d_spec_6_to_8( ((color>>5)&0x3f) );
	const int iBlue = d3d_spec_5_to_8( ((color>>0)&0x1f) );
	return (iRed<<0) | (iGreen<<8) | (iBlue<<16);
}

void GetColors(SDXT_RGB_Pal * pRGB_Pal, const int32 c0, const int32 c1)
{
	const int color0_888 = GetC888_from_C565(c0);
	const int color1_888 = GetC888_from_C565(c1);

	const int r0 = (color0_888>>0)&0xff;
	const int g0 = (color0_888>>8)&0xff;
	const int b0 = (color0_888>>16)&0xff;
	const int r1 = (color1_888>>0)&0xff;
	const int g1 = (color1_888>>8)&0xff;
	const int b1 = (color1_888>>16)&0xff;

	uint32 cA=0, cB=0;
	if(c0 > c1)
	{
		pRGB_Pal->iPalMul = (TWO_X_THREE/3);
		
		cA |= ( (2*r0+r1) << 0 );
		cA |= ( (2*g0+g1) << 10 );
		cA |= ( (2*b0+b1) << 20 );
		cB |= ( (2*r1+r0) << 0 );
		cB |= ( (2*g1+g0) << 10 );
		cB |= ( (2*b1+b0) << 20 );
	}
	else
	{
		pRGB_Pal->iPalMul = (TWO_X_THREE/2);

		cA |= ( (r0+r1) << 0 );
		cA |= ( (g0+g1) << 10 );
		cA |= ( (b0+b1) << 20 );
	}

	const int iMult = pRGB_Pal->iPalMul ^ 1;	// 2 becomes 3 and vice versa
	assert( (iMult==2&&pRGB_Pal->iPalMul==3) || (iMult==3&&pRGB_Pal->iPalMul==2) );
	pRGB_Pal->iColors[0] = ((iMult*r0)<<0) | ((iMult*g0)<<10) | ((iMult*b0)<<20);
	pRGB_Pal->iColors[1] = cA;
	pRGB_Pal->iColors[2] = cB;
	pRGB_Pal->iColors[3] = ((iMult*r1)<<0) | ((iMult*g1)<<10) | ((iMult*b1)<<20);

	if(c0<=c1)
	{
		assert(CurMode()==DXT1_NOALPHA || CurMode()==DXT1_ALPHA || (CurMode()==DXT3_ALPHA && c0==c1));
		pRGB_Pal->iColors[2]=pRGB_Pal->iColors[3];
		if( CurMode()==DXT1_NOALPHA ) pRGB_Pal->iColors[3]=0;		// we get to use black.
	}
}

void CompareResult(const uint32 block[], const int iBlockColors, const SDXT_RGB_Pal &rgb_pal, const int32 color0, const int32 color1, int * piCurDist, int * piBestC0, int * piBestC1)
{
	int iSumSQ = 0;
	int k=0;
	bool bIsWorse = false;
	while(!bIsWorse && k<iBlockColors)
	{
		// iColNr is 0,1,2 or 3
		int iBest_distSq;
		const int iBestColNr = GetBestFit(&iBest_distSq, block[k], rgb_pal);
		const int iRecK = k;
		while((k<(iBlockColors-1)) && ((block[k]&0x00ffffff) == (block[k+1]&0x00ffffff)))
			++k;
		const int iNumInstances = (k-iRecK)+1;
		const int iTotalError = iNumInstances*iBest_distSq;

		iSumSQ += iTotalError;

		// early out.
		if((*piCurDist)!=(-1) && iSumSQ>(*piCurDist))
			bIsWorse = true;

		// next
		if(!bIsWorse) ++k;
	}

	if(!bIsWorse)
	{
		assert(iSumSQ>=0);
		*piCurDist = iSumSQ;
		*piBestC0 = color0;
		*piBestC1 = color1;
	}
}

const int GetBestFit(int * piBestDist, const uint32 col_in, const SDXT_RGB_Pal &rgb_pal)
{
	int iBestIndex = -1;
	int iBestDist;

	const int r0 = TWO_X_THREE*((col_in>>0)&0xff);
	const int g0 = TWO_X_THREE*((col_in>>8)&0xff);
	const int b0 = TWO_X_THREE*((col_in>>16)&0xff);

	for(int i=0; i<4; i++)
	{
		const int dr = rgb_pal.iPalMul*((rgb_pal.iColors[i]>>0)&0x3ff)-r0;
		const int dg = rgb_pal.iPalMul*((rgb_pal.iColors[i]>>10)&0x3ff)-g0;
		const int db = rgb_pal.iPalMul*((rgb_pal.iColors[i]>>20)&0x3ff)-b0;
		const int iCurDist = dr*dr+dg*dg+db*db;
		if(i==0 || iCurDist<iBestDist)
		{
			iBestDist = iCurDist;
			iBestIndex = i;
		}
	}

	*piBestDist = iBestDist;
	return iBestIndex;
}


// a secondary less accurate approach.
int Process3Fast(const uint32 uLocalCols[], const Vec3 vSortedVerts[], const int iNumCols, const bool bIsFourOnLine, const Vec3 &vFinPos, const Mat33 &mFinRot, int * pC0, int * pC1)
{
	// derive an initial worst case error from center of colors.
	int iCenX = 0, iCenY = 0, iCenZ = 0;
	for(int i=0; i<iNumCols; i++)
	{ iCenX += ((uLocalCols[i]>>0)&0xff); iCenY += ((uLocalCols[i]>>8)&0xff); iCenZ += ((uLocalCols[i]>>16)&0xff); }

	iCenX /= iNumCols; iCenY /= iNumCols; iCenZ /= iNumCols;

	int single_color = GetC565_from_C888( (iCenX<<0) | (iCenY<<8) | (iCenZ<<16) );
	iCenX = d3d_spec_5_to_8( d3d_spec_8_to_5( iCenX ) );	// loss by 888 to 565
	iCenY = d3d_spec_6_to_8( d3d_spec_8_to_6( iCenY ) );
	iCenZ = d3d_spec_5_to_8( d3d_spec_8_to_5( iCenZ ) );

	int iSumOfSquaresError = 0;
	for(int c=0; c<iNumCols; c++)
	{
		const int iDr = ((uLocalCols[c]>>0)&0xff) - iCenX;
		const int iDg = ((uLocalCols[c]>>8)&0xff) - iCenY;
		const int iDb = ((uLocalCols[c]>>16)&0xff) - iCenZ;
		const int iLenSQ = iDr*iDr + iDg*iDg + iDb*iDb;
		iSumOfSquaresError += (TWO_X_THREE*TWO_X_THREE)*iLenSQ;
	}

	*pC0 = single_color;
	*pC1 = single_color;


	const float fMinX = vSortedVerts[0].x;
	const float fMaxX = vSortedVerts[iNumCols-1].x;

	for(int c=0; c<2; c++)
	{
		const float fDelta = (fMaxX-fMinX)*0.1f*(c-1);	// 10% of delta

		float fXlocations[4];
		fXlocations[0] = fMinX-fDelta;
		fXlocations[3] = fMaxX+fDelta;

		// determine groups
		int iNrGroups = 4;
		if(bIsFourOnLine)
		{
			fXlocations[1] = (fXlocations[0]*2 + fXlocations[3])/3;
			fXlocations[2] = (fXlocations[0] + fXlocations[3]*2)/3;
		}
		else
		{
			iNrGroups = 3;
			fXlocations[2] = fXlocations[3];
			fXlocations[1] = (fXlocations[0] + fXlocations[2])/2;
		}

		// assign colors from vSortedVerts[]
		int iNrInGroup[] = {0,0,0,0};	// elements assigned to group
		int iGroupIndices[4][16];
		for(int i=0; i<iNumCols; i++)
		{
			float fDists[4];
			for(int j=0; j<iNrGroups; j++)
				fDists[j] = fabs(fXlocations[j]-vSortedVerts[i].x);
			int index_best = 0;
			float dist_best = fDists[0];
			for(int j=1; j<iNrGroups; j++)
				if(fDists[j]<dist_best)
				{ index_best=j; dist_best=fDists[j]; }

				iGroupIndices[index_best][iNrInGroup[index_best]++] = i;
		}

		const int m = iNrInGroup[0];
		const int n = iNrInGroup[1];
		const int p = iNrInGroup[2];
		const int q = bIsFourOnLine ? iNrInGroup[3] : 0;

		int A = bIsFourOnLine ? (9*m + 4*n + p) : (4*m+n);
		int B = bIsFourOnLine ? (2*n + 2*p) : n;
		int C = bIsFourOnLine ? (9*q + 4*p + n) : (4*p+n);

		const int iDeterminant = A*C - B*B;

		if(iDeterminant != 0)
		{
			const int iGroup0_index = 0;	// and m elems from then on
			const int iGroup1_index = m;	// and n elems from then on
			const int iGroup2_index = m+n;	// and p elems from then on
			const int iGroup3_index = m+n+p;// and q elems from then on
			const int iGroup4_index = m+n+p+q;	// end of group3 (and not a real group)

			double dSum1X = 0, dSum2X = 0, dSum3X = 0, dSum4X = 0;
			double dSum1Y = 0, dSum2Y = 0, dSum3Y = 0, dSum4Y = 0;
			double dSum1Z = 0, dSum2Z = 0, dSum3Z = 0, dSum4Z = 0;

			for(int g=0; g<iNrGroups; g++)
			{
				for(int j=0; j<iNrInGroup[g]; j++)
				{
					const int index = iGroupIndices[g][j];
					const Vec3 &vV = vSortedVerts[index];
					if(g==0)
					{ dSum1X += vV.x; dSum1Y += vV.y; dSum1Z += vV.z; }
					else if(g==1)
					{ dSum2X += vV.x; dSum2Y += vV.y; dSum2Z += vV.z; }
					else if(g==2)
					{ dSum3X += vV.x; dSum3Y += vV.y; dSum3Z += vV.z; }
					else
					{ assert(g==3); dSum4X += vV.x; dSum4Y += vV.y; dSum4Z += vV.z; }
				}
			}


			const double C1x = bIsFourOnLine ? (9*dSum1X  +  6*dSum2X  +  3*dSum3X) : (4*dSum1X + 2*dSum2X);
			const double C2x = bIsFourOnLine ? (9*dSum4X  +  6*dSum3X  +  3*dSum2X) : (4*dSum3X + 2*dSum2X);

			const double C1y = bIsFourOnLine ? (9*dSum1Y  +  6*dSum2Y  +  3*dSum3Y) : (4*dSum1Y + 2*dSum2Y);
			const double C2y = bIsFourOnLine ? (9*dSum4Y  +  6*dSum3Y  +  3*dSum2Y) : (4*dSum3Y + 2*dSum2Y);

			const double C1z = bIsFourOnLine ? (9*dSum1Z  +  6*dSum2Z  +  3*dSum3Z) : (4*dSum1Z + 2*dSum2Z);
			const double C2z = bIsFourOnLine ? (9*dSum4Z  +  6*dSum3Z  +  3*dSum2Z) : (4*dSum3Z + 2*dSum2Z);

			// 4 points on line
			// (9*m + 4*n + p)*c0.x + (2*n + 2*p)*c1.x  = 9*SUMi_to_m  (  C[i].x )  +  6*SUMj_to_n  (  C[j].x )  +  3*SUMk_to_p (  C[k].x )
			// (2*p+2*n)*c0.x +  (9*q + 4*p + n)*c1.x   = 9*SUMl_to_q  (  C[l].x )  +  6*SUMk_to_p  (  C[k].x )  +  3*SUMj_to_n (  C[j].x )

			// 3 points on line
			// (4*m+n)*c0.x + n*c1.x = 4 * (SUMi_to_m C[i].x) + 2 * (SUMj_to_n C[j].x)
			// n*c0.x + (4*p+n)*c1.x = 4 * (SUMk_to_p C[k].x) + 2 * (SUMj_to_n C[j].x)

			// solve equation
			// | A B | |x|   |C1|
			// | B C | |y| = |C2|

			// |x|			         | C -B | |C1|
			// |y| = (1/(A*C-B*B)) * |-B  A | |C2|

			const float x0 = (float) ( (C*C1x - B*C2x)/iDeterminant );
			const float x1 = (float) ( ((-B)*C1x + A*C2x)/iDeterminant );

			const float y0 = (float) ( (C*C1y - B*C2y)/iDeterminant );
			const float y1 = (float) ( ((-B)*C1y + A*C2y)/iDeterminant );

			const float z0 = (float) ( (C*C1z - B*C2z)/iDeterminant );
			const float z1 = (float) ( ((-B)*C1z + A*C2z)/iDeterminant );



			// evaluate error
			const Vec3 vC0 = vFinPos + mFinRot*Vec3(x0,y0,z0);
			const Vec3 vC1 = vFinPos + mFinRot*Vec3(x1,y1,z1);

			int color0, color1;
			ComputeColorsAs565(vC0, vC1, &color0, &color1);
			bool bDoFlip = false;
			if( ((color0 > color1) ^ bIsFourOnLine) )
			{	const int iTmp = color0; color0 = color1; color1 = iTmp; bDoFlip=true; }

			// special case: color0 == color1 and bIsFourOnLine
			// this will have to be allowed for now to compute a correct error value.
			// the code will correct it before returning from this function.
			assert(((color0 >= color1)&&bIsFourOnLine) || ((color0 <= color1)&&(!bIsFourOnLine)));

			SDXT_RGB_Pal rgb_pal;
			GetColors(&rgb_pal, color0, color1);
			if(color0==color1 && bIsFourOnLine)	// special case, remove black
			{	rgb_pal.iColors[3] = rgb_pal.iColors[0]; assert(rgb_pal.iColors[0]==rgb_pal.iColors[1]); assert(rgb_pal.iColors[1]==rgb_pal.iColors[2]);	}

			if(bDoFlip)
			{
				for(int i=0; i<(iNrGroups/2); i++)	// swap entries to match group ordering
				{ const int iTmp = rgb_pal.iColors[i]; rgb_pal.iColors[i] = rgb_pal.iColors[iNrGroups-i-1]; rgb_pal.iColors[iNrGroups-i-1] = iTmp; }
			}

			int iSumSqLen = 0;
			for(int index=0; index<iNrGroups; index++)
			{
				// multiplied by 6
				const int iR = rgb_pal.iPalMul*((rgb_pal.iColors[index]>>0)&0x3ff);
				const int iG = rgb_pal.iPalMul*((rgb_pal.iColors[index]>>10)&0x3ff);
				const int iB = rgb_pal.iPalMul*((rgb_pal.iColors[index]>>20)&0x3ff);

				for(int i=0; i<iNrInGroup[index]; i++)
				{
					const int iVertIndex = iGroupIndices[index][i];
					const int iR_in = (uLocalCols[iVertIndex]>>0)&0xff;
					const int iG_in = (uLocalCols[iVertIndex]>>8)&0xff;
					const int iB_in = (uLocalCols[iVertIndex]>>16)&0xff;

					const int iDr = TWO_X_THREE*iR_in - iR;
					const int iDg = TWO_X_THREE*iG_in - iG;
					const int iDb = TWO_X_THREE*iB_in - iB;
					const int iLenSQ = iDr*iDr + iDg*iDg + iDb*iDb;
					iSumSqLen += iLenSQ;
				}
			}

			if(iSumSqLen<iSumOfSquaresError)
			{
				iSumOfSquaresError = iSumSqLen;
				*pC0 = color0;
				*pC1 = color1;
			}		
		}
	}

	// correct special case, force pC0[0] > pC1[0]
	// this only makes a difference when using alpha
	CorrectColors(pC0, pC1, bIsFourOnLine);

	// return result
	return iSumOfSquaresError;
}