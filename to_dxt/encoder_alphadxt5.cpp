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

#include "encoder_alphadxt5.h"
#include <assert.h>


void GetBestAlpha(const unsigned char alpha_block_in[], const int ref_mode, int * piBestA0, int * piBestA1, float * pfCurBestDist);
inline void GetAlphas(const float fA0, const float fA1, float fColorsOut[]);
int FindBestMatch(const unsigned char uAlpha, const float fPalette[]);

void ConvertBlockAlphaDXT5(const uint8 block_in[], const int refinement_value_in, int32 * block_out, uint8 uResultingAlphas[])
{
	int iBestA0, iBestA1;
	float fBestSumSq = -1;	// -1 is no result yet
	GetBestAlpha(block_in, refinement_value_in<0?0:(refinement_value_in>3?3:refinement_value_in), &iBestA0, &iBestA1, &fBestSumSq);	

	float fAlphaPalette[8];
	GetAlphas((const float) iBestA0, (const float) iBestA1, fAlphaPalette);
	((uint8 *) block_out)[0] = iBestA0;		// fill up first short
	((uint8 *) block_out)[1] = iBestA1;
	int best_entries[16];
	for(int i=0; i<16; i++)
	{
		const int ibest_entry = FindBestMatch(block_in[i], fAlphaPalette);
		uResultingAlphas[i] = (int) (fAlphaPalette[ibest_entry]+0.5f);

		// change index from sorted palette to d3d format
		int iD3D_Index;
		if(iBestA0 > iBestA1)
			iD3D_Index = (ibest_entry!=0 && ibest_entry!=7) ? (8-ibest_entry) : (ibest_entry==0?1:0);
		else
			iD3D_Index = (ibest_entry!=0 && ibest_entry!=1 && ibest_entry!=6) ? ibest_entry : (ibest_entry==0?6:(ibest_entry==1?0:1));

		/*
		case: A0 <= A1

		0 --> 6		0
		1 --> 0		A0
		2 --> 2
		3 --> 3
		4 --> 4
		5 --> 5
		6 --> 1		A1
		7 --> 7		255 */


		assert(iD3D_Index>=0 && iD3D_Index<8);
		best_entries[i] = iD3D_Index;
	}

	// insert bits correctly in an endian independant way.
	for(int q=0; q<2; q++)
	{
		const int iByteOffs = q*3;
		const int iIOffs = q*8;
		((uint8 *) block_out)[2+iByteOffs] = ((best_entries[2+iIOffs]&0x3)<<6) | (best_entries[1+iIOffs]<<3) | (best_entries[0+iIOffs]<<0);
		((uint8 *) block_out)[3+iByteOffs] = ((best_entries[5+iIOffs]&0x1)<<7) | (best_entries[4+iIOffs]<<4) | (best_entries[3+iIOffs]<<1) | (best_entries[2+iIOffs]>>2);
		((uint8 *) block_out)[4+iByteOffs] = (best_entries[7+iIOffs]<<5) | (best_entries[6+iIOffs]<<2) | (best_entries[5+iIOffs]>>1);
	}
}

// returns index
int FindBestMatch(const unsigned char uAlpha, const float fPalette[])
{
	// 8 entries
	float fDist = fabsf( uAlpha-fPalette[0] );
	int ibest_entry = 0;
	for(int k=1; k<8; k++)
	{
		const float fTstDst = fabsf( uAlpha-fPalette[k] );
		if(fTstDst<fDist)
		{
			fDist = fTstDst;
			ibest_entry = k;
		}
	}
	return ibest_entry;
}


void BruteForceAlpha(const float alpha_block[], const int num_entry_instances[], const int iNumUniqueColors, int * piBestA0, int * piBestA1, float * pfCurBestDist);
void FastSearchAlpha(const float alpha_block[], const unsigned int alpha_block_integer[], const int num_entry_instances[], const int iNumUniqueColors, int * piBestA0, int * piBestA1, float * pfCurBestDist);

void GetBestAlpha(const unsigned char alpha_block_in[], const int ref_mode, int * piBestA0, int * piBestA1, float * pfCurBestDist)
{
	// sort colors and register duplicates
	unsigned int sorted_alpha_block[16];
	int i;
	for(i=0; i<16; i++) sorted_alpha_block[i] = alpha_block_in[i];	// make a copy
	QuickSort(sorted_alpha_block, 0, 16-1);
	int num_entry_instances[16];

	int entry = 0; num_entry_instances[0] = 1;
	for(i=1; i<16; i++)
	{
		if(sorted_alpha_block[i] == sorted_alpha_block[entry])
			++num_entry_instances[entry];
		else
		{
			++entry; 	// insert new unique entry (in place)
			sorted_alpha_block[entry] = sorted_alpha_block[i];
			num_entry_instances[entry] = 1;
		}
	}

	// amount of unique alphas
	const int iNumUniqueColors = entry+1;
	if(iNumUniqueColors < 3)
	{
		assert(iNumUniqueColors==1 || iNumUniqueColors==2);
		piBestA0[0] = sorted_alpha_block[0];
		piBestA1[0] = iNumUniqueColors==2?sorted_alpha_block[1]:piBestA0[0];
		if(piBestA0[0] < piBestA1[0])	// let's make it an A0 > A1 case
		{ const int iTmp=piBestA0[0]; piBestA0[0]=piBestA1[0]; piBestA1[0]=iTmp; }
		pfCurBestDist[0] = 0;
		return;
	}

	// array cast
	float alpha_block[16];
	for(i=0; i<iNumUniqueColors; i++) alpha_block[i] = (float) sorted_alpha_block[i];

	// check for brute force
	if(ref_mode == 3)
		BruteForceAlpha(alpha_block, num_entry_instances, iNumUniqueColors, piBestA0, piBestA1, pfCurBestDist);
	else
		FastSearchAlpha(alpha_block, sorted_alpha_block, num_entry_instances, iNumUniqueColors, piBestA0, piBestA1, pfCurBestDist);
}

bool ProcessSuggestion(const int iA0, const int iA1, const float unique_alphas_in[], const int num_entry_instances[], const int iNumUniqueColors, int * piBestA0, int * piBestA1, float * pfCurBestDist);

void BruteForceAlpha(const float alpha_block[], const int num_entry_instances[], const int iNumUniqueColors, int * piBestA0, int * piBestA1, float * pfCurBestDist)
{
	//int iEarlyOuts = 0;
	for(int iA0=0; iA0<256; iA0++)
		for(int iA1=0; iA1<256; iA1++)
		{
			/*const int iMax = iA0 < iA1 ? iA1 : iA0;
			const int iMin = iA0 < iA1 ? iA0 : iA1;

			bool bPotentialEarlyOut = false;
			int iTmp;
			if( iMax<iAlphaMin  )
			{
			iTmp = (iAlphaMin-iMax);
			if(iA0 <= iA1) iTmp = Min(iTmp, 255-iAlphaMax);
			bPotentialEarlyOut = true;
			}
			else if( iMin > iAlphaMax )
			{
			iTmp = (iMin-iAlphaMax);
			if(iA0 <= iA1) iTmp = Min(iTmp, iAlphaMin-0);
			bPotentialEarlyOut = true;
			}

			if(bPotentialEarlyOut)
			{
			iTmp *= iTmp;
			iTmp *= 16;	// min error
			if(iTmp > pfCurBestDist[0])
			{
			++iEarlyOuts;
			continue;
			}
			}*/

			// Evaluate Error
			ProcessSuggestion(iA0, iA1, alpha_block, num_entry_instances, iNumUniqueColors, piBestA0, piBestA1, pfCurBestDist);
		}
}


bool ProcessCodeBook(const float fCodeBook_in[], const int iCodeBookSize_in, const bool bPreserveCodeBookEnds, const float fUniqueAlphas[], const int num_entry_instances[], const int iNumUniqueAlphas, int * piBestA0, int * piBestA1, float * pfCurBestDist);

void FastSearchAlpha(const float alpha_block[], const unsigned int alpha_block_integer[], const int num_entry_instances[], const int iNumUniqueColors, int * piBestA0, int * piBestA1, float * pfCurBestDist)
{
	// find non zero minimum, and a max that is not 255
	int q=0;
	while(q<iNumUniqueColors && alpha_block_integer[q]==0) ++q;
	const int iAlphaMinNonZero = (q<iNumUniqueColors) ? alpha_block_integer[q] : 0;
	q=iNumUniqueColors-1;
	while(q>=0 && alpha_block_integer[q]==255) --q;
	const int iAlphaMaxNon255 = (q>=0) ? alpha_block_integer[q] : 255;

	// standard minimum and maximum
	const int iAlphaMin = alpha_block_integer[0];
	const int iAlphaMax = alpha_block_integer[iNumUniqueColors-1];


	// try a little VQ
	for(int z=0; z<2; z++)
		for(int k=2; k<8; k++)	// k segments
		{
			// case A0 > A1 first
			const float fAlphaRange = (const float) (iAlphaMax-iAlphaMin);
			const float fScale = z==0?(k / ((float) (k-1))):1;
			const float fInitialCodeBookSpan = fScale * fAlphaRange;
			float fCodeStart = iAlphaMin - ((fInitialCodeBookSpan-fAlphaRange)/2);
			float fCodeEnd = iAlphaMax + ((fInitialCodeBookSpan-fAlphaRange)/2);

			// out of range
			float fOffs = 0;
			if(fCodeStart<0) fOffs = 0-fCodeStart;
			if(fCodeEnd>255)
			{	const float fTmpOffs = fCodeEnd-255; if(fTmpOffs>fOffs) fOffs=fTmpOffs; }
			fCodeStart += fOffs; fCodeEnd -= fOffs;
			if(fCodeStart<0)fCodeStart=0; if(fCodeEnd>255)fCodeEnd=255;

			// fill up code book
			const int iCodeBookSize = k+1;
			float fCodeBook[8];
			for(int c=0; c<iCodeBookSize; c++)
			{
				const float fT = ((const float) c) / (iCodeBookSize-1);
				fCodeBook[c] = (fCodeEnd-fCodeStart)*fT + fCodeStart;
			}

			// process code book
			const bool bIsBetter = ProcessCodeBook(fCodeBook, iCodeBookSize, false, alpha_block, num_entry_instances, iNumUniqueColors, piBestA0, piBestA1, pfCurBestDist);

			// case A0 <= A1 second
			if(k<6)
			{
				const float fAlphaRange = (const float) (iAlphaMaxNon255-iAlphaMinNonZero);
				const float fInitialCodeBookSpan = fScale * fAlphaRange;
				float fCodeStart = iAlphaMinNonZero - ((fInitialCodeBookSpan-fAlphaRange)/2);
				float fCodeEnd = iAlphaMaxNon255 + ((fInitialCodeBookSpan-fAlphaRange)/2);

				// out of range
				float fOffs = 0;
				if(fCodeStart<0) fOffs = 0-fCodeStart;
				if(fCodeEnd>255)
				{	const float fTmpOffs = fCodeEnd-255; if(fTmpOffs>fOffs) fOffs=fTmpOffs; }
				fCodeStart += fOffs; fCodeEnd -= fOffs;
				if(fCodeStart<0)fCodeStart=0; if(fCodeEnd>255)fCodeEnd=255;

				if(fCodeStart <= fCodeEnd)	// only proceed if input makes sense
				{
					// fill up code book
					fCodeBook[0] = 0; fCodeBook[iCodeBookSize+1] = 255;
					for(int c=0; c<iCodeBookSize; c++)
					{
						const float fT = ((const float) c) / (iCodeBookSize-1);
						fCodeBook[c+1] = (fCodeEnd-fCodeStart)*fT + fCodeStart;
					}

					const int iActualCodeBookSize = k+1+2;	// +2, 0 and 255
					const bool bIsBetter = ProcessCodeBook(fCodeBook, iActualCodeBookSize, true, alpha_block, num_entry_instances, iNumUniqueColors, piBestA0, piBestA1, pfCurBestDist);
				}
			}
		}


		// a naive search
		int y;
		for(y=0; y<(iNumUniqueColors+2); y++)
			for(int x=0; x<(iNumUniqueColors+2); x++)
			{
				const int iA0 = x==0 ? 0 : (x==(iNumUniqueColors+1) ? 255 : alpha_block_integer[x-1]);
				const int iA1 = y==0 ? 0 : (y==(iNumUniqueColors+1) ? 255 : alpha_block_integer[y-1]);

				ProcessSuggestion(iA0, iA1, alpha_block, num_entry_instances, iNumUniqueColors, piBestA0, piBestA1, pfCurBestDist);
			}

			// try a little extra (this seldom finds anything better)
			for(int k=2; k<7; k++)
			{
				// make case A0 > A1 with k segments
				float fA1 = (float) iAlphaMin;
				const float fSegLen = ((float) (iAlphaMax-iAlphaMin))/k;
				float fA0 = fA1 + 7*fSegLen;
				while(fA0>255 && (fA1-fSegLen)>=0) { fA1-=fSegLen; fA0-=fSegLen; }

				if(fA0 <= 255)
					ProcessSuggestion(((int) fA0), ((int) fA1), alpha_block, num_entry_instances, iNumUniqueColors, piBestA0, piBestA1, pfCurBestDist);

				// make case A0 <= A1 with k segments
				if(k<5)
				{
					fA0 = (float) iAlphaMinNonZero;
					float fA1 = fA0 + 5*fSegLen;
					while(fA1>255 && (fA0-fSegLen)>=0) { fA1-=fSegLen; fA0-=fSegLen; }

					if(fA1 <= 255)
						ProcessSuggestion(((int) fA0), ((int) fA1), alpha_block, num_entry_instances, iNumUniqueColors, piBestA0, piBestA1, pfCurBestDist);
				}
			}



			// search for a while in local area
			const int Dim = 5;
			const int iLeftA0 = (piBestA0[0]-Dim)<0 ? 0 : (piBestA0[0]-Dim);
			const int iRightA0 = (piBestA0[0]+Dim)>255 ? 255 : (piBestA0[0]+Dim);
			const int iLeftA1 = (piBestA1[0]-Dim)<0 ? 0 : (piBestA1[0]-Dim);
			const int iRightA1 = (piBestA1[0]+Dim)>255 ? 255 : (piBestA1[0]+Dim);

			const int iOrgBestA0 = piBestA0[0];
			const int iOrgBestA1 = piBestA1[0];

			for(y=iLeftA1; y<=iRightA1; y++)
				for(int x=iLeftA0; x<=iRightA0; x++)
					ProcessSuggestion(x, y, alpha_block, num_entry_instances, iNumUniqueColors, piBestA0, piBestA1, pfCurBestDist);

			// make an additional search for a local minimum
			int iA0State = (piBestA0[0]-iOrgBestA0) > 0 ? 1 : (-1);	// false means left, true means right
			int iA1State = (piBestA1[0]-iOrgBestA1) > 0 ? 1 : (-1);	// false means left, true means right
			int iNextMove = Abs(piBestA1[0]-iOrgBestA1) > Abs(piBestA0[0]-iOrgBestA0) ? 1 : 0;
			bool bKeepGoing = true;
			int iTestA0 = piBestA0[0];
			int iTestA1 = piBestA1[0];
			int iNoProgessCounter = 0;
			while(bKeepGoing)
			{
				if(iNextMove!=0)
					iTestA1 += iA1State;
				else
					iTestA0 += iA0State;

				// clamp
				iTestA0 = iTestA0<0?0:(iTestA0>255?255:iTestA0);
				iTestA1 = iTestA1<0?0:(iTestA1>255?255:iTestA1);

				const bool bIsBetter = ProcessSuggestion(iTestA0, iTestA1, alpha_block, num_entry_instances, iNumUniqueColors, piBestA0, piBestA1, pfCurBestDist);
				if(bIsBetter) iNoProgessCounter = 0;
				else
				{
					iTestA0 = piBestA0[0]; iTestA1 = piBestA1[0];
					++iNoProgessCounter;

					if(iNoProgessCounter == 1) iNextMove=1-iNextMove;
					else if(iNoProgessCounter == 2) { if(iNextMove==0) iA0State = -iA0State; else iA1State = -iA1State; }
					else if(iNoProgessCounter == 3)
					{
						iNextMove=1-iNextMove;
						if(iNextMove==0) iA0State = -iA0State; else iA1State = -iA1State;
					}
					if(iNoProgessCounter == 4) bKeepGoing = false;		
				}
			}
}








bool ProcessSuggestion(const int iA0, const int iA1, const float unique_alphas_in[], const int num_entry_instances[], const int iNumUniqueColors, int * piBestA0, int * piBestA1, float * pfCurBestDist)
{
	bool bIsBetter = false;
	float fPalette[8];
	GetAlphas((const float) iA0, (const float) iA1, fPalette);

	float fSqErrorSum = 0;

	int iCurCandidateLocation = 0;
	for(int iCurAlphaIndex=0; iCurAlphaIndex<iNumUniqueColors; iCurAlphaIndex++)
	{
		const float fCurAlphaValue = unique_alphas_in[iCurAlphaIndex];
		float fBestDistToEntry = (fCurAlphaValue - fPalette[iCurCandidateLocation]);
		fBestDistToEntry *= fBestDistToEntry;
		bool bStillSearching = true;
		while((iCurCandidateLocation+1)<8 && bStillSearching)
		{
			float fDistToEntry = fCurAlphaValue - fPalette[iCurCandidateLocation+1];
			fDistToEntry *= fDistToEntry;
			if(fDistToEntry < fBestDistToEntry)
			{
				fBestDistToEntry = fDistToEntry;
				++iCurCandidateLocation;
			}
			else bStillSearching = false;
		}

		// sum up
		fSqErrorSum += (num_entry_instances[iCurAlphaIndex]*fBestDistToEntry);
	}

	// dunno what the hell the compiler is doing but this makes the compare 100% accurate.
	const unsigned int uVal1 = ((const unsigned int *) &fSqErrorSum)[0] & 0x7fffffff;
	const unsigned int uVal2 = ((const unsigned int *) pfCurBestDist)[0] & 0x7fffffff;
	if((pfCurBestDist[0]==-1) || (uVal1<uVal2) )
		//if(pfCurBestDist[0]==-1 || fSqErrorSum<pfCurBestDist[0])
	{
		piBestA0[0] = iA0; piBestA1[0] = iA1;
		pfCurBestDist[0] = fSqErrorSum;
		bIsBetter = true;
	}

	return bIsBetter;
}

inline void GetAlphas(const float fA0, const float fA1, float fColorsOut[])
{
	if(fA0 > fA1)
	{
		fColorsOut[0] = fA1; fColorsOut[1] = (fA1*6+fA0*1)/7;
		fColorsOut[2] = (fA1*5+fA0*2)/7; fColorsOut[3] = (fA1*4+fA0*3)/7;
		fColorsOut[4] = (fA1*3+fA0*4)/7; fColorsOut[5] = (fA1*2+fA0*5)/7;
		fColorsOut[6] = (fA1*1+fA0*6)/7; fColorsOut[7] = fA0;
	}
	else
	{
		fColorsOut[0] = 0; fColorsOut[1] = fA0;
		fColorsOut[2] = (fA0*4+fA1*1)/5; fColorsOut[3] = (fA0*3+fA1*2)/5;
		fColorsOut[4] = (fA0*2+fA1*3)/5; fColorsOut[5] = (fA0*1+fA1*4)/5;
		fColorsOut[6] = fA1; fColorsOut[7] = 255;
	}
}



bool StraightenOut( const int iLinearAlphas[], const int group_counts[], const bool bParityEven, int * pA0, int * pA1)
{
	bool bFoundMatch = false;
	// fLinearAlphas[] has 16 entries
	// group_counts[] has 8 entries
	const int num_cols = 16;

	const int a = group_counts[0];
	const int b = group_counts[1];
	const int c = group_counts[2];
	const int d = group_counts[3];
	const int e = group_counts[4];
	const int f = group_counts[5];
	const int g = group_counts[6];
	const int h = group_counts[7];

	const int A = bParityEven ? (b+4*c+9*d+16*e+25*f+36*g+49*h) : (25*b+16*c+9*d+4*e+1*f);
	const int B = bParityEven ? (6*b+10*c+12*d+12*e+10*f+6*g) : (4*c+6*d+6*e+4*f);
	const int C = bParityEven ? (g+4*f+9*e+16*d+25*c+36*b+49*a) : (1*c+4*d+9*e+16*f+25*g);

	const int iDeterminant = A*C - B*B;
	if(iDeterminant!=0)
	{
		const int iOffs0 = 0;
		const int iOffs1 = iOffs0+a;
		const int iOffs2 = iOffs1+b;
		const int iOffs3 = iOffs2+c;
		const int iOffs4 = iOffs3+d;
		const int iOffs5 = iOffs4+e;
		const int iOffs6 = iOffs5+f;
		const int iOffs7 = iOffs6+g;
		const int iOffs8 = iOffs7+h;	// group8 doesn't really exist

		int i;
		int iSum0 = 0; for(i=iOffs0; i<iOffs1; i++) iSum0+=iLinearAlphas[i];
		int iSum1 = 0; for(i=iOffs1; i<iOffs2; i++) iSum1+=iLinearAlphas[i];
		int iSum2 = 0; for(i=iOffs2; i<iOffs3; i++) iSum2+=iLinearAlphas[i];
		int iSum3 = 0; for(i=iOffs3; i<iOffs4; i++) iSum3+=iLinearAlphas[i];
		int iSum4 = 0; for(i=iOffs4; i<iOffs5; i++) iSum4+=iLinearAlphas[i];
		int iSum5 = 0; for(i=iOffs5; i<iOffs6; i++) iSum5+=iLinearAlphas[i];
		int iSum6 = 0; for(i=iOffs6; i<iOffs7; i++) iSum6+=iLinearAlphas[i];
		int iSum7 = 0; for(i=iOffs7; i<iOffs8; i++) iSum7+=iLinearAlphas[i];

		const int K1 = bParityEven ? (1*7*iSum1+2*7*iSum2+3*7*iSum3+4*7*iSum4+5*7*iSum5+6*7*iSum6+7*7*iSum7) : (5*5*iSum1+4*5*iSum2+3*5*iSum3+2*5*iSum4+1*5*iSum5);
		const int K2 = bParityEven ? (7*7*iSum0+6*7*iSum1+5*7*iSum2+4*7*iSum3+3*7*iSum4+2*7*iSum5+1*7*iSum6) : (1*5*iSum2+2*5*iSum3+3*5*iSum4+4*5*iSum5+5*5*iSum6);

		float a0 = (float) ( (C*K1 - B*K2)/((double) iDeterminant) );
		float a1 = (float) ( ((-B)*K1 + A*K2)/((double) iDeterminant) );


		// find the correct location
		const float fDiff = a0-a1;
		const float fSegLen = (fDiff<0?(-fDiff):fDiff) / (bParityEven?7:5);
		while(a1>=255.5 && (a0-fSegLen)>=0 ) { a0-=fSegLen; a1-=fSegLen; }
		while(a0>=255.5 && (a1-fSegLen)>=0 ) { a0-=fSegLen; a1-=fSegLen; }
		while(a1<0 && (a0+fSegLen)<255.5 ) { a0+=fSegLen; a1+=fSegLen; }
		while(a0<0 && (a1+fSegLen)<255.5 ) { a0+=fSegLen; a1+=fSegLen; }

		int iA0 = (int) (a0+0.5f);
		int iA1 = (int) (a1+0.5f);
		iA0 = iA0<0 ? 0 : (iA0>255?255:iA0);	// Clamp
		iA1 = iA1<0 ? 0 : (iA1>255?255:iA1);

		// correct value
		if( (iA0 > iA1)^bParityEven )
		{ const int iTmp=iA0; iA0=iA1; iA1=iTmp; }

		if((iA0==iA1) && bParityEven)
		{
			if(iA1 > 0) --iA1;
			else ++iA0;
		}

		pA0[0]=iA0; pA1[0]=iA1;

		bFoundMatch = true;
	}

	return bFoundMatch;
}

//float Squared(const float x) { return x*x; }

bool ProcessCodeBook(const float fCodeBook_in[], const int iCodeBookSize_in, const bool bPreserveCodeBookEnds, const float fUniqueAlphas[], const int num_entry_instances[], const int iNumUniqueAlphas, int * piBestA0, int * piBestA1, float * pfCurBestDist)
{
	// skip entry zero if bPreserveCodeBookEnds
	assert(iCodeBookSize_in>1 || (bPreserveCodeBookEnds && iCodeBookSize_in>2));
	//const float fSegLen = bPreserveCodeBookEnds ? (fCodeBook_in[2]-fCodeBook_in[1]) : (fCodeBook_in[1]-fCodeBook_in[0]);

	assert(iCodeBookSize_in>0 && iCodeBookSize_in<=8);
	int iCodeBookSize = iCodeBookSize_in;
	float fprevCodeBook[8];
	for(int c=0; c<iCodeBookSize; c++) fprevCodeBook[c] = fCodeBook_in[c];

	int entry_group_ass[16];	// only need iNumUniqueAlphas in length

	// make a linear version
	int iLinearAlphas[16];
	int index = 0;
	for(int c=0; c<iNumUniqueAlphas; c++)
		for(int i=0; i<num_entry_instances[c]; i++)
		{
			assert(index<16);
			iLinearAlphas[index++] = (int) (fUniqueAlphas[c]+0.5f);
		}

	assert(index == 16);
	// if bPreserveCodeBookEnds, then A0 <= A1
	int iA0, iA1;
	bool bFoundValues = false;
	for(int k=0; k<3; k++)
	{
		const int iMAX_ITERATIONS = 10;
		int index = 0;
		bool bInProcess = true;
		int prev_contributions[8];
		while(bInProcess)
		{
			float fNextCodeBook[8];
			int contributions[8];
			int c;
			for(c=0; c<iCodeBookSize; c++)
			{
				fNextCodeBook[c] = 0;
				contributions[c] = 0;
			}
			if(bPreserveCodeBookEnds) fNextCodeBook[iCodeBookSize-1] = 255;

			int iCurCandidateLocation = 0;
			for(int iCurAlphaIndex=0; iCurAlphaIndex<iNumUniqueAlphas; iCurAlphaIndex++)
			{
				const float fCurAlphaValue = fUniqueAlphas[iCurAlphaIndex];
				float fBestDistToEntry = (fCurAlphaValue - fprevCodeBook[iCurCandidateLocation]);
				fBestDistToEntry *= fBestDistToEntry;
				bool bStillSearching = true;
				while((iCurCandidateLocation+1)<iCodeBookSize && bStillSearching)
				{
					float fDistToEntry = fCurAlphaValue - fprevCodeBook[iCurCandidateLocation+1];
					fDistToEntry *= fDistToEntry;
					if(fDistToEntry < fBestDistToEntry)
					{
						fBestDistToEntry = fDistToEntry;
						++iCurCandidateLocation;
					}
					else bStillSearching = false;
				}

				// update codebook entry
				if((!bPreserveCodeBookEnds) || (iCurCandidateLocation>0 && iCurCandidateLocation<(iCodeBookSize-1)))
					fNextCodeBook[iCurCandidateLocation] += (fCurAlphaValue*num_entry_instances[iCurAlphaIndex]);
				contributions[iCurCandidateLocation] += num_entry_instances[iCurAlphaIndex];
				entry_group_ass[iCurAlphaIndex] = iCurCandidateLocation;
			}

			// average
			const int iStart = bPreserveCodeBookEnds ? 1 : 0;
			const int iEnd = iCodeBookSize - (bPreserveCodeBookEnds ? 1 : 0);
			for(int s=iStart; s<iEnd; s++)
				if(contributions[s]!=0)
					fNextCodeBook[s] /= contributions[s];
				else
					fNextCodeBook[s] = fprevCodeBook[s];


			for(int s=0; s<iCodeBookSize; s++)
				assert( fNextCodeBook[s] >= 0 );

			// bubble sort is approximately linear for almost sorted arrays
			for(int i=1; i<iCodeBookSize; i++)
			{
				int j=i;
				while(j>0 && (((unsigned int *) fNextCodeBook)[j-1]>((unsigned int *) fNextCodeBook)[j]) )
				{
					// swap entries
					const float fTmp = fNextCodeBook[j-1]; fNextCodeBook[j-1]=fNextCodeBook[j]; fNextCodeBook[j]=fTmp;
					const int iTmp = contributions[j-1]; contributions[j-1]=contributions[j]; contributions[j]=iTmp;
					--j;	// next
				}
			}

			for(int s=1; s<iCodeBookSize; s++)
				assert( fNextCodeBook[s-1] <= fNextCodeBook[s] );




			++index;
			if(index == iMAX_ITERATIONS) bInProcess = false;

			// check for termination
			if(bInProcess && index>1)
			{
				int i=0;
				bool bMightBeSame = true;
				while(i<iCodeBookSize && bMightBeSame)	// look for no cluster change
				{
					if(prev_contributions[i] != contributions[i]) bMightBeSame = false;
					else ++i;
				}
				if(bMightBeSame) bInProcess = false;
			}

			// copy
			for(int c=0; c<iCodeBookSize; c++)
			{
				fprevCodeBook[c] = fNextCodeBook[c];
				prev_contributions[c] = contributions[c];
			}
		}

		// straighten out code book
		const bool bParityEven = !bPreserveCodeBookEnds;

		// convert to 8 element group
		int group_counts[8];
		const int iEnd = bParityEven ? iCodeBookSize : (iCodeBookSize-1);
		for(int c=0; c<iEnd; c++) group_counts[c] = prev_contributions[c];
		for(int c=iEnd; c<8; c++) group_counts[c] = 0;	// clear groups
		if(!bParityEven) group_counts[7] = prev_contributions[iEnd];	// 255 case

		bFoundValues = StraightenOut( iLinearAlphas, group_counts, bParityEven, &iA0, &iA1);
		if(bFoundValues)
		{
			iCodeBookSize = 8;
			GetAlphas((float) iA0, (float) iA1, fprevCodeBook);

			/*float fAlphaPal[8];
			GetAlphas((float) iA0, (float) iA1, fAlphaPal);

			unsigned char uPalEntriesWeKeep = bPreserveCodeBookEnds?0x81:0;
			for(int c=0; c<iCodeBookSize; c++)
			{
			float fEntryDists[] = {0,0,0,0,0,0,0,0};
			for(int i=0; i<iNumUniqueAlphas; i++)
			{
			if(entry_group_ass[i]==c)
			for(int p=0; p<8; p++)
			fEntryDists[p] += (num_entry_instances[i]*Squared(fAlphaPal[p]-fUniqueAlphas[i]));
			}

			float fMin = fEntryDists[0];
			int iBestEntry = 0;
			for(int p=1; p<8; p++)
			if(fEntryDists[p]<fMin)
			{ fMin=fEntryDists[p]; iBestEntry=p; }

			uPalEntriesWeKeep |= (1<<iBestEntry);
			}

			int taken = 0;
			int p;
			for(p=0; p<8; p++)
			if((uPalEntriesWeKeep&(1<<p))!=0)
			++taken;
			//assert(taken<=iCodeBookSize);
			if(taken>iCodeBookSize) iCodeBookSize = taken;
			int entries_missing = iCodeBookSize-taken;
			int index = 0;
			for(p=0; p<8; p++)
			{
			const bool bTest = (uPalEntriesWeKeep&(1<<p))!=0;
			if(bTest || entries_missing>0)
			{
			fprevCodeBook[index++] = fAlphaPal[p];
			if(!bTest) --entries_missing;
			}
			}
			assert(index == iCodeBookSize);*/

			for(int s=1; s<iCodeBookSize; s++)
				assert( fprevCodeBook[s-1] <= fprevCodeBook[s] );
		}
	}

	if(!bFoundValues) return false;
	else return ProcessSuggestion(iA0, iA1, fUniqueAlphas, num_entry_instances, iNumUniqueAlphas, piBestA0, piBestA1, pfCurBestDist);
}