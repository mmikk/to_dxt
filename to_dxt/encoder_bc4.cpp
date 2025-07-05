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

#include "encoder_bc4.h"
#include <assert.h>

float EvalError(unsigned short uLeft, unsigned short uRight, unsigned char Indices[], const uint16 block_in[], bool bTestingEndPointsBW);
bool BuildIndices(unsigned short uLeft, unsigned short uRight, unsigned char Indices[], const uint16 block_in[], bool bTestingEndPointsBW);
void CorrectEndPoints(unsigned short * puLeft, unsigned short *puRight);
void DecodeBC4(uint16 uResultingPixels[], const uint64 block_bc4);
int LinearIndexFromBC4(const int idx0, const bool bTestingEndPointsBW);

void ConvertBlockToBC4(const uint16 block_in[], uint64 block_out[], uint16 uResultingPixels[])
{
	// build bound
	uint16 uMinReg=block_in[0], uMaxReg=block_in[0];
	for(int i=1; i<16; i++)
	{
		if(uMinReg>block_in[i]) uMinReg=block_in[i];
		else if(uMaxReg<block_in[i]) uMaxReg=block_in[i];
	}

	// build bound BW
	uint16 uMinBW=0xffff, uMaxBW=0;
	for(int i=0; i<16; i++)
	{
		const uint16 uX = block_in[i];
		
		if(uX!=0 && uX!=0xffff)
		{
			if(uMinBW>uX) uMinBW=uX;
			if(uMaxBW<uX) uMaxBW=uX;
		}
	}
	const bool isValidCaseBW = uMinBW<=uMaxBW;
	 

	//
	int iStepReg = (15*(uMaxReg-uMinReg)+50)/100;
	int iStepBW = (15*(uMaxBW-uMinBW)+50)/100;

	// record best candidate
	int iBestItt = -1;
	float fSmallestError = -1;
	uint16 uLeftBest=0, uRightBest=0;
	bool bBestCandidateHasEndPointsBW = false;
	unsigned char best_indices[16];

	// search for suiter
	for(int q=0; q<(isValidCaseBW ? 2 : 1); q++)
	{
		bool bTestingEndPointsBW = q!=0;
		const int fix_pnt = bTestingEndPointsBW ? 6 : 8;//bTesting8BitEndPoints ? 128 : 64;
		const int fix_pnt_dec = fix_pnt-1;

		const int iStep = bTestingEndPointsBW ? iStepBW : iStepReg;
		const uint16 uMin = bTestingEndPointsBW ? uMinBW : uMinReg;
		const uint16 uMax = bTestingEndPointsBW ? uMaxBW : uMaxReg;

		for(int j=-1; j<=1; j++)
		{
			int iMinBound = ((int) uMin)+j*iStep;
			for(int i=-1; i<=1; i++)
			{
				if(fSmallestError==0.0f) break;

				int iMaxBound = ((int) uMax)+i*iStep;
				assert(iMaxBound>=iMinBound);
				if(iMaxBound<iMinBound)
				{
					const int iTmp = iMaxBound;
					iMaxBound = iMinBound;
					iMinBound = iTmp;
				}

				// clamp to 16 bit
				uint16 uLeft = iMinBound>0xffff ? 0xffff : (iMinBound<0 ? 0 : iMinBound);
				uint16 uRight = iMaxBound>0xffff ? 0xffff : (iMaxBound<0 ? 0 : iMaxBound);

				unsigned char Indices[16];

				/*if(bTesting8BitEndPoints)*/ CorrectEndPoints(&uLeft, &uRight);
				BuildIndices(uLeft, uRight, Indices, block_in, bTestingEndPointsBW);

				// initial naive result
				float fError = EvalError(uLeft, uRight, Indices, block_in, bTestingEndPointsBW);	
				if( fError>=0 && (fError<fSmallestError || fSmallestError<0) )
				{
					uLeftBest = uLeft; uRightBest = uRight;
					for(int p=0; p<16; p++) best_indices[p]=Indices[p];
					bBestCandidateHasEndPointsBW = bTestingEndPointsBW;
					fSmallestError = fError;
					iBestItt = -1;
				}
				
				// proceed with iterations
				bool earlyOut = fSmallestError==0.0f;
				//for(int k=0; k<20; k++)
				for(int k=0; k<12; k++)
				{
					if(earlyOut) break;
					//if(fSmallestError==0.0f) break;

					//uLeft  * SUM_l (63-Indices[l])*(63-Indices[l]) + uRight * SUM_l (63-Indices[l])*Indices[l] =	SUM_l (63-Indices[l])*63*block_in[l]) -(63-Indices[l])*31
					//uLeft  * SUM_l Indices[l](63-Indices[l])		 + uRight * SUM_l Indices[l]*Indices[l]	     =	SUM_l  Indices[l]*63*block_in[l] -Indices[l]*31

					// | A B | |l| | k0 |
					// | B C |*|r|=| k1 |
					int iA = 0, iB=0, iC=0;		// used to be doubles
					int iK0 = 0, iK1 = 0;		// used to be doubles
					int64 iK064 = 0, iK164 = 0;		// using these for validation only
					for(int l=0; l<16; l++)
					{
						const int idx0 = (const int) Indices[l];
						//if(!bTestingEndPointsBW || idx0<fix_pnt)
						if(idx0<fix_pnt)
						{
							const int idx = LinearIndexFromBC4(idx0, bTestingEndPointsBW);		// used to be double
							const int iInvidx = fix_pnt_dec-idx;		// used to be double

							iA += (iInvidx*iInvidx);
							iB += (idx*iInvidx);
							iC += (idx*idx);
							const auto tmp0 = (-iInvidx*(fix_pnt_dec>>1) + iInvidx*fix_pnt_dec*((int) block_in[l]));
							const auto tmp1 = (-idx*(fix_pnt_dec>>1) + idx*fix_pnt_dec*((int) block_in[l]));
							iK0 += tmp0; iK1 += tmp1;
							iK064 += tmp0; iK164 += tmp1;		// using these for validation only
						}
					}
					assert(iK0==iK064 && iK1==iK164);

					const int iDet = iA*iC - iB*iB;		// used to be double
					if(iDet!=0)
					{
						// need int64 bit or we get overflow
						const int64 iMinBound = ((iC*((int64) iK0) - iB*((int64) iK1))/iDet);
						const int64 iMaxBound = ((-iB*((int64) iK0) + iA*((int64) iK1))/iDet);

						uint16 uLeftOrg = uLeft, uRightOrg = uRight;
						// snap to
						uLeft = iMinBound>0xffff ? 0xffff : (iMinBound<0 ? 0 : iMinBound);
						uRight = iMaxBound>0xffff ? 0xffff : (iMaxBound<0 ? 0 : iMaxBound);
						assert(uLeft<=uRight);
						
						/*if(bTesting8BitEndPoints)*/ CorrectEndPoints(&uLeft, &uRight);

						// determine error
						bool allIndicesMatched = BuildIndices(uLeft, uRight, Indices, block_in, bTestingEndPointsBW);
						fError = EvalError(uLeft, uRight, Indices, block_in, bTestingEndPointsBW);
						const bool earlyOutOrg = earlyOut;
						earlyOut = k>0 && uLeft==uLeftOrg && uRight==uRightOrg;
						assert((uLeft==uLeftOrg && uRight==uRightOrg) || (!earlyOut));
						assert(earlyOut || earlyOut==earlyOutOrg);
						assert(!(uLeft == uLeftOrg && uRight == uRightOrg) || allIndicesMatched);
					
						if( fError>=0 && (fError<fSmallestError || fSmallestError<0) )
						{
							uLeftBest = uLeft; uRightBest = uRight;
							for(int p=0; p<16; p++) best_indices[p]=Indices[p];
							bBestCandidateHasEndPointsBW = bTestingEndPointsBW;
							fSmallestError = fError;
							iBestItt = k;
							earlyOut |= fSmallestError==0.0f;
						}
					}
					else earlyOut=true;
				}
			}
		}
	}

	const int fix_pnt = bBestCandidateHasEndPointsBW ? 6 : 8;
	const int fix_pnt_dec = fix_pnt-1;

	const bool mustSwap = (uLeftBest>uRightBest && bBestCandidateHasEndPointsBW) || (uLeftBest<=uRightBest && (!bBestCandidateHasEndPointsBW));
	if(mustSwap)
	{
		for(int p=0; p<16; p++)
		{
			const unsigned char idx = best_indices[p];
			if(idx<fix_pnt)
			{
				best_indices[p] = idx==1 ? 0 : (idx==0 ? 1 : ((fix_pnt_dec+2)-idx));
			}
		}
		const unsigned short uTmp = uLeftBest;
		uLeftBest=uRightBest; uRightBest=uTmp;
	}

	// write out encoded block
#if 0
	((uint8 *) block_out)[0] = (uint8) (uLeftBest>>8);		// fill up first short
	((uint8 *) block_out)[1] = (uint8) (uRightBest>>8);

	// insert bits
	for(int q=0; q<2; q++)
	{
		const int iByteOffs = q*3;
		const int iIOffs = q*8;
		((uint8 *) block_out)[2+iByteOffs] = ((best_indices[2+iIOffs]&0x3)<<6) | (best_indices[1+iIOffs]<<3) | (best_indices[0+iIOffs]<<0);
		((uint8 *) block_out)[3+iByteOffs] = ((best_indices[5+iIOffs]&0x1)<<7) | (best_indices[4+iIOffs]<<4) | (best_indices[3+iIOffs]<<1) | (best_indices[2+iIOffs]>>2);
		((uint8 *) block_out)[4+iByteOffs] = (best_indices[7+iIOffs]<<5) | (best_indices[6+iIOffs]<<2) | (best_indices[5+iIOffs]>>1);
	}
#else
	block_out[0] = (uRightBest&0xff00) | ((uLeftBest>>8)&0xff);

	for(int i=0; i<16; i++)
	{
		assert(best_indices[i]>=0 && best_indices[i]<8);
		block_out[0] |= ((((uint64) best_indices[i])&0x7)<<(i*3+16));
	}
#endif

	if(uResultingPixels!=NULL) 
	{
		DecodeBC4(uResultingPixels, block_out[0]);
		if(fSmallestError==0.0f)
		{
			for(int q=0; q<16; q++) assert(uResultingPixels[q]==block_in[q]);
		}
	}
}

void DecodeBC4(uint16 uResultingPixels[], const uint64 block_bc4)
{
	unsigned short uPnt0, uPnt1;

	uPnt0 = (block_bc4&0xff); uPnt0 = (uPnt0<<8) | uPnt0;
	uPnt1 = ((block_bc4>>8)&0xff); uPnt1 = (uPnt1<<8) | uPnt1;

	const bool hasEndPointsBW = !(uPnt0>uPnt1);
	const int fix_pnt = hasEndPointsBW ? 6 : 8;
	const int fix_pnt_dec = fix_pnt-1;

	unsigned char indices[16];
	for(int i=0; i<16; i++)
		indices[i]=(unsigned char) ((block_bc4>>(i*3+16))&0x7);

	// use decompressed indices to build tile colors
	for(int p=0; p<16; p++)
	{
		const int idx0 = indices[p];
		assert(idx0>=0 && idx0<8);

		if(hasEndPointsBW && (idx0==6 || idx0==7))
		{
			uResultingPixels[p] = idx0==6 ? 0 : 0xffff;
		}
		else
		{
			// reverse bc4 index order to a sorted ordering
			const int idx = LinearIndexFromBC4(idx0, hasEndPointsBW);
			uResultingPixels[p] = (uint16) ((uPnt1*idx + uPnt0*(fix_pnt_dec-idx) + (fix_pnt_dec>>1))/fix_pnt_dec);
		}
	}
}

int LinearIndexFromBC4(const int idx0, const bool bTestingEndPointsBW)
{
	const int fix_pnt = bTestingEndPointsBW ? 6 : 8;
	const int fix_pnt_dec = fix_pnt-1;

	const int idx = idx0==0 ? 0 : (idx0==1 ? fix_pnt_dec : (idx0<fix_pnt ? (idx0-1) : idx0));
	assert(idx0<fix_pnt || bTestingEndPointsBW);

	return idx;
}

int IndexBC4FromLinearIndex(const int idx0, const bool bTestingEndPointsBW)
{
	const int fix_pnt = bTestingEndPointsBW ? 6 : 8;
	const int fix_pnt_dec = fix_pnt-1;

	// remap to silly bc4 index layout order: red0, red1, ...
	const int idx = idx0==0 ? 0 : (idx0==fix_pnt_dec ? 1 : (idx0<fix_pnt_dec ? (idx0+1) : idx0));
	assert(idx<fix_pnt || bTestingEndPointsBW);

	return idx;
}

float EvalError(unsigned short uLeft, unsigned short uRight, unsigned char Indices[], const uint16 block_in[], bool bTestingEndPointsBW)
{
	// place new min and max bound to minimize the error
	//(  err = (uLeft*(63-Indices[l]) + uRight*Indices[l] + 31)  -  63*block_in[l]   )^2

	double dError = 0;

	const int fix_pnt = bTestingEndPointsBW ? 6 : 8;//bTesting16BitEndPoints ? 64 : 128;
	const int fix_pnt_dec = fix_pnt-1;
	//(  err = (uLeft*(63-Indices[l]) + uRight*Indices[l] + 31)  -  63*block_in[l]   )^2
	for(int l=0; l<16; l++)
	{
		const int idx0 = Indices[l];

		double delta;
		if(bTestingEndPointsBW && (idx0==6 || idx0==7))
		{
			delta = idx0==6 ? block_in[l] : (0xffff-block_in[l]);
		}
		else
		{
			// reverse bc4 index order to a sorted ordering
			//const int idx = idx0==0 ? 0 : (idx0==1 ? fix_pnt_dec : (idx0-1));
			const int idx = LinearIndexFromBC4(idx0, bTestingEndPointsBW);
			
			const int iInvidx = fix_pnt_dec-idx;
			//delta = ( ((double) uLeft)*iInvidx + ((double) uRight)*idx + (fix_pnt_dec>>1)  -  fix_pnt_dec*((int) block_in[l]) )/fix_pnt_dec;
			const int iDelta = ( ((int) uLeft)*iInvidx + ((int) uRight)*idx + (fix_pnt_dec>>1)  -  fix_pnt_dec*((int) block_in[l]) )/fix_pnt_dec;
			delta = (double) iDelta;
		}

		dError += delta*delta;
	}
	return (float) dError;
}

bool BuildIndices(unsigned short uLeft, unsigned short uRight, unsigned char Indices[], const uint16 block_in[], bool bTestingEndPointsBW)
{
	bool indicesUnchanged = true;

	const int fix_pnt = bTestingEndPointsBW ? 6 : 8;//bTesting16BitEndPoints ? 64 : 128;
	const int fix_pnt_dec = fix_pnt-1;
	const int iDiff = (int) (uRight - uLeft);

	unsigned char idx;
#if 0
	for(int l=0; l<16; l++)
	{
		if(bTestingEndPointsBW && (block_in[l]==0 || block_in[l]==0xffff))
		{
			idx = block_in[l]==0xffff ? 7 : 6;
		}
		else
		{
			int Index=iDiff>0 ? ((fix_pnt_dec*(((int) block_in[l])-uLeft) + (iDiff/2)) / iDiff) : 0;
			idx = Index>fix_pnt_dec ? fix_pnt_dec : (Index<0 ? 0 : Index);

			// remap to silly bc4 index layout order: red0, red1, ...
			idx = idx==0 ? 0 : (idx==fix_pnt_dec ? 1 : (idx+1));
			assert(idx<fix_pnt);
		}

		indicesUnchanged = indicesUnchanged && Indices[l] == idx;
		Indices[l] = idx;
	}
#else
	for(int l=0; l<16; l++)
	{
		int Index=iDiff>0 ? ((fix_pnt_dec*(((int) block_in[l])-uLeft) + (iDiff/2)) / iDiff) : 0;
		idx = Index>fix_pnt_dec ? fix_pnt_dec : (Index<0 ? 0 : Index);

		if(bTestingEndPointsBW)
		{
			const int iInvidx = fix_pnt_dec-idx;
			const int iDelta = ( ((int) uLeft)*iInvidx + ((int) uRight)*idx + (fix_pnt_dec>>1)  -  fix_pnt_dec*((int) block_in[l]) )/fix_pnt_dec;

			const unsigned short uDelta = (const unsigned short) (iDelta<0 ? (-iDelta) : iDelta);
			const unsigned short uWhiteError = (const unsigned short) (0xffff-block_in[l]);
			const unsigned short uBlackError = block_in[l];

			if(uWhiteError<uDelta && uWhiteError<uBlackError) idx = 7;
			else if(uBlackError<uDelta && uBlackError<uWhiteError) idx = 6;
		}

		idx = IndexBC4FromLinearIndex(idx, bTestingEndPointsBW);

		indicesUnchanged = indicesUnchanged && Indices[l] == idx;
		Indices[l] = idx;
	}
#endif

	return indicesUnchanged;
}

void CorrectEndPoints(unsigned short * puLeft, unsigned short *puRight)
{
	unsigned short uL = *puLeft;
	unsigned short uR = *puRight;

	uL >>= 8; uR >>= 8;
	if(((uL<<8)|uL)>puLeft[0] && uL>0) --uL;
	if(((uR<<8)|uR)<puRight[0] && uR<0xff) ++uR;

	*puLeft = ((uL<<8)|uL);
	*puRight = ((uR<<8)|uR);
}