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

#include "tgaloader.h"
#include <assert.h>


// bitmaps are received and given in the traditional zig zag pattern.
// first 3 bytes represent the upper-left pixel and the last 3 bytes
// represent the lower-right pixel. The pixels are ordered as r, g, b, r, ...

bool ReadTGA(unsigned char ** ppbitmap, int * piWidth, int * piHeight, int * piNrChannels, FILE * fptr_in)
{
	bool bSuccess = false;
	piWidth[0] = 0; piHeight[0] = 0; piNrChannels[0] = 0;
	unsigned char header[18];
	fread(header, 18, 1, fptr_in);

	// 

	// read ID field.
	const int iDescSize = header[0];
	for(int c=0; c<iDescSize; c++)
	{ unsigned char tmp; fread(&tmp, 1, 1, fptr_in); }

	// fetch dimensions
	const int dimX = ((header[13]&0xff)<<8) | (header[12]&0xff);
	const int dimY = ((header[15]&0xff)<<8) | (header[14]&0xff);

	const bool bXFlip = (header[17]&0x10)!=0;
	const bool bYFlip = (header[17]&0x20)!=0;

	const int bits_per_pixel = header[16];
	//assert(bits_per_pixel==8 || bits_per_pixel==24);

	const int iNrChannelsIn = bits_per_pixel / 8;

	const int iNrBytes = dimX*dimY*iNrChannelsIn;
	ppbitmap[0] = new unsigned char[iNrBytes];
	if(ppbitmap[0]!=NULL)
	{
		piWidth[0] = dimX; 	piHeight[0] = dimY;
		piNrChannels[0] = iNrChannelsIn;

		unsigned char * tmp_buf = new unsigned char[iNrChannelsIn*dimX*dimY];

		// speed read
		if(tmp_buf!=NULL) fread(tmp_buf, 1, iNrChannelsIn*dimX*dimY, fptr_in);

		int index = 0;
		for(int y=0; y<dimY; y++)
			for(int x=0; x<dimX; x++)
			{
				const int iY = bYFlip ? y : ((dimY-1)-y);
				const int iX = bXFlip ? ((dimX-1)-x) : x;
				unsigned char * pBuf = ppbitmap[0]+iNrChannelsIn*(iY*dimX+iX);

				if(tmp_buf==NULL)
				{
					// slow read
					for(int c=0; c<iNrChannelsIn; c++)
					{
						const int index = c==2 ? 0 : (c==0 ? 2 : c);
						fread(&pBuf[index], 1, 1, fptr_in);
					}
				}
				else
				{
					// cached read
					for(int c=0; c<iNrChannelsIn; c++)
					{
						const int index2 = c==2 ? 0 : (c==0 ? 2 : c);
						pBuf[index2]=tmp_buf[index++];
					}
				}
				//assert(iNrChannelsIn>0 && iNrChannelsIn<=3);
			}		

		bSuccess = true;
		if(tmp_buf!=NULL) { delete [] tmp_buf; tmp_buf = NULL; }
	}

	return bSuccess;
}

void WriteTGA(FILE * fptr_out, const unsigned char bitmap[], const int iWidth, const int iHeight, const int iNrChannels);

void WriteTGA_24Bit(FILE * fptr_out, const unsigned char bitmap[], const int iWidth, const int iHeight)
{
	WriteTGA(fptr_out, bitmap, iWidth, iHeight, 3);
}

void WriteTGA_32Bit(FILE * fptr_out, const unsigned char bitmap[], const int iWidth, const int iHeight)
{
	WriteTGA(fptr_out, bitmap, iWidth, iHeight, 4);
}


void WriteTGA(FILE * fptr_out, const unsigned char bitmap[], const int iWidth, const int iHeight, const int iNrChannels)
{
	unsigned char header[18] = {
		0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, (iWidth>>0)&0xff, (iWidth>>8)&0xff, (iHeight>>0)&0xff, (iHeight>>8)&0xff, 8*iNrChannels, 8*0};
	fwrite(header, 18, 1, fptr_out);

	unsigned char * tmp_buf = new unsigned char[iNrChannels*iWidth*iHeight];
	if(tmp_buf!=NULL)
	{
		// speed write
		int index = 0;
		for(int y=0; y<iHeight; y++)
			for(int x=0; x<iWidth; x++)
			{
				const int iY = (iHeight-1)-y;	// flip
				const int iOffs = iNrChannels*(iY*iWidth+x);
				for(int c=0; c<iNrChannels; c++)
				{
					const int index2 = c==2 ? 0 : (c==0 ? 2 : c);
					tmp_buf[index++] = bitmap[iOffs+index2];
				}
			}
		fwrite(tmp_buf, 1, iNrChannels*iWidth*iHeight, fptr_out);
		delete [] tmp_buf; tmp_buf = NULL;
	}
	else
	{
		// slow write
		for(int y=0; y<iHeight; y++)
			for(int x=0; x<iWidth; x++)
			{
				const int iY = (iHeight-1)-y;	// flip
				const int iOffs = iNrChannels*(iY*iWidth+x);
				for(int c=0; c<iNrChannels; c++)
				{
					const int index = c==2 ? 0 : (c==0 ? 2 : c);
					fwrite(&bitmap[iOffs+index], 1, 1, fptr_out);
				}
			}
	}
}