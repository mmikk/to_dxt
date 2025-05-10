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

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include "tgaloader.h"
#include "miniutils.h"
#include "encoder_dxt.h"
#include "encoder_alphadxt5.h"


// these should be correct
#define DDSD_CAPS					0x00000001
#define DDSD_HEIGHT					0x00000002
#define DDSD_WIDTH					0x00000004
#define DDSD_PITCH					0x00000008
#define DDSD_PIXELFORMAT			0x00001000
#define DDSD_MIPMAPCOUNT			0x00020000
#define DDSD_LINEARSIZE				0x00080000
#define DDSD_DEPTH					0x00800000

typedef unsigned int uint32;

// Direct3D 9
#define DDPF_ALPHAPIXELS			0x00000001
#define DDPF_FOURCC					0x00000004
#define DDPF_RGB					0x00000040

struct DDPIXELFORMAT
{
    uint32 dwSize;
    uint32 dwFlags;
    uint32 dwFourCC;
    uint32 dwRGBBitCount;
    uint32 dwRBitMask, dwGBitMask, dwBBitMask;
    uint32 dwRGBAlphaBitMask;
};


// must be used on any file that contains more than one surface (a mipmap, a cubic environment map, or volume texture).
#define DDSCAPS_COMPLEX				0x00000008

// Required
#define DDSCAPS_TEXTURE				0x00001000

// should be used for a mipmap.
#define DDSCAPS_MIPMAP				0x00400000

#define DDSCAPS2_CUBEMAP			0x00000200
#define DDSCAPS2_CUBEMAP_POSITIVEX	0x00000400
#define DDSCAPS2_CUBEMAP_NEGATIVEX	0x00000800
#define DDSCAPS2_CUBEMAP_POSITIVEY	0x00001000
#define DDSCAPS2_CUBEMAP_NEGATIVEY	0x00002000
#define DDSCAPS2_CUBEMAP_POSITIVEZ	0x00004000
#define DDSCAPS2_CUBEMAP_NEGATIVEZ	0x00008000
#define DDSCAPS2_VOLUME				0x00200000

struct DDSCAPS2
{
    uint32 dwCaps1;
    uint32 dwCaps2;
    uint32 Reserved[2];
};


struct DDSURFACEDESC2
{
    uint32 dwSize;
    uint32 dwFlags;
    uint32 dwHeight;
    uint32 dwWidth;
    uint32 dwPitchOrLinearSize;
    uint32 dwDepth;
    uint32 dwMipMapCount;
    uint32 dwReserved1[11];
    DDPIXELFORMAT ddpfPixelFormat;
    DDSCAPS2 ddsCaps;
    uint32 dwReserved2;
};

//const float fmaxf(const float a, const float b) { return a>b ? a : b; }
//const float fminf(const float a, const float b) { return a<b ? a : b; }
//const float clamp(const float x, const float a, const float b) { return fminf(b, fmaxf(x, a)); }
//const float log2(const float x) { return log(x)/log(2.0f); }

//#define PACK_DXT5


void main()
{
	// read file
	//FILE * fptr_in = fopen("indirect.tga", "rb");
	//FILE * fptr_in = fopen("shayan.tga", "rb");
	//FILE * fptr_in = fopen("nmap_alpha.tga", "rb");
	//FILE * fptr_in = fopen("face.tga", "rb");
	FILE * fptr_in = fopen("colors_small.tga", "rb");

	if(fptr_in!=NULL)
	{
		unsigned char * pbitmap = NULL; 
		int iW, iH;
		int iNrChannels = 0;
		ReadTGA(&pbitmap, &iW, &iH, &iNrChannels, fptr_in);
		if(pbitmap!=NULL)
		{
			// the apparently magic DWORD
			uint32 uMagicDWORD;
			strncpy((char *) &uMagicDWORD, "DDS ", 4);

			// surface desc
			DDSURFACEDESC2 dds_header;
			memset(&dds_header, 0, sizeof(dds_header));
			dds_header.dwSize = 124;	// always
#ifndef PACK_DXT5
			const int iByteSizePer4x4Block = 8;		// 8 for dxt1 and 16 for dxt5
#else
			const int iByteSizePer4x4Block = 16;
#endif

			dds_header.dwFlags = DDSD_CAPS | DDSD_WIDTH | DDSD_HEIGHT | DDSD_PIXELFORMAT | DDSD_LINEARSIZE;
			dds_header.dwWidth = iW;
			dds_header.dwHeight = iH;
			dds_header.dwPitchOrLinearSize = ((dds_header.dwWidth+3)/4) * iByteSizePer4x4Block;	

			// pixel format
			dds_header.ddpfPixelFormat.dwSize = 32;
			dds_header.ddpfPixelFormat.dwFlags = DDPF_FOURCC;
#ifndef PACK_DXT5
			strncpy((char *) &dds_header.ddpfPixelFormat.dwFourCC, "DXT1", 4);
#else
			strncpy((char *) &dds_header.ddpfPixelFormat.dwFourCC, "DXT5", 4);
#endif
			
			// Caps required
			dds_header.ddsCaps.dwCaps1 = DDSCAPS_TEXTURE;

			// output texels
			const int iNrBlocks = ((dds_header.dwWidth+3)/4) * ((dds_header.dwHeight+3)/4);

			FILE * fptr_out = fopen("bla4.dds", "wb");
			fwrite(&uMagicDWORD, sizeof(uMagicDWORD), 1, fptr_out);
			fwrite(&dds_header, sizeof(dds_header), 1, fptr_out);

			// read block
			for(int y=0; y<iH; y+=4)
			{
				for(int x=0; x<iW; x+=4)
				{
					uint64 uBlockOut;
					uint32 uResultingPixels[16];
					uint8 uResultingAlphas[16];
					uint32 block_in[16];
					uint8 block_in_alpha[16];

					for(int y2=0; y2<4; y2++)
					{
						int j = y+y2; if(j>(iH-1)) j=iH-1;

						for(int x2=0; x2<4; x2++)
						{
							int i = x+x2; if(i>(iW-1)) i=iW-1;
							const int index = j*iW+i;
							assert(iNrChannels>=3);
							uint8 uR = pbitmap[index*iNrChannels+0];
							uint8 uG = pbitmap[index*iNrChannels+1];
							uint8 uB = pbitmap[index*iNrChannels+2];

#ifdef PACK_DXT5
							assert(iNrChannels==4);
							block_in_alpha[y2*4+x2] = pbitmap[index*4+3];
#endif
							block_in[y2*4+x2] = uR | (((uint32) uG)<<8) | (((uint32) uB)<<16) | 0x80000000;
						}
					}
					
#ifndef PACK_DXT5
					const int iMode = 0;	// dxt1 (fuck alpha)
#else
					ConvertBlockAlphaDXT5(block_in_alpha, 2, (int32 *) &uBlockOut, uResultingAlphas);
					fwrite(&uBlockOut, 8, 1, fptr_out);

					const int iMode = 2;	// dxt3/dxt5
#endif
					ConvertBlockHeuristic3(block_in, iMode, 2, (int32 *) &uBlockOut, uResultingPixels);
					fwrite(&uBlockOut, 8, 1, fptr_out);
				}
			}

			fclose(fptr_out);



			// clean up
			delete [] pbitmap; pbitmap = NULL;
		}

		// clean up
		fclose(fptr_in);
	}
}