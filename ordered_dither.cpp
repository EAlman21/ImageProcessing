// ========================================================
// ordered_dither.cpp - dithering, halftoning, and neighborhood operations program.
//
// Written by: Christyn Vasquez
// ========================================================

#include "IP.h"

using namespace std;

void ordered_dither(imageP, int, double, imageP);

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// main:
//
// Main routine for quantization of image.
//
int
main(int argc, char** argv)
{
	int 	m, gamma;
	imageP 	inImage, outImage;
	
	// Checks image version
	if (argc != 5) {
        cerr << "Usage: ordered_dither imageP int double imageP\n";
        exit(1);
    }
	
	// Read input image and reserve space for output
	inImage = IP_readImage(argv[1]);
	outImage = NEWIMAGE;

	// Read quantization levels
	m = atoi(argv[2]);

	// Read gamma number
	gamma = atof(argv[3]);
	
	// Quantization image and save result in file
	ordered_dither(inImage, m, gamma, outImage);
	IP_saveImage(outImage, argv[4]);

	// Free up image structures/memory
	IP_freeImage(inImage);
	IP_freeImage(outImage);

	return 1;
}



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ordered_dither:
//
// 
//
void
ordered_dither(imageP inImage, int m, double gamma, imageP outImage){

	int		i, j, total, result, x, y;
	int 	jitter, bias, scale;
	int		width, height;
	uchar	*in, *out, lut[MXGRAY];

	scale = MXGRAY / pow(m,2);
	
	// Total number of pixels in image
	total = inImage->width * inImage->height;
		
	// Init out image dimensions and buffer
	outImage->width  = inImage->width;
	outImage->height = inImage->height;
	outImage->image  = (uchar *) malloc(total);
	
	if(outImage->image == NULL) {
		cerr << "ordered_dither: Insufficient memory\n";
		exit(1);
	}
	
	in  = inImage->image;	// input  image buffer
	out = outImage->image;	// output image buffer
		
	// Init gamma correction lut
	for(i = 0; i < MXGRAY; i++){
		result = pow((double)(i/MaxGray), 1.0/gamma)*MXGRAY;
		in[i] = (uchar)result;
	}

	// Init lookup table
	for(i = 0; i < MXGRAY; i++) 
			lut[i] = i/scale;
	
	// Visit all input pixels and apply lut to threshold
	for(i = 0; i < total; i++) out[i] = lut[in[i]];	
	
	// Init 3x3 dither matrix
	int D3by3[9]  = { 6, 8, 4, 
                    1, 0, 3, 
                    5, 2, 7 };

	// Init 4x4 dither matrix
	int D4by4[16] = { 0 , 8 , 2 , 10,
                    12, 4 , 14, 6 ,
                    3 , 11, 1 , 9 , 
                    15, 7 , 13, 5 };
  
	// Init 8by8 dither matrix
	int D8by8[64] = { 0 , 48, 12, 60, 3 , 51, 15, 63,
                    32, 16, 44, 28, 35, 19, 47, 31,
                    8 , 56, 4 , 52, 11, 59, 7 , 55,
                    40, 24, 36, 20, 43, 27, 39, 23,
                    2 , 50, 14, 62, 1 , 49, 13, 61,
                    34, 18, 46, 30, 33, 17, 45, 29,
                    10, 58, 6 , 54, 9 , 57, 5 , 53,
                    42, 26, 38, 22, 41, 25, 37, 21 };
	
	width = outImage->width;
	height = outImage->height;
	
    for(x = 0; x < width; x++) {
        for(y = 0; y < height; y++) {
            i = x % m;
            j = y % m;

            if(m == 3)      
                out[y*width+x] = (out[y*width+x] > D3by3[j*m+i]) ? 255 : 0;
            else if(m == 4) 
                out[y*width+x] = (out[y*width+x] > D4by4[j*m+i]) ? 255 : 0;
            else if(m == 8) 
                out[y*width+x] = (out[y*width+x] > D8by8[j*m+i]) ? 255 : 0;
			
        }
    }
}
