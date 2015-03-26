//error_diffusion.cpp
#include "IP.h"

using namespace std;

// function prototype
void error_diffusion(imageP, int, int, double, imageP);

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// main:
int
main(int argc, char** argv)
{
	double gamma;
	int mtd, serpentine;
	imageP	I1, I2;

	// error checking: proper usage
	if(argc != 6) {
		cerr << "Usage: thr infile thr_low thr_high outfile\n";
		exit(1);
	}

	// read input image (I1) and reserve space for output (I2)
	I1 = IP_readImage(argv[1]);
	I2 = NEWIMAGE;

	//mtd
	mtd =atoi(argv[2]);

	//serpentine
	serpentine =atoi(argv[3]);

	// read gamma
	gamma =atof(argv[4]);

	// threshold image and save result in file
	error_diffusion(I1, mtd, serpentine, gamma, I2);
	IP_saveImage(I2, argv[5]);

	// free up image structures/memory
	IP_freeImage(I1);
	IP_freeImage(I2);

	return 1;
}


void
error_diffusion(imageP I1, int mtd, int serpentine, double gamma, imageP I2)
{
	int	 i, total;
	uchar *in, *out, lut[256], newlut[256];

	// total number of pixels in image
	total = I1->width * I1->height;

	// init I2 dimensions and buffer
	I2->width  = I1->width;
	I2->height = I1->height;
	I2->image  = (uchar *) malloc(total);
	if(I2->image == NULL) {
		cerr << "thr: Insufficient memory\n";
		exit(1);
	}

	
	// init lookup table
	for(i=0; i<256 ; i++){
	lut[i] = i;
	//cout<<int(lut[i])<<'\n';	
	}


	//gamma correction
	gamma = 1/gamma;
	for(i = 0; i< 256; i++){
	lut[i] = (unsigned char) (255.0 * pow((i/255.0), gamma));	
	//cout << lut[i]<<'\n';
	}
    
	//floyd steinberg error diffusion fractions
	float d1 = 7.0/16.0;
	float d2 = 3.0/16.0;
	float d3 = 5.0/16.0;
	float d4 = 1.0/16.0;

	//Jarvis-Judice-Ninke fractions
	float j7 = 7.0/48.0;
	float j5 = 5.0/48.0;
	float j3 = 3.0/48.0;
	float j1 = 1.0/48.0;

	int h = I1->height;
	int w = I1->width;

	int e;
	int thr = MXGRAY /2; // init threshold value

	if(mtd == 0)
	{
	if(serpentine ==0 )
	{
	//Apply Floyd Steinberg Algorithm
		for(int y=0; y<h; y++){ // visit all input rows
		for(int x=0; x<w; x++) { // visit all input cols
		if(lut[in[(x+y)*w]]<thr){
			e = lut[in[(x+y)*w]];
		}
		else{
			e = ABS(lut[in[(x+y)*w]] - 255);
		}
		
		lut[in[(x+y)*w+1]] += int(e*d1); // add error to E nbr
		lut[in[(x+(y+1))*w-1]] += int(e*d2); // add error to SW nbr
		lut[in[(x+(y+1))*w]] += int(e*d3); // add error to S nbr
		lut[in[(x+(y+1))*w+1]] += int(e*d4); // add error to SE nbr
		
		}
	}
	}
	else if (serpentine ==1)
	{
		//Apply Floyd Steinberg Algorithm with serpentine
		for(int y=0; y<h; y+=2){ // visit all input rows
		for(int x=0; x<w; x++) { // visit all input cols
		if(lut[in[(x+y)*w]]<thr){
			e = lut[in[(x+y)*w]];
		}
		else{
			e = ABS(lut[in[(x+y)*w]] - 255);
		}
		
		lut[in[(x+y)*w+1]] += int(e*d1); // add error to E nbr
		lut[in[(x+(y+1))*w-1]] += int(e*d2); // add error to SW nbr
		lut[in[(x+(y+1))*w]] += int(e*d3); // add error to S nbr
		lut[in[(x+(y+1))*w+1]] += int(e*d4); // add error to SE nbr
		
		}
	}
		for(int y=1; y<h; y+=2){ // visit all input rows
		for(int x=w; x>0; x--) { // visit all input cols
		if(lut[in[(x+y)*w]]<thr){
			e = lut[in[(x+y)*w]];
		}
		else{
			e = ABS(lut[in[(x+y)*w]] - 255);
		}
		
		lut[in[(x+y)*w+1]] += int(e*d1); // add error to E nbr
		lut[in[(x+(y+1))*w-1]] += int(e*d2); // add error to SW nbr
		lut[in[(x+(y+1))*w]] += int(e*d3); // add error to S nbr
		lut[in[(x+(y+1))*w+1]] += int(e*d4); // add error to SE nbr
		
		}
	}
	}
	}
	else if(mtd ==1)
	{
		if(serpentine ==0 )
	{
	//Apply Jarvis-Judice-Ninke Algorithm
		for(int y=0; y<h; y++){ // visit all input rows
		for(int x=0; x<w; x++) { // visit all input cols
		if(lut[in[(x+y)*w]]<thr){
			e = lut[in[(x+y)*w]];
		}
		else{
			e = ABS(lut[in[(x+y)*w]] - 255);
		}
		
		lut[in[(x+y)*w+1]] += int(e*j7); // add error to E nbr
		lut[in[(x+y)*w+2]] += int(e*j5); // add error to EE nbr
		lut[in[(x+(y+1))*w-2]] += int(e*j3); // add error to SWW nbr
		lut[in[(x+(y+1))*w-1]] += int(e*j5); // add error to SW nbr
		lut[in[(x+(y+1))*w]] += int(e*j7); // add error to S nbr
		lut[in[(x+(y+1))*w+1]] += int(e*j5); // add error to SE nbr
		lut[in[(x+(y+1))*w+1]] += int(e*j3); // add error to SEE nbr
		lut[in[(x+(y+2))*w-2]] += int(e*j1); // add error to SSWW nbr
		lut[in[(x+(y+2))*w-1]] += int(e*j3); // add error to SSW nbr
		lut[in[(x+(y+2))*w]] += int(e*j5); // add error to SS nbr
		lut[in[(x+(y+2))*w+1]] += int(e*j3); // add error to SSE nbr
		lut[in[(x+(y+2))*w+2]] += int(e*j1); // add error to SSEE nbr
		}
	}
	}
	else if (serpentine ==1)
	{
		//Apply Jarvis-Judice-Ninke Algorithm with serpentine
		for(int y=0; y<h; y+=2){ // visit all input rows
		for(int x=w; x>0; x--) { // visit all input cols
		if(lut[in[(x+y)*w]]<thr){
			e = lut[in[(x+y)*w]];
		}
		else{
			e = ABS(lut[in[(x+y)*w]] - 255);
		}
		
		lut[in[(x+y)*w+1]] += int(e*j7); // add error to E nbr
		lut[in[(x+y)*w+2]] += int(e*j5); // add error to EE nbr
		lut[in[(x+(y+1))*w-2]] += int(e*j3); // add error to SWW nbr
		lut[in[(x+(y+1))*w-1]] += int(e*j5); // add error to SW nbr
		lut[in[(x+(y+1))*w]] += int(e*j7); // add error to S nbr
		lut[in[(x+(y+1))*w+1]] += int(e*j5); // add error to SE nbr
		lut[in[(x+(y+1))*w+1]] += int(e*j3); // add error to SEE nbr
		lut[in[(x+(y+2))*w-2]] += int(e*j1); // add error to SSWW nbr
		lut[in[(x+(y+2))*w-1]] += int(e*j3); // add error to SSW nbr
		lut[in[(x+(y+2))*w]] += int(e*j5); // add error to SS nbr
		lut[in[(x+(y+2))*w+1]] += int(e*j3); // add error to SSE nbr
		lut[in[(x+(y+2))*w+2]] += int(e*j1); // add error to SSEE nbr
		}
	}

		for(int y=1; y<h; y+=2){ // visit all input rows
		for(int x=w; x>0; x--) { // visit all input cols
		if(lut[in[(x+y)*w]]<thr){
			e = lut[in[(x+y)*w]];
		}
		else{
			e = ABS(lut[in[(x+y)*w]] - 255);
		}
		
		lut[in[(x+y)*w+1]] += int(e*j7); // add error to E nbr
		lut[in[(x+y)*w+2]] += int(e*j5); // add error to EE nbr
		lut[in[(x+(y+1))*w-2]] += int(e*j3); // add error to SWW nbr
		lut[in[(x+(y+1))*w-1]] += int(e*j5); // add error to SW nbr
		lut[in[(x+(y+1))*w]] += int(e*j7); // add error to S nbr
		lut[in[(x+(y+1))*w+1]] += int(e*j5); // add error to SE nbr
		lut[in[(x+(y+1))*w+1]] += int(e*j3); // add error to SEE nbr
		lut[in[(x+(y+2))*w-2]] += int(e*j1); // add error to SSWW nbr
		lut[in[(x+(y+2))*w-1]] += int(e*j3); // add error to SSW nbr
		lut[in[(x+(y+2))*w]] += int(e*j5); // add error to SS nbr
		lut[in[(x+(y+2))*w+1]] += int(e*j3); // add error to SSE nbr
		lut[in[(x+(y+2))*w+2]] += int(e*j1); // add error to SSEE nbr
		}
	}

	}
	}

	// visit all input pixels and apply lut to threshold
	in  = I1->image;	// input  image buffer
	out = I2->image;	// output image buffer
	for(i=0; i<total; i++) out[i] = lut[ in[i] ];
}