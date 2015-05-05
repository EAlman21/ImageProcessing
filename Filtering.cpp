#include "IP.h"
#include <math.h>

using namespace std;
double PI = 3.14159265358979323846;
typedef struct{
	int len;
	float *real;
	float *imag;
} complexS, *complexP;

int cel(int x,int y){
	//return  (x % y == 0) ? x / y + 1 : x / y);
	if (x%y == 0) return x/y;
	else return x/y+1;
}

bool poweroftwo(int n){
	if (n == 1) return true;
	if(n%2 == 0) return poweroftwo(n/2);
	return false;
}

void DFT(complexP in,bool dir, complexP out){
	out->len = in->len;
	out->real = new float[out->len];
	out->imag = new float[out->len];
	float real, imag,c,s;
	int d;
	if(dir == 0) d = 1;
	else d = -1;
		for(int u=0; u < out->len; u++){
			real = imag = 0;
			for(int x = 0; x < out->len; x++){
				c = cos(d*-2 * PI * u *x/out->len);
				s = sin(d*-2 * PI * u *x/out->len);
				real += in->real[x]*c - in->imag[x]*s;
				imag += in->real[x]*s + in->imag[x]*c;

			}

			if(dir){
				out->real[u] = real/in->len;
				out->imag[u] = imag/in->len;
			}else{
				out->real[u] = real;
				out->imag[u] = imag;				
			}
		}
}

void FFT(complexP in,bool dir, complexP out){
        if(!poweroftwo(in->len)){
		DFT(in,dir,out);
		return;
	}
	out->len = in->len;
	out->real = new float[out->len];
	out->imag = new float[out->len];
	int d;
	if(dir == 0) d = 1;
	else d = -1;
		if(in->len == 1) //bottom of tree, do nothing.
		{

			for(int i = 0; i<in->len; i++){
				out->real[i] = in->real[i];
				out->imag[i] = in->imag[i];
			}

		}else{
			complexP even, odd;
			even = new complexS;
			odd = new complexS;
			even->len = cel(in->len,2); // even will have pairs when there is an odd number of elements.
			even->real = new float[even->len];
			even->imag = new float[even->len];		
			odd->len = in->len/2;
			odd->real = new float[odd->len];
			odd->imag = new float[odd->len];
			//build even
			for(int i = 0; i < cel(in->len,2); i++){
				even->real[i] = in->real[i*2];
				even->imag[i] = in->imag[i*2];
			}
			//build odd
			for(int i = 0; i < in->len/2; i++){
				odd->real[i] = in->real[i*2 + 1];
				odd->imag[i] = in->imag[i*2 + 1];
			}
			complexS reven,rodd; 
			FFT(even,dir,&reven);
			FFT(odd,dir,&rodd);
			for(int i = 0; i < rodd.len; i++){
				float eIMAG;
				float eREAL;
				eIMAG = sin(d*-2*PI*(i)/(out->len));
				eREAL = cos(d* -2*PI*(i)/(out->len));

				out->real[i] = eREAL*rodd.real[i] - eIMAG* rodd.imag[i] + reven.real[i];
				out->imag[i] = eREAL*rodd.imag[i] + eIMAG*rodd.real[i] + reven.imag[i];

			}
			for(int i=rodd.len; i < out->len; i++){
				float eIMAG;
				float eREAL;
				eIMAG = sin(d*-2*PI*(rodd.len-i)/(out->len));
				eREAL = cos(d*-2*PI*(rodd.len-i)/(out->len));
				out->real[i] = reven.real[i - rodd.len] - eREAL*rodd.real[i - rodd.len] - eIMAG* rodd.imag[i - rodd.len];
				out->imag[i] = reven.imag[i- rodd.len] - eREAL*rodd.imag[i - rodd.len] + eIMAG*rodd.real[i - rodd.len];
			}
			delete(even);//Not
			delete(odd);//Sure
			delete(reven.imag);//If
			delete(reven.real);//This
			delete(rodd.imag);//Is
			delete(rodd.real);//Necessary.
			
		if(dir){
			for(int i = 0; i < out->len; i++){
				out->real[i] = out->real[i]/2;
				out->imag[i] = out->imag[i]/2;
			}
		}
	}


}

complexS CompfromImag(imageP image){
	int total = image->width * image->height;	
	complexS temp;
	temp.len = total;
	temp.real = new float[total];
	temp.imag = new float[total];
	
	for(int i =0; i <total ; ++i){
		temp.real[i] = image->image[i];
		temp.imag[i] = 0;		
	}
	
	return temp;
}

complexP TwoDFFT(complexP Compimage, int dir, int width, int height){

	
	complexS Column;
	Column.len = height;
	Column.real = new float[height];
	Column.imag = new float[height];
	

	complexS DFTColumn;
	DFTColumn.len = height;
	DFTColumn.real = new float[height];
	DFTColumn.imag = new float[height];

	complexS DFTRow;
	DFTRow.len  = width;
	DFTRow.imag = new float[height]; 
	DFTRow.real = new float[height];

	complexS Row;
	Row.len  = width;
	Row.imag = new float[height]; 
	Row.real = new float[height];
	

	complexP temp;
	temp = new complexS;
	temp->len = width * height;
	temp->real = new float[width * height];
	temp->imag = new float[width * height];


	for(int i = 0; i < width; i++){
		//construct a list of conplex numbers from column i
		for(int j = 0; j<height; j++){
			Column.real[j] = Compimage->real[i + width*j];
			Column.imag[j] = Compimage->imag[i + width*j];
		}
		//now we get FFT/DFT

		FFT(&Column,dir,&DFTColumn);
		//new move DFTED to outImage

		for(int j = 0; j<height; j++){
			temp->real[i + width*j] = DFTColumn.real[j];
			temp->imag[i + width*j] = DFTColumn.imag[j];
			//temp->image[i + inImage->width*j] = sqrt(DFTED.real[j]*DFTED.real[j] + DFTED.imag[j]*DFTED.imag[j]);
		}
		
	}

	// now we do all rows of temp.
	for(int i = 0; i < height; i++){
		//construct list of complex number from row i
		for(int j = 0; j<width; j++){
			Row.real[j] = temp->real[i*width + j];
			Row.imag[j] = temp->imag[i*width + j];
		}
		
		FFT(&Row,dir,&DFTRow);

	
		for(int j = 0; j<width; j++){
			temp->real[i*width +j] = DFTRow.real[j];
			temp->imag[i*width +j] = DFTRow.imag[j];
		}
	}		


	return temp;

}

void MagfromComplex(complexP Compimage, imageP image,bool scale = true){
	int max = 255;
	int min = 0;
	int current;
	
	image->image = (uchar *) malloc(Compimage->len);

	//first find max and min.
	for(int i =0; i < Compimage->len; i++){
		current = sqrt(Compimage->real[i] * Compimage->real[i] + Compimage->imag[i] *Compimage->imag[i]);
		if (max<current) max = current;
		if (min > current) min = current;
	}

	// logarithmic scaling.

	double c = 255/log(1+(max + min));
		
	// add values to the image;
	for(int i=0; i <Compimage->len; i++){
		if(scale){
		image->image[i] = c * log(1 + sqrt(Compimage->real[i] * Compimage->real[i] + Compimage->imag[i] *Compimage->imag[i]) + min );
		}else{
		image->image[i] = sqrt(Compimage->real[i] * Compimage->real[i] + Compimage->imag[i] *Compimage->imag[i]);
		}
	}
}

void PhafromComplex(complexP Compimage, imageP image, bool scale = true){
	int max = 255;
	int min = 0;
	int current;
	//first find max and min.
	for(int i = 0; i < Compimage->len; i++){
		current = atan2(Compimage->imag[i],Compimage->real[i]);
		if (max<current) max = current;
		if (min > current) min = current;
	}

	// add values to the image;
	for(int i=0; i <Compimage->len; i++){
		if(scale){
		image->image[i] = (atan2(Compimage->imag[i],Compimage->real[i]) + min) * 255 / (max+min);
		}else{
		image->image[i] = atan2(Compimage->imag[i],Compimage->real[i]);
		}	
	}

}

void Filtering(imageP im1, imageP im2, imageP outm){
	
	imageP gawbage;
	gawbage = NEWIMAGE;
	gawbage->width = im1->width;
	gawbage->height = im1->height;
	gawbage->image = (uchar *) malloc(im1->width * im1->height);
	outm->width = im1->width;
	outm->height = im1->height;
	outm->image = (uchar *) malloc(im1->width * im1->height);

	for(int i =0; i<im1->width * im1->height; i++){
		gawbage->image[i] = 0;		
	}

	gawbage->image[125*im1->width + 125] = 255;
	gawbage->image[125*im1->width + 126] = 255;
	gawbage->image[125*im1->width + 124] = 255;
	gawbage->image[124*im1->width + 125] = 255;
	gawbage->image[126*im1->width + 125] = 255;
	gawbage->image[126*im1->width + 126] = 255;
	gawbage->image[126*im1->width + 124] = 255;
	gawbage->image[124*im1->width + 126] = 255;
	gawbage->image[124*im1->width + 124] = 255;

	complexS s1 = CompfromImag(im1);
	complexS s2 = CompfromImag(gawbage);

	complexP P1 = TwoDFFT(&s1,0,im1->width,im1->height);
	complexP P2 = TwoDFFT(&s2,0,im1->width,im1->height);
	complexS Pout;
	Pout.len = P1->len;
	Pout.real = new float[P1->len];
	Pout.imag = new float[P1->len];
	
	for(int i = 0; i<im1->width*im1->height; ++i){
		Pout.real[i] = P1->real[i] * P2->real[i] - P1->imag[i]*P2->imag[i];
		Pout.imag[i] = P1->real[i] *P2->imag[i] + P1->imag[i] * P2->real[i];
	}
	
	complexP PFFT = TwoDFFT(&Pout,1,im1->width,im1->height);
	
	MagfromComplex(PFFT,outm);

}



int
main(int argc, char** argv)
{
	int 	n;
	imageP 	inImage,inImage2, outImage, outImage2;
	
	// Checks image version
	//if (argc != 4) {
        //cerr << "Usage: qntz infile n outfile\n";
        //exit(1);
    //}
	
	// Read input image and reserve space for output
	inImage = IP_readImage(argv[1]);
	inImage2 = IP_readImage(argv[2]);
	outImage = NEWIMAGE;
	outImage2 = NEWIMAGE;

	//spec(inImage,outImage,outImage2);
	//PhaseSwap(inImage,inImage2,outImage,outImage2);	
	Filtering(inImage,inImage2,outImage);	
	cout<<"so far so good"<<endl;
	//IP_saveImage(outImage2, argv[4]);	
	IP_saveImage(outImage, argv[3]);
	

	// Free up image structures/memory
	IP_freeImage(inImage);
	IP_freeImage(outImage);
	IP_freeImage(outImage2);
	return 1;
}
