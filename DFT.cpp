#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
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

int main(int argc, char** argv){
	int numbers;

	complexS a;
	string line;
	ifstream myfile (argv[1]);
	if (myfile.is_open())
	{
		getline(myfile,line);
		istringstream iss(line, istringstream::in);
		iss >> numbers;//this number will always be 2.
		iss >> numbers;//this number represents the height.
		a.len = numbers;
		a.real = new float[a.len];
		a.imag = new float[a.len];
		for(int i = 0; i < a.len; i++){
			getline(myfile,line);
			istringstream iss(line, istringstream::in);
			iss >> numbers;
			a.real[i] = numbers;
			iss >> numbers;
			a.imag[i] = numbers;
		}
		myfile.close();
	}
	else cout << "Unable to open file";

	complexS b;
	string dir = argv[2];
	if(dir.compare("0") == 0){
		DFT(&a,0,&b);
	}else{
		DFT(&a,1,&b);
	}
	ofstream outfile;
	outfile.open (argv[3]);
	outfile<<2<<" "<<a.len<<endl;
	for(int i = 0; i < a.len; i++){
		outfile <<b.real[i]<<" "<<b.imag[i]<<endl;
	}
	outfile.close();

return 0;

}
