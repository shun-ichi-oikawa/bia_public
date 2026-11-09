#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define N_L	11
void URAND (const int n, double *u)
{	for(int i = 0; i < n; i++)
		u[i] = ( (double)rand() + 1 ) / ( (double)RAND_MAX + 2 ) - 0.5;
	return;
}
void output(const int N_P, double *YZ)
{	char filename[20]; sprintf(filename,"%07d.ccc",N_P);
	FILE *f;
	f = fopen(filename, "wb");
	fwrite(YZ,sizeof(double),6*N_P,f);
	return;
}
void print(double *YZ, const int BASE)
{	double *yz; yz = YZ + BASE;
	if(BASE==0) printf("%16s%7s%7s%7s%7s%7s\n","X","Y","Z","U","V","W");
	printf("YZ[%3d] =",BASE);
	for(int i=0;i<6;i++) printf("%7.3f",yz[i]);
	printf("\n");
	return;
}
int main (void)
{	const int N_P = ( 2 * N_L + 1 )*( 2 * N_L + 1 )*( 2 * N_L + 1 ) + 1;
	const int NEQ = 6 * N_P;
	printf("N_P = %d, NEQ = %d\n",N_P,NEQ);
	double YZ[NEQ];
	int BASE = 0;
	for(int i=0; i<6; i++) YZ[i]=0;
	print(YZ,BASE);
	srand( (unsigned)time(NULL) );
	for(int k = - N_L; k <= N_L; k++)
	{	//int K = k + N_L;
		for(int j = - N_L; j <= N_L; j++)
		{	int J = J + N_L;
			for(int i = - N_L; i <= N_L; i++)
			{	//int I = i + N_L;
				BASE += 6;
				URAND(6,YZ+BASE);
				print(YZ,BASE);
			}
		}
	}
	output(N_P,YZ);
	return 0;
}
