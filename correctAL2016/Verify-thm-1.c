/****************************************************************/
/* This file is part of an auxiliary package to the paper:  	*/
/* Rotational Cryptanalysis in the Presence of Constants 	*/
/* by Tomer Ashur and Yunwen Liu, KU Leuven.	 		*/
/* It is to empirically verify the correctness of  		*/
/* Theorem 1, over all 2^{32} possible inputs			*/
/* This program is provided without any warranty and you 	*/
/* are using it at your responsibility.				*/
/* If you ever use this program in an official publication,	*/
/* please cite:     						*/
/* Tomer Ashur, Yunwen Liu: 					*/
/* Rotational Cryptanalysis in the Presence of Constants. 	*/
/* FSE 2017: (to be published)					*/
/* For any questions please contact tashur@esat.kuleuven.be	*/
/****************************************************************/

#include <stdint.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#define u64 unsigned long long int

#define LCS(in,amount) ((((in&0x0000ffff) << (amount)) | ((in&0x0000ffff) >> (16 - amount))) & 0x0000ffff)
#define RCS(in,amount) ((((in&0x0000ffff) >> (amount)) | ((in&0x0000ffff) << (16 - amount))) & 0x0000ffff)

int main(int argc, char **argv)
{

	u64 counter0 = 0;
	
	int a1 = 0;
	int b1 = 0;
	int Delta1 = 0;

	int a2, b2, Delta2;

	#ifdef MISMATCH
		printf("Transition with incorrect fractional part (not .0 or .415)\n");
		a2 = LCS(0x2556, 7);
		b2 = 0x0b21;
		Delta2 = 0x2e76;
	#else
		printf("Transition with incorrect probability (factor x2, example in the paper)\n");
		a2 = LCS(0xe013, 7);
		b2 = 0x13a1;
		Delta2 = 0x0c4d;
	#endif
	
	int counter = 0;
	
	#pragma omp parallel for default(none) collapse(2) reduction (+:counter0) firstprivate(a1,a2,b1,b2,Delta1,Delta2)
	for (int i=0;i<= 0xffff;i++)
		for (int j=0;j<=0xffff;j++)
		{
			int x = i;
			int y = j;
			int xp = LCS(x,1);
			int yp = LCS(y,1);
			if (LCS(((RCS((x ^ a1),7) + (y ^ b1)) & 0xffff) ^ Delta1,1) == (((RCS((xp ^ a2),7) + (yp ^ b2)) & 0xffff) ^ Delta2))
				counter0++;
		}

	printf ("counter0: %lld %lf\n",counter0, (double)log((double)counter0)/(double)log((double)2));
	/*printf ("counter1: %lld\n",counter1);
	printf ("counter2: %lld\n",counter2);
	printf ("counter3: %lld\n",counter3);
	printf ("counter4: %lld\n",counter4);
	printf ("counter5: %lld\n",counter5);*/
	//printf ("counter6: %lld\n",counter6);
	//printf ("counter7: %lld\n",counter7);
	printf("\n");



}
