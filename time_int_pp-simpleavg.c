/*
  This program is to obtain a time weighted average of the nuFnu values using 
  the files obtained from the post_process.c code. The input file for this 
  program must contain the file name prefix, starting, ending, and one-step-before 
  starting time, file name suffix, integration time window, and resulting filename.

*/

#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define ASZ 50
#define FREQ_GRID 50

int main() {

  FILE *fpin, *fpout, *fpvar;
  char *varfile = "variables.txt";
  int j = 0, i = 0, k = 0, t_start = 0, t_end = 0, window = 0, t_begin = 0;
  double sumdiff = 0.0;
  double freq[ASZ], nufnu[ASZ], sum[ASZ], avgsum[ASZ], prod[ASZ];
  double diff[15000];
  double nufnusyn[ASZ], nufnussc[ASZ], nufnuecd[ASZ], nufnublr[ASZ]; 
  double nufnudt[ASZ], prodsyn[ASZ], prodssc[ASZ], prodecd[ASZ], prodblr[ASZ]; 
  double proddt[ASZ], sumsyn[ASZ], sumssc[ASZ], sumecd[ASZ], sumblr[ASZ]; 
  double sumdt[ASZ], avgsumsyn[ASZ], avgsumssc[ASZ], avgsumecd[ASZ]; 
  double avgsumblr[ASZ], avgsumdt[ASZ];
  double nufnuup[ASZ], nufnudo[ASZ], produp[ASZ], proddo[ASZ];
  double sumup[ASZ], sumdo[ASZ], avgsumup[ASZ], avgsumdo[ASZ]; 

  char* filename = (char*) malloc(sizeof(char) * 1024);
  char* file_prefix = (char*) malloc(sizeof(char) * 1024);
  char* file_suffix = (char*) malloc(sizeof(char) * 1024);
  char* outfile = (char*) malloc(sizeof(char) * 1024);

  // Initialize all the arrays.
  //
  for (i = 0; i < FREQ_GRID; i++) {
    
    freq[i] = nufnu[i] = sum[i] = avgsum[i] = prod[i] = 0.0;
    nufnusyn[i] = nufnussc[i] = nufnuecd[i] = nufnublr[i] = nufnudt[i] = 0.0;
    nufnuup[i] = nufnudo[i] = 0.0;
    
  }

  // Read input file and get parameters
  //
  fpvar = fopen(varfile, "r");
  if(!fpvar)
    {
      printf ("\nCould not open Var file.\n");
      return -1;
    }

  fscanf(fpvar, "%s %d %d %d %s %d %s\n", file_prefix, &t_start, &t_end, &t_begin, file_suffix, &window, outfile);
  printf("%s %d %d %d %s %d %s\n", file_prefix, t_start, t_end, t_begin, file_suffix, window, outfile);

  fclose(fpvar);

  // Sum the time intervals here.
  //
  k = 0;

  for (j = t_start; j <= t_end; j += window) {

    sprintf(filename, "%s%d%s", file_prefix, j, file_suffix);

    fpin = fopen(filename, "r");
    if (fpin != NULL) {

      printf("File exists %s\n", filename);

      diff[k] = ((double) (j)) - ((double) (t_begin));
      sumdiff += diff[k];
      printf("j = %d, t_begin = %d, diff[%d] = %e, sumdiff = %e\n", j, t_begin, k, diff[k], sumdiff);

      t_begin = j;
      k++;

    }
    
    fclose(fpin);
    
  }

  printf("sumdiff = %e, k = %d\n", sumdiff, k);

  fpout = fopen(outfile, "w");
  if (!fpout)
    {
      printf("\nCould not open Int file. \n");
      return -1;
    }
     
  // Calculate the time-weighted average of the nuFnu values here.
  //
  k = 0;

  for (i = t_start; i <= t_end; i += window) {

    sprintf(filename, "%s%d%s", file_prefix, i, file_suffix);
    
    fpin = fopen(filename, "r");
    if (fpin != NULL) {
      
      for (j = 0; j < FREQ_GRID; j++) {
	
	fscanf(fpin, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &freq[j], &nufnu[j], 
	       &nufnusyn[j], &nufnussc[j], &nufnuecd[j], &nufnublr[j], 
	       &nufnudt[j], &nufnuup[j], &nufnudo[j]); 
	  	  		
	//prod[j] = nufnu[j] * diff[k];
	prod[j] = nufnu[j];
	sum[j] += prod[j];

        prodsyn[j] = nufnusyn[j] * diff[k];
        sumsyn[j] += prodsyn[j];

        prodssc[j] = nufnussc[j] * diff[k];
        sumssc[j] += prodssc[j];

        prodecd[j] = nufnuecd[j] * diff[k];
        sumecd[j] += prodecd[j];

        prodblr[j] = nufnublr[j] * diff[k];
        sumblr[j] += prodblr[j];

        proddt[j] = nufnudt[j] * diff[k];
        sumdt[j] += proddt[j];

        produp[j] = nufnuup[j] * diff[k];
        sumup[j] += produp[j];

        proddo[j] = nufnudo[j] * diff[k];
        sumdo[j] += proddo[j];
	
	if (i == t_end) {
	  
	  //avgsum[j] = sum[j] / sumdiff;
	  avgsum[j] = sum[j] / k;

	  avgsumsyn[j] = sumsyn[j] / sumdiff;
	  avgsumssc[j] = sumssc[j] / sumdiff;
	  avgsumecd[j] = sumecd[j] / sumdiff;
	  avgsumblr[j] = sumblr[j] / sumdiff;
	  avgsumdt[j] = sumdt[j] / sumdiff;
	  avgsumup[j] = sumup[j] / sumdiff;
	  avgsumdo[j] = sumdo[j] / sumdiff;

	  fprintf(fpout, "%e %e %e %e %e %e %e %e %e\n", freq[j], avgsum[j], 
		  avgsumsyn[j], avgsumssc[j], avgsumecd[j], avgsumblr[j], 
		  avgsumdt[j], avgsumup[j], avgsumdo[j]); 
	  
	}
	
      }
      
      k++;

    }

    fclose(fpin);
    
  }
  
  fclose(fpout);
  
  printf("\nEnd of Program\n");
  
  return 0;
  
}
