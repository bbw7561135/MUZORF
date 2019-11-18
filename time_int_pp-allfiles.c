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
#include <unistd.h>

#define ASZ 150
#define FREQ_GRID 150

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
   int t_current = 0, t_previous = 0;
   double t_diff = 0.0;

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
   
   t_current = t_previous = i = t_start;
   while(1)
     {
	sprintf(filename, "%s%d%s", file_prefix, i, file_suffix);
	if ( access(filename, R_OK) != -1 )
	  {
	     t_previous = i;
	   
	     for(j = t_previous+1; j <= t_end; j++) 
	       {
		  sprintf(filename, "%s%d%s", file_prefix, j, file_suffix);
		  if( access(filename, R_OK) != -1 ) 
		    {
		       t_current = j;
		       
		       printf("\nPrevious = %d, Current = %d", t_previous, t_current);
		       
		       // we now have i and (i-1) elements identified. now we calculate the time average
		       // 
		       t_diff = (double) (t_current - t_previous);
		       sumdiff += t_diff;
		       
		       sprintf(filename, "%s%d%s", file_prefix, t_current, file_suffix);
		       
		       fpin = fopen(filename, "r");
		       if (fpin != NULL) {
			  
			  for (k = 0; k < FREQ_GRID; k++) {
			     
			     fscanf(fpin, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &freq[k], &nufnu[k],
				    &nufnusyn[k], &nufnussc[k], &nufnuecd[k], &nufnublr[k], 
				    &nufnudt[k], &nufnuup[k], &nufnudo[k]); 
			     
			     prod[k] = nufnu[k] * t_diff;
			     sum[k] += prod[k];
			     
			     prodsyn[k] = nufnusyn[k] * t_diff;
			     sumsyn[k] += prodsyn[k];
			     
			     prodssc[k] = nufnussc[k] * t_diff;
			     sumssc[k] += prodssc[k];
			     
			     prodecd[k] = nufnuecd[k] * t_diff;
			     sumecd[k] += prodecd[k];
			     
			     prodblr[k] = nufnublr[k] * t_diff;
			     sumblr[k] += prodblr[k];
			     
			     proddt[k] = nufnudt[k] * t_diff;
			     sumdt[k] += proddt[k];
			     
			     produp[k] = nufnuup[k] * t_diff;
			     sumup[k] += produp[k];
			     
			     proddo[k] = nufnudo[k] * t_diff;
			     sumdo[k] += proddo[k];
			  }

			  fclose(fpin);
			  
		       } else {
			  
			  printf("\n\nUnable to open file (%s) even though it exists. Please check the permissions on the file.\n\n", filename);
			  exit(-1);
		       }
		       
		       t_previous = t_current;
		       i = t_current;
		       
		    } // end IF file j exists
		  		
	       } // end for j=t_previous+1, looking the current file
	     
	     if ( j > t_end) 
	       {
		  break;
	       }
	   
	  } else {
	     
	     // file previous doesn't exist. go to next i
	     // 
	     i++;
	  }
	
	if (i > t_end) 
	  {
	     break;
	  }
	
	
     } // end while(1)
   
   
  printf("\n\nsumdiff = %e\n", sumdiff);

  fpout = fopen(outfile, "w");
  if (!fpout)
    {
      printf("\nCould not open Int file. \n");
      return -1;
    }
     
   
   for(j = 0; j < FREQ_GRID; j++) 
     {
	
	avgsum[j] = sum[j] / sumdiff;
	
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

  fclose(fpout);
  
  printf("\nEnd of Program\n");
  
  return 0;
  
}
