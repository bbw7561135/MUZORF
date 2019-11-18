/*
  Post-processing of the files obtained from the shock_rad code. The original 
  nuFnu files are binned and combined into new files according to the specfied 
  value of the bin size. The new files are normalized in the last part of the 
  code according to the total number of entries in each bin. The bin size is 
  passed as a command line parameter.

  Use of command line arguments: 
  while compiling this program, use the following command:

  ./exec_file bin_size
  
  */

#include <stdio.h>
#include <dirent.h>             /* Directory information.*/
#include <stdlib.h>
#include <malloc.h>
#include <strings.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include "sortfile-mfa.h"

#define ASZ 150
#define FREQ_GRID 150                           // Number of grid points
#define PPCSIZE 6


double max(int x, int y)
{

  if (x > y) return (x);
  else return (y);
}


double min(int x, int y)
{
  if (x > y) return (y);
  else return (x);
}



int main(int argc, char *argv[]) {
   
   DIR  *dir_p, *dir_test, *dir_test_norm;
   struct dirent *dir_entry_p;
   char *filename, *file_prefix, *new_file, *old_file, *filename_part;
   char *filename_copy, *file_sim_count, *file_time_count, *window_filename;
   char *window_dir, *window_dir_norm, *outlc, *outhic1, *outal, *outhic2;
   char *outhic3, *outhic4, *outhic5;
   int match_count = 0, file_prefix_size = -1, sim_count, i, j, obs_time;
   int n = 0, window;
   double window_existing, dont_care_1, dont_care_2;
   double nuFnu_zn[ASZ], freq[ASZ], lc1, lc2, lc3, nu[PPCSIZE];
   double lc4, lc5, nu1, nu2, nu3, nu4, nu5, D_obstime;
   double alpha1, alpha2, alpha3, alpha4, alpha5;
   FILE *winfp, *fpin, *fplc, *fphic1, *fpal, *fphic2, *fphic3;
   FILE *fphic4, *fphic5;
   double nuFnusyn_zn[ASZ], nuFnussc_zn[ASZ], nuFnuecd_zn[ASZ], nuFnublr_zn[ASZ], nuFnudt_zn[ASZ];
   double syn_existing, ssc_existing, ecd_existing, blr_existing, dt_existing;
   int obstime_min = 6000000, obstime_max = 10, winobs_time, obstime_nextbinmin, obstime_nextbinmax;
   int FileStartBin, HoleStartBin, Last, Current, HoleRange, LastRatio, CurrentRatio, run_id;
   int n_nextbin = 0, last_bin = 0;
   char *winfilename_copy, *winfilename_part, *winfile_sim_count; 
   char *winfile_time_count, *nextbin_file, *nextbin_minfile; 
   FILE *nextbinfp, *nextbin_minfp;
   double nextbinnuFnu_zn[ASZ], nextbinnuFnusyn_zn[ASZ], nextbinnuFnussc_zn[ASZ], nextbinnuFnuecd_zn[ASZ]; 
   double nextbinnuFnublr_zn[ASZ], nextbinnuFnudt_zn[ASZ], nextbinmin_nuFnu_zn[ASZ]; 
   double nextbinmin_nuFnussc_zn[ASZ], nextbinmin_nuFnuecd_zn[ASZ], nextbinmin_nuFnublr_zn[ASZ]; 
   double nextbinmin_nuFnudt_zn[ASZ], nextbinmin_nuFnusyn_zn[ASZ];
   double nuFnuup_zn[ASZ], nuFnudo_zn[ASZ], up_existing, do_existing, nextbinnuFnuup_zn[ASZ]; 
   double nextbinnuFnudo_zn[ASZ], nextbinmin_nuFnuup_zn[ASZ], nextbinmin_nuFnudo_zn[ASZ];
   double nu6, nu7, nu8, nu9, nu10, lc6, lc7, lc8, lc9, lc10, alpha6, alpha7, alpha8, alpha9, alpha10;

   filename = (char*) malloc(sizeof(char) * 1024);
   file_prefix = (char*) malloc(sizeof(char) * 1024);
   file_prefix = "nuFnu_com_cyl";
   filename_part = (char*) malloc(sizeof(char) * 1024);
   filename_copy = (char*) malloc(sizeof(char) * 1024);
   file_sim_count = (char*) malloc(sizeof(char) * 1024);
   file_time_count = (char*) malloc(sizeof(char) * 1024);
   window_filename = (char*) malloc(sizeof(char) * 1024);
   old_file = (char*) malloc(sizeof(char) * 1024);
   new_file = (char*) malloc(sizeof(char) * 1024);
   window_dir = (char*) malloc(sizeof(char) * 1024);
   window_dir_norm = (char*) malloc(sizeof(char) * 1024);
   outlc = (char*) malloc(sizeof(char) * 1024);
   outal = (char*) malloc(sizeof(char) * 1024);
   outhic1 = (char*) malloc(sizeof(char) * 1024);
   outhic2 = (char*) malloc(sizeof(char) * 1024);
   outhic3 = (char*) malloc(sizeof(char) * 1024);
   outhic4 = (char*) malloc(sizeof(char) * 1024);
   outhic5 = (char*) malloc(sizeof(char) * 1024);
   winfilename_part = (char*) malloc(sizeof(char) * 1024);
   winfilename_copy = (char*) malloc(sizeof(char) * 1024);
   winfile_sim_count = (char*) malloc(sizeof(char) * 1024);
   winfile_time_count = (char*) malloc(sizeof(char) * 1024);
   nextbin_file = (char*) malloc(sizeof(char) * 1024);
   nextbin_minfile = (char*) malloc(sizeof(char) * 1024);

   // OVERWRITE DEFAULTS IF COMMAND LINE PARAMETERS PROVIDED
   //
   file_prefix_size = strlen(file_prefix);
   

   // PROCESS FILES
   //
   // Open the current directory
   // 
   dir_p = opendir(".");
   
   // Read each entry until NULL
   // 
   while( NULL != (dir_entry_p = readdir(dir_p))) {

     match_count = strspn(dir_entry_p->d_name, file_prefix);
     
     if(match_count == file_prefix_size) {
       
       filename = dir_entry_p->d_name;
       
     } else {
       
       continue;
     }
     
     sim_count = -1;
     
     for (i = 0; i < FREQ_GRID; i++) {
       
       nuFnu_zn[i] = 0.0;
       freq[i] = 0.0;
       nuFnusyn_zn[i] = nuFnussc_zn[i] = nuFnuecd_zn[i] = 0.0;
       nuFnublr_zn[i] = nuFnudt_zn[i] = 0.0;
       nuFnuup_zn[i] = nuFnudo_zn[i] = 0.0;
     }
     
     strcpy(filename_copy, filename);
     filename_part = strtok(filename_copy, "_");
     filename_part = strtok(NULL, "_");
     filename_part = strtok(NULL, "_");
     memcpy(file_sim_count, filename_part, strlen(filename_part)+1);
     file_sim_count[0] = 's';
     file_sim_count[1] = 'i';
     file_sim_count[2] = 'm';
     filename_part = strtok(NULL, ".");
     file_time_count = filename_part;
     file_time_count[0] = ' ';
     sscanf(file_time_count, "%d", &obs_time);
     
     window = atoi(argv[1]);

     i = obs_time % window;
     sprintf(window_dir, "nuFnu_%s", file_sim_count);
     sprintf(window_filename, "%s/nuFnu_obs_%s_t%d.dat", window_dir, file_sim_count, (obs_time - i));
     
     printf("%d %s %s\n", obs_time, filename, window_filename);
     
     dir_test = opendir(window_dir);
     if(dir_test == NULL) {
       
       // CREATE DIR
       //
       mkdir(window_dir, 0777);

     } else {
       
       closedir(dir_test);
     }
     
     winfp = fopen(filename, "r");
     if(winfp != NULL) {
       
       for (j = 0; ((j < FREQ_GRID) && (!feof(winfp))); j++) {	 
	 
	 fscanf(winfp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
		&freq[j], &nuFnu_zn[j], &nuFnusyn_zn[j], &nuFnussc_zn[j], 
		&nuFnuecd_zn[j], &nuFnublr_zn[j], &nuFnudt_zn[j], 
		&nuFnuup_zn[j], &nuFnudo_zn[j]);
       }
       
       fclose(winfp);
       
       winfp = fopen(window_filename, "r");
       if(winfp != NULL) {
	 
	 // File already exists so read the counter value from the first line
	 //
	 fscanf(winfp, "%d %d %d\n", &n, &obstime_min, &obstime_max);
	 
	 // File already exists so add to it
	 //
	 for (j = 0; ((j < FREQ_GRID) && (!feof(winfp))); j++) {	    
	   
           fscanf(winfp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &freq[j], 
		  &window_existing, &syn_existing, &ssc_existing, 
		  &ecd_existing, &blr_existing, &dt_existing, &up_existing, 
		  &do_existing);

           nuFnu_zn[j] = nuFnu_zn[j] + window_existing;
	   nuFnusyn_zn[j] = nuFnusyn_zn[j] + syn_existing;
	   nuFnussc_zn[j] = nuFnussc_zn[j] + ssc_existing;
	   nuFnuecd_zn[j] = nuFnuecd_zn[j] + ecd_existing;
	   nuFnublr_zn[j] = nuFnublr_zn[j] + blr_existing;
	   nuFnudt_zn[j] = nuFnudt_zn[j] + dt_existing;
	   nuFnuup_zn[j] += up_existing;
	   nuFnudo_zn[j] += do_existing;

	 }
	 
	 fclose(winfp);
	 
       }
       
       
       // NOW WRITE THE VALUES INTO THE WINDOW FILE
       //
       winfp = fopen(window_filename, "w");
       
       if(winfp == NULL) {
	 
	 printf("\n\n3 Unable to open file to write to: %s\n\n", window_filename);
	 exit(3);
       }
       
       n++;
       obstime_min = min(obstime_min, obs_time);
       obstime_max = max(obstime_max, obs_time);

       fprintf(winfp, "%d %d %d\n", n, obstime_min, obstime_max);
       
       for (j = 0; j < FREQ_GRID; j++) {	 
	 
	 fprintf(winfp, "%e %e %e %e %e %e %e %e %e\n", freq[j], nuFnu_zn[j], 
		 nuFnusyn_zn[j], nuFnussc_zn[j], nuFnuecd_zn[j], 
		 nuFnublr_zn[j], nuFnudt_zn[j], nuFnuup_zn[j], nuFnudo_zn[j]);

       }
       
       fclose(winfp);
       
     } else {
       
       printf("\n\n1 Unable to open file: %s\n\n", filename);
       printf("\n%s", strerror(errno));
       exit (1);
     }
     
     n = 0;
     obstime_min = 6000000;
     obstime_max = 10;
     
     last_bin = max(last_bin, (obs_time - i));
     
   }
   
   // Tidy up
   //
   closedir(dir_p);
   
   
   // Normalize the final values to annul the effect of binning size.
   //
   sprintf(window_dir_norm, "%s_norm", window_dir);
   
   dir_test_norm = opendir(window_dir_norm);
   if(dir_test_norm == NULL) {
     
     // CREATE DIR
     //
     mkdir(window_dir_norm, 0777);
     
   } else {
     
     closedir(dir_test_norm);
   }
   
   dir_test = opendir(window_dir);
   printf("\nDirectory opened is %s\n", window_dir);
   
   if(dir_test == NULL) {
     printf("Could not open the nuFnu_sim sub-directory.\n");
     exit (0);

   } else {
     
     while(NULL != (dir_entry_p = readdir(dir_test))) {
       
       if ((strcmp(dir_entry_p->d_name, ".") == 0) || (strcmp(dir_entry_p->d_name, "..") == 0)) continue;
       
       window_filename = dir_entry_p->d_name;
       sprintf(old_file, "%s/%s", window_dir, window_filename);

       strcpy(winfilename_copy, window_filename);
       winfilename_part = strtok(winfilename_copy, "_");
       winfilename_part = strtok(NULL, "_");
       winfilename_part = strtok(NULL, "_");
       winfile_sim_count = winfilename_part;
       winfile_sim_count[0] = ' ';
       winfile_sim_count[1] = ' ';
       winfile_sim_count[2] = ' ';
       sscanf(winfile_sim_count, "%d", &run_id);
       winfilename_part = strtok(NULL, ".");
       winfile_time_count = winfilename_part;
       winfile_time_count[0] = ' ';
       sscanf(winfile_time_count, "%d", &winobs_time);

       if ((winobs_time + window) < last_bin) {

       sprintf(nextbin_file, "%s/nuFnu_obs_%s_t%d.dat", window_dir, file_sim_count, (winobs_time + window));

       winfp = fopen(old_file, "r");
       if(winfp == NULL) {

	 printf("\n\n3 Unable to open file for reading: %s\n\n", old_file);
	 exit(3);
       }
       
       fscanf(winfp, "%d %d %d\n", &n, &obstime_min, &obstime_max);

       for (j = 0; ((j < FREQ_GRID) && (!feof(winfp))); j++) {

         fscanf(winfp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &freq[j], 
		&nuFnu_zn[j], &nuFnusyn_zn[j], &nuFnussc_zn[j], 
		&nuFnuecd_zn[j], &nuFnublr_zn[j], &nuFnudt_zn[j], 
		&nuFnuup_zn[j], &nuFnudo_zn[j]);

       }

       fclose(winfp);

       nextbinfp = fopen(nextbin_file, "r");
       if(nextbinfp == NULL) {

	 printf("\n\n3 Unable to open the next bin file for reading: %s\n\n", nextbin_file);
	 exit(3);
       }

	 fscanf(nextbinfp, "%d %d %d\n", &n_nextbin, &obstime_nextbinmin, &obstime_nextbinmax);

       for (j = 0; ((j < FREQ_GRID) && (!feof(nextbinfp))); j++) {

         fscanf(nextbinfp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &freq[j], 
		&nextbinnuFnu_zn[j], &nextbinnuFnusyn_zn[j], 
		&nextbinnuFnussc_zn[j], &nextbinnuFnuecd_zn[j], 
		&nextbinnuFnublr_zn[j], &nextbinnuFnudt_zn[j], 
		&nextbinnuFnuup_zn[j], &nextbinnuFnudo_zn[j]);

       }       

       fclose(nextbinfp);


       // Looking for a hole in the bin and if it exists then adding the appropriate
       // fraction of the first file of the next bin to the sum total values of nuFnu
       // of the current bin before dividing them by binsize and subtracting that 
       // fraction from the sum total values of nuFnu of the next bin.
       //
       HoleStartBin = floor(obstime_max / window) * window;
       FileStartBin = floor(obstime_nextbinmin / window) * window;
       
       if ((FileStartBin - HoleStartBin) == window) {
	 	 		
	 Last = (winobs_time + window) - obstime_max;
	 HoleRange = obstime_nextbinmin - obstime_max;
	 LastRatio = Last / HoleRange;

	   sprintf(nextbin_minfile, "./nuFnu_com_cyl%d_t%d.dat", run_id, obstime_nextbinmin);

	 nextbin_minfp = fopen(nextbin_minfile, "r");
	 if(nextbin_minfp == NULL) {
	   
	   printf("\n\n3 Unable to open file that exists in output folder for reading: %s\n\n", nextbin_minfile);
	   exit(3);
	 }

	 for (j = 0; ((j < FREQ_GRID) && (!feof(nextbin_minfp))); j++) {

	   fscanf(nextbin_minfp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
		  &freq[j], &nextbinmin_nuFnu_zn[j], 
		  &nextbinmin_nuFnusyn_zn[j], &nextbinmin_nuFnussc_zn[j], 
		  &nextbinmin_nuFnuecd_zn[j], &nextbinmin_nuFnublr_zn[j], 
		  &nextbinmin_nuFnudt_zn[j], &nextbinmin_nuFnuup_zn[j], 
		  &nextbinmin_nuFnudo_zn[j]);

	 }   

	 fclose(nextbin_minfp);

	 nextbinfp = fopen(nextbin_file, "w");
	 if(nextbinfp == NULL) {
	   
	   printf("\n\n3 Unable to open file to write to: %s\n\n", nextbin_file);
	   exit(3);
	 }

	 winfp = fopen(old_file, "w");
	 if(winfp == NULL) {
	   
	   printf("\n\n3 Unable to open file to write to: %s\n\n", old_file);
	   exit(3);
	 }
	 	 
	   fprintf(nextbinfp, "%d %d %d\n", n_nextbin, obstime_nextbinmin, obstime_nextbinmax);
	   
	   fprintf(winfp, "%d %d %d\n", n, obstime_min, obstime_max);
	   
	 for (j = 0; j < FREQ_GRID; j++) {

	   nextbinnuFnu_zn[j] = nextbinnuFnu_zn[j] - (nextbinmin_nuFnu_zn[j] * LastRatio);
	   nextbinnuFnusyn_zn[j] = nextbinnuFnusyn_zn[j] - (nextbinmin_nuFnusyn_zn[j] * LastRatio);
	   nextbinnuFnussc_zn[j] = nextbinnuFnussc_zn[j] - (nextbinmin_nuFnussc_zn[j] * LastRatio);
	   nextbinnuFnuecd_zn[j] = nextbinnuFnuecd_zn[j] - (nextbinmin_nuFnuecd_zn[j] * LastRatio);
	   nextbinnuFnublr_zn[j] = nextbinnuFnublr_zn[j] - (nextbinmin_nuFnublr_zn[j] * LastRatio);
	   nextbinnuFnudt_zn[j] = nextbinnuFnudt_zn[j] - (nextbinmin_nuFnudt_zn[j] * LastRatio);
	   nextbinnuFnuup_zn[j] -= (nextbinmin_nuFnuup_zn[j] * LastRatio);
	   nextbinnuFnudo_zn[j] -= (nextbinmin_nuFnudo_zn[j] * LastRatio);
	   
	   nuFnu_zn[j] = nuFnu_zn[j] + (nextbinmin_nuFnu_zn[j] * LastRatio);
	   nuFnusyn_zn[j] = nuFnusyn_zn[j] + (nextbinmin_nuFnusyn_zn[j] * LastRatio);
	   nuFnussc_zn[j] = nuFnussc_zn[j] + (nextbinmin_nuFnussc_zn[j] * LastRatio);
	   nuFnuecd_zn[j] = nuFnuecd_zn[j] + (nextbinmin_nuFnuecd_zn[j] * LastRatio);
	   nuFnublr_zn[j] = nuFnublr_zn[j] + (nextbinmin_nuFnublr_zn[j] * LastRatio);
	   nuFnudt_zn[j] = nuFnudt_zn[j] + (nextbinmin_nuFnudt_zn[j] * LastRatio);
	   nuFnuup_zn[j] += (nextbinmin_nuFnuup_zn[j] * LastRatio);
	   nuFnudo_zn[j] += (nextbinmin_nuFnudo_zn[j] * LastRatio);

	     
	   fprintf(nextbinfp, "%e %e %e %e %e %e %e %e %e\n", freq[j], nextbinnuFnu_zn[j], nextbinnuFnusyn_zn[j], 
		   nextbinnuFnussc_zn[j], nextbinnuFnuecd_zn[j], nextbinnuFnublr_zn[j], nextbinnuFnudt_zn[j], 
		   nextbinnuFnuup_zn[j], nextbinnuFnudo_zn[j]);

	   fprintf(winfp, "%e %e %e %e %e %e %e %e %e\n", freq[j], nuFnu_zn[j], nuFnusyn_zn[j], 
		   nuFnussc_zn[j], nuFnuecd_zn[j], nuFnublr_zn[j], nuFnudt_zn[j], nuFnuup_zn[j], nuFnudo_zn[j]);
    	 }

	 fclose(winfp);
	 fclose(nextbinfp);

       } else {

	 printf("\n\nHole is bigger than one bin... Code is not designed to handle this...\n\n");
	 exit(3);

       }
       
       }
       
       winfp = fopen(old_file, "r");
       if(winfp == NULL) {

	 printf("\n\n3 Unable to open file to write to: %s\n\n", old_file);
	 exit(3);
       }

       fscanf(winfp, "%d %d %d\n", &n, &obstime_min, &obstime_max);

       for (j = 0; ((j < FREQ_GRID) && (!feof(winfp))); j++) {

         fscanf(winfp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &freq[j], &nuFnu_zn[j], 
		&nuFnusyn_zn[j], &nuFnussc_zn[j], &nuFnuecd_zn[j], &nuFnublr_zn[j], &nuFnudt_zn[j], 
		&nuFnuup_zn[j], &nuFnudo_zn[j]);

	 /*
         nuFnu_zn[j] = nuFnu_zn[j] / ((double) (window));
	 nuFnusyn_zn[j] = nuFnusyn_zn[j] / ((double) (window));
	 nuFnussc_zn[j] = nuFnussc_zn[j] / ((double) (window));
	 nuFnuecd_zn[j] = nuFnuecd_zn[j] / ((double) (window));
	 nuFnublr_zn[j] = nuFnublr_zn[j] / ((double) (window));
	 nuFnudt_zn[j] = nuFnudt_zn[j] / ((double) (window));
	 nuFnuup_zn[j] = nuFnuup_zn[j] / ((double) (window));
	 nuFnudo_zn[j] = nuFnudo_zn[j] / ((double) (window));
	 */

         nuFnu_zn[j] = nuFnu_zn[j] / ((double) (n));
         nuFnusyn_zn[j] = nuFnusyn_zn[j] / ((double) (n));
         nuFnussc_zn[j] = nuFnussc_zn[j] / ((double) (n));
	 nuFnuecd_zn[j] = nuFnuecd_zn[j] / ((double) (n));
	 nuFnublr_zn[j] = nuFnublr_zn[j] / ((double) (n));
         nuFnudt_zn[j] = nuFnudt_zn[j] / ((double) (n));
         nuFnuup_zn[j] = nuFnuup_zn[j] / ((double) (n));
         nuFnudo_zn[j] = nuFnudo_zn[j] / ((double) (n));

       }

       fclose(winfp);
       
       strcpy(filename_copy, old_file);
       filename_part = strtok(filename_copy, "/");
       filename_part = strtok(NULL, ".");
       sprintf(new_file, "%s/%s%s", window_dir_norm, filename_part, ".dat");
       printf("%s %d %s\n", old_file, n, new_file);

       winfp = fopen(new_file, "w");
       if(winfp == NULL) {
	    
	 printf("\n\n3 Unable to open file to write to: %s\n\n", window_filename);
	 exit(3);
       }

       for (j = 0; j < FREQ_GRID; j++) {

	 fprintf(winfp, "%e %e %e %e %e %e %e %e %e\n", freq[j], nuFnu_zn[j], nuFnusyn_zn[j], 
		 nuFnussc_zn[j], nuFnuecd_zn[j], nuFnublr_zn[j], nuFnudt_zn[j], nuFnuup_zn[j], nuFnudo_zn[j]);
	 
       }

       fclose(winfp);

     }

   }

   closedir(dir_test);

   
   // Interpolate the nuFnu values of the normalized files to obtain the 
   // lightcurves at given frequencies.
   //
   // Read the input file to store the values of frequencies at which 
   // lightcurves need to be calculated
   //
   fpin = fopen("input_freq-mfa.txt", "r");
   if (!fpin) {
     
     printf("\nCould not open the input_freq file\n");
     return -1;
     
   }
   
   for (i = 0; i < PPCSIZE && !feof(fpin); i++) {
     
     while(!feof(fpin) && fgetc(fpin) != ':');
     while(!feof(fpin) && fgetc(fpin) == ' ');
     fseek(fpin, -1, SEEK_CUR);
     fscanf(fpin, "%lf\n", &nu[i]);

   }

   fclose(fpin);

   nu1 = nu[0];
   nu2 = nu[1];
   nu3 = nu[2];
   nu4 = nu[3];
   nu5 = nu[4];
   nu6 = nu[5];
   
  printf("\nCalculating lightcurves at frequencies: %e, %e, %e, %e, %e, %e\n", 
	 nu1, nu2, nu3, nu4, nu5, nu6);
 
   dir_test_norm = opendir(window_dir_norm);
   printf("\nDirectory opened is %s\n", window_dir_norm);
   
   if(dir_test_norm == NULL) {
     printf("Could not open the nuFnu_sim_norm sub-directory.\n");
     exit (0);

   } else {
     
     while(NULL != (dir_entry_p = readdir(dir_test_norm))) {
       
       if ((strcmp(dir_entry_p->d_name, ".") == 0) || (strcmp(dir_entry_p->d_name, "..") == 0)) continue;
       
       window_filename = dir_entry_p->d_name;
       sprintf(new_file, "%s/%s", window_dir_norm, window_filename);

       for (i = 0; i < FREQ_GRID; i++) {
	 
	 nuFnu_zn[i] = 0.0;
	 freq[i] = 0.0;
       }
       
       strcpy(filename_copy, window_filename);
       filename_part = strtok(filename_copy, "_");
       filename_part = strtok(NULL, "_");
       filename_part = strtok(NULL, "_");
       file_sim_count = filename_part;
       filename_part = strtok(NULL, ".");
       file_time_count = filename_part;
       file_time_count[0] = ' ';
       sscanf(file_time_count, "%d", &obs_time);

       fpin = fopen(new_file, "r");
       if(fpin == NULL) {

	 printf("\n\n3 Unable to open file to read from: %s\n\n", new_file);
	 exit(3);
       }
       
       for (j = 0; j < FREQ_GRID && !feof(fpin); j++) {
	   
	   fscanf(fpin , "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &freq[j], &nuFnu_zn[j], &nuFnusyn_zn[j], 
		  &nuFnussc_zn[j], &nuFnuecd_zn[j], &nuFnublr_zn[j], &nuFnudt_zn[j], &nuFnuup_zn[j], &nuFnudo_zn[j]);

	 }

       fclose(fpin);

       for (j = 1; j < FREQ_GRID; j++) {

	 if ((nu1 > freq[j - 1]) && (nu1 < freq[j])) {

	   lc1 = exp(((log(nu1 / freq[j - 1]) * log(nuFnu_zn[j] / nuFnu_zn[j - 1])) / log(freq[j] / freq[j - 1]))) * nuFnu_zn[j - 1];
	   alpha1 = 1.0 - (log(nuFnu_zn[j] / nuFnu_zn[j - 1]) / log(freq[j] / freq[j - 1]));

	 }
	 
	 if ((nu2 > freq[j - 1]) && (nu2 < freq[j])) {
	  
	   lc2 = exp(((log(nu2 / freq[j - 1]) * log(nuFnu_zn[j] / nuFnu_zn[j - 1])) / log(freq[j] / freq[j - 1]))) * nuFnu_zn[j - 1];
	   alpha2 = 1.0 - (log(nuFnu_zn[j] / nuFnu_zn[j - 1]) / log(freq[j] / freq[j - 1]));

	 }
	 
	 if ((nu3 > freq[j - 1]) && (nu3 < freq[j])) {
	  
	   lc3 = exp(((log(nu3 / freq[j - 1]) * log(nuFnu_zn[j] / nuFnu_zn[j - 1])) / log(freq[j] / freq[j - 1]))) * nuFnu_zn[j - 1];
	   alpha3 = 1.0 - (log(nuFnu_zn[j] / nuFnu_zn[j - 1]) / log(freq[j] / freq[j - 1]));

	 }
	 
	 if ((nu4 > freq[j - 1]) && (nu4 < freq[j])) {
	   
	   lc4 = exp(((log(nu4 / freq[j - 1]) * log(nuFnu_zn[j] / nuFnu_zn[j - 1])) / log(freq[j] / freq[j - 1]))) * nuFnu_zn[j - 1];
	   alpha4 = 1.0 - (log(nuFnu_zn[j] / nuFnu_zn[j - 1]) / log(freq[j] / freq[j - 1]));

	 }
	 
	 if ((nu5 > freq[j - 1]) && (nu5 < freq[j])) {
	   
	   lc5 = exp(((log(nu5 / freq[j - 1]) * log(nuFnu_zn[j] / nuFnu_zn[j - 1])) / log(freq[j] / freq[j - 1]))) * nuFnu_zn[j - 1]; 
	   alpha5 = 1.0 - (log(nuFnu_zn[j] / nuFnu_zn[j - 1]) / log(freq[j] / freq[j - 1]));

	 }

	 if ((nu6 > freq[j - 1]) && (nu6 < freq[j])) {
	   
	   lc6 = exp(((log(nu6 / freq[j - 1]) * log(nuFnu_zn[j] / nuFnu_zn[j - 1])) / log(freq[j] / freq[j - 1]))) * nuFnu_zn[j - 1]; 
	   alpha6 = 1.0 - (log(nuFnu_zn[j] / nuFnu_zn[j - 1]) / log(freq[j] / freq[j - 1]));

	 }

       }

       sprintf(outlc, "lightcurve_%s.dat", file_sim_count);
       fplc = fopen(outlc, "a");
       if(fplc == NULL) {

	 printf("\n\n3 Unable to open file to write to: %s\n\n", outlc);
	 exit(3);
       }
     
       fprintf(fplc, "%e %e %e %e %e %e %e\n", ((double) (obs_time)), 
	       lc1, lc2, lc3, lc4, lc5, lc6);
       fclose(fplc);

       sprintf(outal, "spec_index_%s.dat", file_sim_count);
       fpal = fopen(outal, "a");
       if(fpal == NULL) {

	 printf("\n\n3 Unable to open file to write to: %s\n\n", outal);
	 exit(3);
       }
     
       fprintf(fpal, "%e %e %e %e %e %e %e\n", ((double) (obs_time)), alpha1, 
	       alpha2, alpha3, alpha4, alpha5, alpha6);
       fclose(fpal);

     }

   }

   closedir(dir_test_norm);

   if(sortFile(outlc) != 0) {
      
      printf("\nUnable to sort file.\n");
   
   } else {
      
      printf("\nFile %s sorted.\n", outlc);
   }
   
   if(sortFile(outal) != 0) {
      
      printf("\nUnable to sort file.\n\n");
   
   } else {
      
      printf("\nFile %s sorted.\n\n", outal);
   }

   fplc = fopen(outlc, "r");
   if(fplc == NULL) {
     
     printf("\n\n3 Unable to open file to write to: %s\n\n", outlc);
     exit(3);
   }

   fpal = fopen(outal, "r");
   if(fpal == NULL) {
     
     printf("\n\n3 Unable to open file to write to: %s\n\n", outlc);
     exit(3);
   }

   sprintf(outhic1, "hyst_lc1_%s.dat", file_sim_count);
   fphic1 = fopen(outhic1, "w");
   if(fphic1 == NULL) {
     
     printf("\n\n3 Unable to open file to write to: %s\n\n", outhic1);
     exit(3);
   }

   sprintf(outhic2, "hyst_lc2_%s.dat", file_sim_count);
   fphic2 = fopen(outhic2, "w");
   if(fphic2 == NULL) {
     
     printf("\n\n3 Unable to open file to write to: %s\n\n", outhic2);
     exit(3);
   }

   sprintf(outhic3, "hyst_lc3_%s.dat", file_sim_count);
   fphic3 = fopen(outhic3, "w");
   if(fphic3 == NULL) {
     
     printf("\n\n3 Unable to open file to write to: %s\n\n", outhic3);
     exit(3);
   }

   sprintf(outhic4, "hyst_lc4_%s.dat", file_sim_count);
   fphic4 = fopen(outhic4, "w");
   if(fphic4 == NULL) {
     
     printf("\n\n3 Unable to open file to write to: %s\n\n", outhic4);
     exit(3);
   }

   sprintf(outhic5, "hyst_lc5_%s.dat", file_sim_count);
   fphic5 = fopen(outhic5, "w");
   if(fphic5 == NULL) {
     
     printf("\n\n3 Unable to open file to write to: %s\n\n", outhic5);
     exit(3);
   }

   while(!feof(fplc) && !feof(fpal)) {
     
     fscanf(fplc, "%lf %lf %lf %lf %lf %lf %lf\n", &D_obstime, &lc1, &lc2, &lc3, &lc4, &lc5, &lc6);
     fscanf(fpal, "%lf %lf %lf %lf %lf %lf %lf\n", &D_obstime, &alpha1, &alpha2, &alpha3, &alpha4, &alpha5, &alpha6);
     fprintf(fphic1, "%e %e\n", lc1, alpha1);
     fprintf(fphic2, "%e %e\n", lc2, alpha2);
     fprintf(fphic3, "%e %e\n", lc3, alpha3);
     fprintf(fphic4, "%e %e\n", lc4, alpha4);
     fprintf(fphic5, "%e %e\n", lc5, alpha5);
   
   }

   fclose(fplc);
   fclose(fpal);
   fclose(fphic1);
   fclose(fphic2);
   fclose(fphic3);
   fclose(fphic4);
   fclose(fphic5);

   printf("\n\nComplete.\n\n");
   
   return 0;
}
