#include "uthash.h"

typedef struct {

	double dFreq;
	double dFlux;
	double dSYNFlux;
	double dSSCFlux;
	double dECDFlux;
	double dBLRFlux;
	double dDTFlux;
        double dUPPrevFlux;
        double dDOPrevFlux;

	UT_hash_handle hh;

} flux;

typedef struct {

	double dObsTime;
	flux *aFreqFlux;

	UT_hash_handle hh;

} shockrad;

// the hash map for storing freq and flux that is referenced by observer time
//
//shockrad *map = NULL;

extern shockrad *map;
extern double *syn_nu, D_Z;
extern int Ntot;
extern char startbuffer[25];
extern int f;
extern int r;


inline int photonDetector (double oTime, double dFreq, double dFlux, double dSYNFlux, double dSSCFlux, double dECDFlux, double dBLRFlux, double dDTFlux, double dUPPrevFlux, double dDOPrevFlux) {

	shockrad *item = NULL, *item2 = NULL;
	flux *ff = NULL, *ff2 = NULL;

	// check if the node for this dObsTime already exists
	//
	HASH_FIND(hh, map, &oTime, sizeof(double), item);

	if(item != NULL) {

		// update existing node for this dObsTime
		//
		//printf("\n\tupdating existing node for dObsTime = %e", oTime);

		HASH_FIND(hh, item->aFreqFlux, &dFreq, sizeof(double), ff2);

		if(ff2 != NULL) {

		  //printf("\n\t\tupdating existing node for dObsTime = %e, dFreq = %e, dFluz = %e", oTime, dFreq, dFlux);

		  // this freq exists; add to existing value
		  //
		  ff2->dFlux += dFlux;
		  ff2->dSSCFlux += dSSCFlux;
		  ff2->dECDFlux += dECDFlux;
		  ff2->dBLRFlux += dBLRFlux;
		  ff2->dDTFlux += dDTFlux;
		  ff2->dSYNFlux += dSYNFlux;
		  ff2->dUPPrevFlux += dUPPrevFlux;
		  ff2->dDOPrevFlux += dDOPrevFlux;

		} else {

		  //printf("\n\t\tadding new item to the map with dObsTime = %e, dFreq = %e, dFlux = %e", oTime, dFreq, dFlux);

			// this freq is new; create new node
			//
			ff2 = (flux*) malloc( sizeof(flux) );
			ff2->dFreq = dFreq;
			ff2->dFlux = dFlux;
			ff2->dSSCFlux = dSSCFlux;
			ff2->dSYNFlux = dSYNFlux;
			ff2->dECDFlux = dECDFlux;
			ff2->dBLRFlux = dBLRFlux;
			ff2->dDTFlux = dDTFlux;
			ff2->dUPPrevFlux = dUPPrevFlux;
			ff2->dDOPrevFlux = dDOPrevFlux;

			HASH_ADD(hh, item->aFreqFlux, dFreq, sizeof(double), ff2);
		}

	} else {

		// create new node for this dObsTime
		//
	        //printf("\nadding new item to the map with dObsTime = %e, dFreq = %e, dFlux = %e", oTime, dFreq, dFlux);

		item = (shockrad*) malloc( sizeof(shockrad) );
		item->dObsTime = oTime;
		item->aFreqFlux = NULL;

		HASH_ADD(hh, map, dObsTime, sizeof(double), item);

		ff = (flux*) malloc( sizeof(flux) );
		ff->dFreq = dFreq;
		ff->dFlux = dFlux;
		ff->dSSCFlux = dSSCFlux;
		ff->dSYNFlux = dSYNFlux;
		ff->dECDFlux = dECDFlux;
		ff->dBLRFlux = dBLRFlux;
		ff->dDTFlux = dDTFlux;
		ff->dUPPrevFlux = dUPPrevFlux;
		ff->dDOPrevFlux = dDOPrevFlux;

		HASH_ADD(hh, item->aFreqFlux, dFreq, sizeof(double), ff);
	}

	return(0);
}


inline int sendToDetector (int zone, int oTime, double aaPhotonDensity[FREQ_GRID], double dVolCoeff, int acclerationFlag, double *nuFnu_max, char *sOutputDirectory, char *nfnfile, double aaSYNDensity[FREQ_GRID], double aaSSCDensity[FREQ_GRID], double aaECDDensity[FREQ_GRID], double aaBLRDensity[FREQ_GRID], double aaDTDensity[FREQ_GRID], double aaUPPrevDensity[FREQ_GRID], double aaDOPrevDensity[FREQ_GRID], double angle, int timecountfsrs) {


	int j;
	double aFlux[FREQ_GRID];
	double aSYNFlux[FREQ_GRID], aSSCFlux[FREQ_GRID], aECDFlux[FREQ_GRID], aBLRFlux[FREQ_GRID], aDTFlux[FREQ_GRID];
	double aUPFlux[FREQ_GRID], aDOFlux[FREQ_GRID];

	for (j = 0; j < FREQ_GRID; j++) {

	        aFlux[j] = dVolCoeff * SQR(syn_nu[j]) * aaPhotonDensity[j];

		aSYNFlux[j] = dVolCoeff * SQR(syn_nu[j]) * aaSYNDensity[j];
		aSSCFlux[j] = dVolCoeff * SQR(syn_nu[j]) * aaSSCDensity[j];
		aECDFlux[j] = dVolCoeff * SQR(syn_nu[j]) * aaECDDensity[j];
		aBLRFlux[j] = dVolCoeff * SQR(syn_nu[j]) * aaBLRDensity[j];
		aDTFlux[j] = dVolCoeff * SQR(syn_nu[j]) * aaDTDensity[j];

		if (timecountfsrs == 0) {

		  aUPFlux[j] = aDOFlux[j] = 0.0;

		} else {

		  aUPFlux[j] = dVolCoeff * SQR(syn_nu[j]) * aaUPPrevDensity[j];
		  aDOFlux[j] = dVolCoeff * SQR(syn_nu[j]) * aaDOPrevDensity[j];

		}

		photonDetector(oTime, (D_Z * syn_nu[j]), max(1.0e-40, aFlux[j]), max(1.0e-40, aSYNFlux[j]), max(1.0e-40, aSSCFlux[j]), max(1.0e-40, aECDFlux[j]), max(1.0e-40, aBLRFlux[j]), max(1.0e-40, aDTFlux[j]), max(1.0e-40, aUPFlux[j]), max(1.0e-40, aDOFlux[j]));		
		
		if (acclerationFlag == FALSE) *nuFnu_max = max(aFlux[j], *nuFnu_max);
	}
	
	return (0);
}


inline int sendToDetector2 (int zone, int oTime, double aaSiPhotonDensity[FREQ_GRID], double aaPhotonDensity[FREQ_GRID], double dVolCoeff, int acclerationFlag, double *nuFnu_max, char *sOutputDirectory, char *nfnfile, double iZoneNumber, double aaSiSYNDensity[FREQ_GRID], double aaSYNDensity[FREQ_GRID], double aaSiSSCDensity[FREQ_GRID], double aaSSCDensity[FREQ_GRID], double aaSiECDDensity[FREQ_GRID], double aaECDDensity[FREQ_GRID], double aaSiBLRDensity[FREQ_GRID], double aaBLRDensity[FREQ_GRID], double aaSiDTDensity[FREQ_GRID], double aaDTDensity[FREQ_GRID], double aaSiUPPrevDensity[FREQ_GRID], double aaUPPrevDensity[FREQ_GRID], double aaSiDOPrevDensity[FREQ_GRID], double aaDOPrevDensity[FREQ_GRID], int timecountfsrs) {

	int j;
	double aFlux[FREQ_GRID], aEscapeFlux[FREQ_GRID];
	double aSYNFlux[FREQ_GRID], aSSCFlux[FREQ_GRID], aECDFlux[FREQ_GRID], aBLRFlux[FREQ_GRID], aDTFlux[FREQ_GRID];
	double aUPFlux[FREQ_GRID], aDOFlux[FREQ_GRID], aEscSYNFlux[FREQ_GRID], aEscSSCFlux[FREQ_GRID], aEscECDFlux[FREQ_GRID];
	double aEscBLRFlux[FREQ_GRID], aEscDTFlux[FREQ_GRID];
	double aEscUPFlux[FREQ_GRID], aEscDOFlux[FREQ_GRID];

	for (j = 0; j < FREQ_GRID; j++) {

		aFlux[j] = dVolCoeff * SQR(syn_nu[j]) * aaSiPhotonDensity[j];

                aSYNFlux[j] = dVolCoeff * SQR(syn_nu[j]) * aaSiSYNDensity[j];
                aSSCFlux[j] = dVolCoeff * SQR(syn_nu[j]) * aaSiSSCDensity[j];
                aECDFlux[j] = dVolCoeff * SQR(syn_nu[j]) * aaSiECDDensity[j];
                aBLRFlux[j] = dVolCoeff * SQR(syn_nu[j]) * aaSiBLRDensity[j];
                aDTFlux[j] = dVolCoeff * SQR(syn_nu[j]) * aaSiDTDensity[j];

		if (timecountfsrs == 0) {

                  aUPFlux[j] = aDOFlux[j] = 0.0;

                } else{

                  aUPFlux[j] = dVolCoeff * SQR(syn_nu[j]) * aaSiUPPrevDensity[j];
                  aDOFlux[j] = dVolCoeff * SQR(syn_nu[j]) * aaSiDOPrevDensity[j];

		}

		if (zone == (Ntot - 1) && (acclerationFlag == FALSE)) {

			aEscapeFlux[j] = dVolCoeff * iZoneNumber * SQR(syn_nu[j]) * aaPhotonDensity[j];

			aEscSYNFlux[j] = dVolCoeff * iZoneNumber * SQR(syn_nu[j]) * aaSYNDensity[j];
			aEscSSCFlux[j] = dVolCoeff * iZoneNumber * SQR(syn_nu[j]) * aaSSCDensity[j];
			aEscECDFlux[j] = dVolCoeff * iZoneNumber * SQR(syn_nu[j]) * aaECDDensity[j];
			aEscBLRFlux[j] = dVolCoeff * iZoneNumber * SQR(syn_nu[j]) * aaBLRDensity[j];
			aEscDTFlux[j] = dVolCoeff * iZoneNumber * SQR(syn_nu[j]) * aaDTDensity[j];
			
			aEscUPFlux[j] = dVolCoeff * iZoneNumber * SQR(syn_nu[j]) * aaUPPrevDensity[j];
			aEscDOFlux[j] = dVolCoeff * iZoneNumber * SQR(syn_nu[j]) * aaDOPrevDensity[j];

			// No cummulative addition as the radiation is coming only from the last zone,
			// also no time delay required.

			aFlux[j] += aEscapeFlux[j];

			aSYNFlux[j] += aEscSYNFlux[j];
			aSSCFlux[j] += aEscSSCFlux[j];
			aECDFlux[j] += aEscECDFlux[j];
			aBLRFlux[j] += aEscBLRFlux[j];
			aDTFlux[j] += aEscDTFlux[j];

			aUPFlux[j] += aEscUPFlux[j];
			aDOFlux[j] += aEscDOFlux[j];
		}

		photonDetector(oTime, (D_Z * syn_nu[j]), max(1.0e-40, aFlux[j]), max(1.0e-40, aSYNFlux[j]), max(1.0e-40, aSSCFlux[j]), max(1.0e-40, aECDFlux[j]), max(1.0e-40, aBLRFlux[j]), max(1.0e-40, aDTFlux[j]), max(1.0e-40, aUPFlux[j]), max(1.0e-40, aDOFlux[j]));

		if (acclerationFlag == FALSE) *nuFnu_max = max(aFlux[j], *nuFnu_max);
	}

	return (0);
}

void displayTelescopeData () {

	shockrad *item;
	flux *ff;

    for(item = map; item != NULL; item = item->hh.next) {

        printf("\n\nobs time = %e", item->dObsTime);

		for(ff = item->aFreqFlux; ff != NULL; ff = ff->hh.next) {

                  printf("\n\tdFreq = %e & dFlux = %e, dSYNFlux = %e, dSSCFlux = %e, dECDFlux = %e, dBLRFlux = %e, dDTFlux = %e, dUPPrevFlux = %e, dDOPrevFlux = %e", ff->dFreq, ff->dFlux, ff->dSYNFlux, ff->dSSCFlux, ff->dECDFlux, ff->dBLRFlux, ff->dDTFlux, ff->dUPPrevFlux, ff->dDOPrevFlux);

		}

		printf("\n\n");
    }

}


void writeTelescopeData (char *sOutputDirectory, char *nfnfile) {

	shockrad *item;
	flux *ff;
	FILE *fp;
	char* sOutfile = (char*) malloc(sizeof(char) * 1024);

    for(item = map; item != NULL; item = item->hh.next) {

		sprintf(sOutfile, "%s/%s%d%s", sOutputDirectory, nfnfile, (int) item->dObsTime, ".dat");
		fp = fopen(sOutfile, "w");

		if(fp) {

			for(ff = item->aFreqFlux; ff != NULL; ff = ff->hh.next) {

			  fprintf(fp, "%e %e %e %e %e %e %e %e %e\n", ff->dFreq, ff->dFlux, ff->dSYNFlux, ff->dSSCFlux, ff->dECDFlux, 
				  ff->dBLRFlux, ff->dDTFlux, ff->dUPPrevFlux, ff->dDOPrevFlux);

			}

			fclose(fp);

		} else {

			printf("\nERROR: Unable to open output file: %s", sOutfile);
		}

    }

}

inline int sendToRustyTelescope (int zone, int oTime, double aaPhotonDensity[][FREQ_GRID], double dVolCoeff, int acclerationFlag, double *nuFnu_max, char *sOutputDirectory, char *nfnfile) {

	char* sOutfile = (char*) malloc(sizeof(char) * 1024);
	int j;
	FILE *fp_nfn;
	double freq, nuFnu[FREQ_GRID], aFlux[FREQ_GRID];

	sprintf(sOutfile, "%s/%s%d%s", sOutputDirectory, nfnfile, oTime, ".dat");

	fp_nfn = fopen(sOutfile, "r+");
	if (fp_nfn) {

		for (j = 0; ((j < FREQ_GRID) && (!feof(fp_nfn))); j++) {

			fscanf(fp_nfn, "%lf %lf\n", &freq, &nuFnu[j]);
		}

		fseek(fp_nfn, 0, SEEK_SET);

	} else {

		fp_nfn = fopen(sOutfile, "w");
		if(!fp_nfn) {

			printf("\nCould not open the output file nfn in nuFnu_up case. \n");
			return -1;

		}

		for (j = 0; j < FREQ_GRID; j++) {

			nuFnu[j] = 0.0;
		}

	}

	for (j = 0; j < FREQ_GRID; j++) {

		aFlux[j] = dVolCoeff * SQR(syn_nu[j]) * aaPhotonDensity[zone][j];

		nuFnu[j] += aFlux[j];

		fprintf(fp_nfn, "%e %e\n", (D_Z * syn_nu[j]), max(1.0e-40, nuFnu[j]));

		if (acclerationFlag == FALSE) *nuFnu_max = max(nuFnu[j], *nuFnu_max);
	}

	fclose(fp_nfn);

	return (0);
}


inline int sendToRustyTelescope2 (int zone, int oTime, double aaSiPhotonDensity[][FREQ_GRID], double aaPhotonDensity[][FREQ_GRID], double dVolCoeff, int acclerationFlag, double *nuFnu_max, char *sOutputDirectory, char *nfnfile) {

	char* sOutfile = (char*) malloc(sizeof(char) * 1024);
	int j;
	FILE *fp_nfn;
	double freq, nuFnu[FREQ_GRID], aFlux[FREQ_GRID], aEscapeFlux[FREQ_GRID];

	sprintf(sOutfile, "%s/%s%d%s", sOutputDirectory, nfnfile, oTime, ".dat");

	fp_nfn = fopen(sOutfile, "r+");
	if (fp_nfn) {

		for (j = 0; ((j < FREQ_GRID) && (!feof(fp_nfn))); j++) {

			fscanf(fp_nfn, "%lf %lf\n", &freq, &nuFnu[j]);
		}

		fseek(fp_nfn, 0, SEEK_SET);

	} else {

		fp_nfn = fopen(sOutfile, "w");
		if(!fp_nfn) {

			printf("\nCould not open the output file nfn in nuFnu_up case. \n");
			return -1;

		}

		for (j = 0; j < FREQ_GRID; j++) {

			nuFnu[j] = 0.0;
		}

	}

	for (j = 0; j < FREQ_GRID; j++) {

		aFlux[j] = dVolCoeff * SQR(syn_nu[j]) * aaSiPhotonDensity[zone][j];


		if (zone == (Ntot - 1) && (acclerationFlag == FALSE)) {

			aEscapeFlux[j] = dVolCoeff * SQR(syn_nu[j]) * aaPhotonDensity[zone][j];
			// No cummulative addition as the radiation is coming only from the last zone,
			// also no time delay required.

			nuFnu[j] += (aFlux[j] + aEscapeFlux[j]);

		} else {

			nuFnu[j] += aFlux[j];

		}


		nuFnu[j] += aFlux[j];

		fprintf(fp_nfn, "%e %e\n", (D_Z * syn_nu[j]), max(1.0e-40, nuFnu[j]));

		if (acclerationFlag == FALSE) *nuFnu_max = max(nuFnu[j], *nuFnu_max);
	}

	fclose(fp_nfn);

	return (0);
}

inline int exitSimulation (time_t start, char *logfile, int iStatus, char *sOutputDirectory, char *nfnfile) {

  time_t end, end2;
  double elapsed, elapsed2;
  FILE *fp_log;

  time_t timer;
  char buffer[25];
  struct tm* tm_info;


  end = time(NULL);
  elapsed = end - start;
  printf("Total elapsed time = %e seconds\n", elapsed);

  fp_log = fopen(logfile, "a");
  if(!fp_log) {

    printf("\nCould not open the output file log. \n");
    return -1;
  }

  fprintf(fp_log, "\n\nTotal elapsed time = %e seconds\n\n", elapsed);

  fclose(fp_log);
  printf("\nEnd of Simulation\n");

  printf("\nCreating output files\n");

  //displayTelescopeData();
  writeTelescopeData(sOutputDirectory, nfnfile);

  end2 = time(NULL);
  elapsed2 = end2 - end;

  time(&timer);
  tm_info = localtime(&timer);

  strftime(buffer, 25, "%Y:%m:%d %H:%M:%S", tm_info);

  printf("\n\n\n");
  printf("\nSimulation started: %s", startbuffer);
  printf("\nSimulation ended: %s", buffer);
  printf("\nTotal elapsed time = %e seconds", elapsed);
  printf("\nTotal time for writing output file = %e seconds", elapsed2);
  printf("\n\nExiting program.");
  printf("\n\n\n");

  exit(iStatus);
}
