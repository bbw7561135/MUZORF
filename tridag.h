void tridag(double D_Aarr[], double D_Barr[], double D_Carr[], double D_Rarr[], double D_Uarr[], int n)
{
	int j;
	double bet, gam[500];
	
	
	if (D_Barr[0] == 0.0) {

	  printf("Error 1 in tridag");
	  exit(0);

	}

	D_Uarr[0] = D_Rarr[0] / (bet = D_Barr[0]);
	
	for (j = 1; j < n; j++) {
	  
	  gam[j] = D_Carr[j-1]/bet;
	  bet = D_Barr[j] - (D_Aarr[j] * gam[j]);

	  if (bet == 0.0) {

	    printf("Error 2 in tridag");
	    exit(0);

	  }

	  D_Uarr[j] = (D_Rarr[j] - (D_Aarr[j] * D_Uarr[j-1]))/bet;
	  
	}
	
	
	for (j = (n-2); j >= 0; j--) {
	  
	  D_Uarr[j] -= (gam[j+1] * D_Uarr[j+1]);
	  
	}
}

/* (C) Copr. 1986-92 Numerical Recipes Software 0@.1Y.. */
