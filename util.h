extern clock_t start;
extern char* logfile;
extern char* sOutputDirectory;
extern char* nfnfile;

double max(double x, double y)
{

  if (x > y) return (x);
  else return (y);
}


double min(double x, double y)
{
  if (x > y) return (y);
  else return (x);
}


double beta(double g)
{

  if (g < 1.0000001) return (0.0);
  else if (g < 1.0001) return (sqrt((2.0*(g - 1.0))));
  else return (sqrt((1.0 - (1.0/SQR(g)))));

}


inline int logMsg (char* logfile, char *label1, char *label2, int zone, int time_count, double oTime, char *sOutfile) 
{

	FILE *fp_log;

	fp_log = fopen(logfile, "a");
	if(!fp_log) {

		printf("\nCould not open the output file log. \n");
		return -1;
	}

	fprintf(fp_log, "%s = %d, time counter = %d, %s = %e, sOutfile = %s\n", label1, zone, time_count, label2, oTime, sOutfile);

	fclose(fp_log);

	return (0);
}


double myround(x)
	double x;
{
	double y, z, t_mid;

	y = floor(x);
	z = ceil(x);
	t_mid = (z - y)/2.0;
	if ((x - y) > t_mid) return(z);
	else return (y);
}


double dL_calc(double z, double H_0) 
{
	double x = 0.0, xs, sum = 0.0, sd_1, dx, s;
	double Omega0 = 0.3, OmegaL = 0.7;

	dx = 0.01 * min(z, 1.0);
	xs = s = 0.5 * dx;

	do {
		sd_1 = sqrt(((1.0 + (Omega0 * xs)) * SQR((1.0 + xs)) - 
			(OmegaL * xs * (2.0 + xs))));
		sum += (dx/sd_1);
		x += dx;
		xs = x + s;
	} while (xs < z);

	return (3.0857e19 * c * (1.0 + z) * sum / H_0);
}


double ne_num(double *e_engy, double *negam, int grid)
{
  double sum = 0.0, yofn[grid];

	int j;

	 yofn[0] = e_engy[0] * negam[0];

	for (j = 1; j < grid; j++) {

	  yofn[j] = e_engy[j] * negam[j];

	  sum += (0.5 * (yofn[j] + yofn[j - 1]) * (e_engy[j] - e_engy[j - 1]));

	}

	return(sum);

}


double ne_numtot(double *xarr, double *yarr, int grid) 
{
  double sum = 0.0;

  int j;

  for (j = 1; j < grid; j++) {

    sum += (0.5 * (yarr[j] + yarr[j - 1]) * (xarr[j] - xarr[j - 1]));

  }

  return(sum);

}


void handleInterrupt( int signal) {

  printf("\n\n\nCtrl-C pressed. Program will exit after saving the data in file. Please wait...\n\n\n");

  exitSimulation(start, logfile, -1, sOutputDirectory, nfnfile);
}

/*
double epsummain(double BLF, double betaphterm, double *phengy) {

  double sum = 0.0, expterm, denomterm, blfterm, yofep[CONTGRID];
  int j;

  if (betaphterm == 0.0) {
    blfterm = 1.0;
  } else {
    blfterm = BLF * betaphterm;
  }

  expterm = phengy[0] / AVTEMP;
  denomterm = exp(expterm) - 1.0;

  if (expterm >= 150.0) {
    yofep[0] = 0.0;

  } else if ((expterm >= 24.0) && (expterm < 150.0)) {
    yofep[0] = pow(phengy[0], 3.0) / exp(expterm);

  } else if ((expterm < 24.0) && (expterm > 1.0e-4)) {
    yofep[0] = pow(phengy[0], 3.0) / denomterm;

  } else {
    yofep[0] = pow(phengy[0], 3.0) / expterm;

  }

  if (yofep[0] < 1.0e-100)
    yofep[0] = 0.0;

  for (j = 1; j < CONTGRID; j++) {

    expterm = phengy[j] / AVTEMP;
    denomterm = exp(expterm) - 1.0;

    if (expterm >= 150.0) {
      yofep[j] = 0.0;

    } else if ((expterm >= 24.0) && (expterm < 150.0)) {
      yofep[j] = pow(phengy[j], 3.0) / exp(expterm);

    } else if ((expterm < 24.0) && (expterm > 1.0e-4)) {
      yofep[j] = pow(phengy[j], 3.0) / denomterm;

    } else {
      yofep[j] = pow(phengy[j], 3.0) / expterm;

    }

    if (yofep[j] < 1.0e-100) {
      yofep[j] = 0.0;
      continue;
    }

    sum += (0.5 * (yofep[j] + yofep[j - 1]) * (phengy[j] - phengy[j - 1]));

  }

  return (sum);

}
*/
