#define LL 3.21754e-6       // Value of (h*nu_lower)/(m_e*c^2)
                            // nu_lower = 10^(14.6) Hz (Liu & Bai 2006)
#define UL 2.55579e-4       // Value of (h*nu_upper)/(m_e*c^2)
                            // nu_upper = 10^(16.5) Hz (Liu & Bai 2006)
#define ULGRID 60           // Value of grid points between LL & UL
#define AVTEMP 1.68404e-05  // Value of THETA_D for T = 10^5 K (Liu & Bai 2006)
#define LINECOEFF 7.31563e-08 // Value of 1/(c*mc^2*555.77)
#define CONTCOEFF 4.06581e-05 // Value of 1/(c*mc^2)
#define GDOTBLRCOEFF 4.70061e-14  // Value of 3*c*pi*SIGMA_T/4
//#define NDOTBLRCOEFF 7.48125e-15  // Value of 3*c*SIGMA_T/8
#define NDOTBLRTHCOEFF 1.995e-14 // Value of c*SIGMA_T
#define THCONTCOEFF 5.22297e-19  // Value of (pi*THETA_D)^4/15.0
#define THLINECOEFF 555.77   // Value of total NVLINE
#define GDOTBLRTHCOEFF 1.2535e-13  // Value of 2*pi*c*SIGMA_T
#define GDOTTHCOMPCOEFF 1.67133e-13 // Value of 8*pi*c*SIGMA_T/3
#define ETASTEP 51
#define LINEGRID 35
#define CONTGRID 30
#define absaccur 3.0e-11


// Functions to make the grid for eta_ph to be used in the calculation of
// BLR intensity profile, EC scattering & electron energy loss due to BLR
// photons. The grid has been defined in the comoving frame of the plasma.
// Subroutine taken from Numerical Recipes in C, 2nd ed., 1992, pg 148.
// The abscissa & weights have been taken from Abr. & Stegun, for W(x) = 1, &
// N = 48.
//
void lingrid(double etamin, double etamax, double *etaarr)
{

  int j;
  double step;

  step = (etamax - etamin) / ((double) (ETASTEP));

  for (j = 0; j < ETASTEP; j++) {

    etaarr[j] = etamin;

    etamin += step;

  }

  return;

}


void gausgrid(double etmin, double etmax, int grpt, double *x, double *aetagausxmpdx, double *aetagausxmmdx, int *gausindex)
{

  FILE *fp_etaarrval;
  int j, k = 0;
  double xr, xm, dx;

  xm = 0.5 * (etmax + etmin);
  xr = 0.5 * (etmax - etmin);

  fp_etaarrval = fopen("etaarrval_48-10.dat", "a");
  if(!fp_etaarrval) {

    printf("\nCould not open the output file etaarrval_48-10. \n");
    return;
  }


  for (j = 1; j <= grpt; j++) {

    dx = xr * x[j];

    aetagausxmpdx[k] = xm + dx;
    fprintf(fp_etaarrval, "%e\n", (xm + dx));

    aetagausxmmdx[k] = xm - dx;
    fprintf(fp_etaarrval, "%e\n", (xm - dx));

    k++;

  }

  *gausindex = k;

  fclose(fp_etaarrval);

  return;

}


void gauleg(double lowlim, double uplim, double *absci, double *w, int grpt)
{

  int m, j, i;
  double z1, z, xm, xl, pp, p3, p2, p1, grterm;

  m = (grpt + 1) / 2;

  grterm = grpt + 0.5;

  xm = 0.5 * (uplim + lowlim);
  xl = 0.5 * (uplim - lowlim);

  for (i = 1; i <= m; i++) {

    z = cos((3.141592654 * (i - 0.25)) / grterm);

    do {
      p1 = 1.0;
      p2 = 0.0;

      for (j = 1; j <= grpt; j++) {
	p3 = p2;
	p2 = p1;
	p1 = ((((2.0 * j) - 1.0) * z * p2) - ((j - 1.0) * p3)) / j;

      }

      pp = (grpt * ((z * p1) - p2)) / (SQR(z) - 1.0);
      z1 = z;
      z = z1 - (p1 / pp);

    } while (fabs(z - z1) > absaccur);

    absci[i] = xm - (xl * z);
    absci[grpt + 1 - i] = xm + (xl * z);
    w[i] = (2.0 * xl) / ((1.0 - SQR(z)) * SQR(pp));
    w[grpt + 1 - i] = w[i];

  }

}


double integ_exam(double xval)
{

  double expr;

  expr = exp(-xval) / sqrt(1.0 - SQR(xval));

  return(expr);

}


double qgaus(double (*fungaus)(double), double etmin, double etmax, double *x, double *w)
{

  int j;
  double xr, xm, dx, sum = 0.0, R1, R2;

  xm = 0.5 * (etmax + etmin);
  xr = 0.5 * (etmax - etmin);

  for (j = 1; j <= 24; j++) {

    dx = xr * x[j];
    sum += (w[j] * ((*fungaus)(xm + dx) + (*fungaus)(xm - dx)));

  }

  return(sum *= xr);

}


// Start of calculation of BLR intensity as a function of emission distance,
// zht, and angle at which BLR photons enter the emission region. Both
// quantities are in the AGN frame.
//
void blrintens(double LBLR, double optdep, double covfac, double zht, double rin_blr, double rout_blr, double phcos, double Gamma_bulk, double lineinteg, double continteg, double *lineI, double *contI)
{
  double lmin = 0., lmax = 0., lmax2 = 0., mustar = 0., lineI2, contI2, contI3;
  double logstep1 = 0., Beta_bulk = 0., lineterm = 0., contterm = 0., lineI2a;
  double contI2a, r = 0., logstep3b = 0., logstep2b = 0., linep1 = 0.;
  double lmin2 = 1e-3, logstep2 = 0., logstep2a = 0., logstep3a = 0.;
  double lineI2b, contI2b, lineI3a, contI3a, lineI3b, contI3b, q = 0.;
  double logstep3d = 0., lineI3d, contI3d;
  double contp1 = 0., term1 = 0., term2 = 0., lstep = 0., rootl1 = 0.;
  double p = 1.5, s = 1.0, lineI1, contI1, mustaroutcr = 0., rootl2 = 0.;
  double rootmin1 = 0., rootmin2 = 0., rootmax1 = 0., rootmax2 = 0.;
  double len[GS], jline[GS], jcont[GS];
  double loutcr = 0., mustarincr = 0., lincr = 0., lmin_min = 0., lineI3;
  double lmin_max = 0., lmax_min = 0., lmax_max = 0., lstep3a = 0.;
  double lstep3b = 0., lstep3d = 0.;

  int j;

  Beta_bulk = beta(Gamma_bulk);

  q = 1.0 / 3.0;
  linep1 = (2.0 * q) - p - 2.0;
  contp1 = s + 2.0;

  if ((q == 1.0/3.0) && (p == 1.5) && (s == 1.0)) {
    linep1 = -2.83333;
    contp1 = 3.0;
  }

  mustar = (phcos + Beta_bulk)/(1.0 + (Beta_bulk * phcos));

  if (mustar > .99999999) mustar = .99999999;
  if (mustar < -.99999999) mustar = -.99999999;

  lineI1 = lineI2 = lineI2a = lineI2b = lineI3 = lineI3a = 0.0;
  lineI3b = lineI3d = 0.0;
  contI1 = contI2 = contI2a = contI2b = contI3 = contI3a = 0.0;
  contI3b = contI3d = 0.0;

  if (zht <= rin_blr) {
    term1 = SQR(mustar) + SQR(rin_blr / zht) - 1.0;
    rootmin1 = (zht * mustar) - (zht * sqrt(term1));
    rootmin2 = (zht * mustar) + (zht * sqrt(term1));
    lmin = max(rootmin1, rootmin2);

    term2 = SQR(mustar) + SQR(rout_blr / zht) - 1.0;
    rootmax1 = (zht * mustar) - (zht * sqrt(term2));
    rootmax2 = (zht * mustar) + (zht * sqrt(term2));
    lmax = max(rootmax1, rootmax2);

    if ((term1 < 0.0) || (term2 < 0.0)) {
      printf("\nSomething wrong in blrintens, Pos.1:\n");
      printf("term1 = %e, term2 = %e\n", term1, term2);
      exit(Exit_success);
    }

    if (lmax <= lmin) {
      printf("\nSomething wrong in blrintens, Pos.1:\n");
      printf("lmax = %e, lmin = %e\n", lmax, lmin);
      exit(Exit_success);
    }

    logstep1 = pow((lmax / lmin), 1.0/((double) (GS - 1)));
    lstep = lmin;

    for (j = 0; j < GS; j++) {

      len[j] = lstep;
      r = sqrt(SQR(zht) - (2.0 * zht * len[j] * mustar) + SQR(len[j]));
      lineterm = LBLR * pow(r, linep1);
      jline[j] = lineterm / (JCOEFF * lineinteg);

      contterm = LBLR * optdep / (pow(r, contp1) * covfac);
      jcont[j] = contterm / (JCOEFF * continteg);

      lstep *= logstep1;

    }

    for (j = 1; j < GS; j++) {

      lineI1 += (0.5 * (jline[j] + jline[j - 1]) * (len[j] - len[j - 1]));
      contI1 += (0.5 * (jcont[j] + jcont[j - 1]) * (len[j] - len[j - 1]));

    }

    *lineI = lineI1;
    *contI = contI1;

  } else if ((zht > rin_blr) && (zht <= rout_blr)) {
    lincr = sqrt((SQR(zht) - SQR(rin_blr)));
    mustarincr = lincr / zht;

    if (mustar > mustarincr) {
      term1 = SQR(mustar) + SQR(rin_blr / zht) - 1.0;
      rootl1 = (zht * mustar) - (zht * sqrt(term1));
      rootl2 = (zht * mustar) + (zht * sqrt(term1));
      lmax2 = max(rootl1, rootl2);
      lmin = min(rootl1, rootl2);

      if (term1 < 0.0) {
	printf("\nSomething wrong in blrintens, Pos.2, 1st integral\n");
	printf("term1 = %e, zht = %e, mustar = %e, phcos = %e\n",
	       term1, zht, mustar, phcos);
	printf("sqmu = %e, sqrinz = %e\n", SQR(mustar), SQR(rin_blr/zht));
	exit(Exit_success);
      }

      if (lmin <= lmin2) {
	printf("\nSomething wrong in blrintens, Pos.2, 1st integral\n");
	printf("rootl1 = %e, rootl2 = %e, lmin = %e, lmin2 = %e\n",
	       rootl1, rootl2, lmin, lmin2);
	exit(Exit_success);
      }

      logstep2 = pow((lmin / lmin2), 1.0/((double) (GS - 1)));
      lstep = lmin2;

      for (j = 0; j < GS; j++) {

	len[j] = lstep;
	r = sqrt(SQR(zht) - (2.0 * zht * len[j] * mustar) + SQR(len[j]));
	lineterm = LBLR * pow(r, linep1);
	jline[j] = lineterm / (JCOEFF * lineinteg);

	contterm = LBLR * optdep / (pow(r, contp1) * covfac);
	jcont[j] = contterm / (JCOEFF * continteg);

	lstep *= logstep2;

      }

      for (j = 1; j < GS; j++) {

	lineI2 += (0.5 * (jline[j] + jline[j - 1]) * (len[j] - len[j - 1]));
	contI2 += (0.5 * (jcont[j] + jcont[j - 1]) * (len[j] - len[j - 1]));

      }
    } else {
      lmax2 = lmin2;
      lineI2 = 0.0;
      contI2 = 0.0;
    }

    term2 = SQR(mustar) + SQR(rout_blr / zht) - 1.0;
    rootmax1 = (zht * mustar) - (zht * sqrt(term2));
    rootmax2 = (zht * mustar) + (zht * sqrt(term2));
    lmax = max(rootmax1, rootmax2);

    if (term2 < 0.0) {
      printf("\nSomething wrong in blrintens, Pos.2, 2nd integral\n");
      printf("term1 = %e, term2 = %e\n", term1, term2);
      printf("sqmu = %e, sqrinz = %e, sqrouz = %e\n",
	     SQR(mustar), SQR(rin_blr/zht), SQR(rout_blr/zht));
      exit(Exit_success);
    }

    if ((lmax <= lmax2) || (lmax2 < 0.0)) {
      printf("\nSomething wrong in blrintens, Pos.2, 2nd integral\n");
      printf("rootmax1 = %e, rootmax2 = %e\n", rootmax1, rootmax2);
      printf("lmax = %e, lmax2 = %e\n", lmax, lmax2);
      exit(Exit_success);
    }

    logstep2a = pow((lmax / lmax2), 1.0/((double) (GS - 1)));
    lstep = lmax2;

    for (j = 0; j < GS; j++) {

      len[j] = lstep;
      r = sqrt(SQR(zht) - (2.0 * zht * len[j] * mustar) + SQR(len[j]));
      lineterm = LBLR * pow(r, linep1);
      jline[j] = lineterm / (JCOEFF * lineinteg);

      contterm = LBLR * optdep / (pow(r, contp1) * covfac);
      jcont[j] = contterm / (JCOEFF * continteg);

      lstep *= logstep2a;

    }

    for (j = 1; j < GS; j++) {

      lineI2a += (0.5 * (jline[j] + jline[j - 1]) * (len[j] - len[j - 1]));
      contI2a += (0.5 * (jcont[j] + jcont[j - 1]) * (len[j] - len[j - 1]));

    }

    *lineI = lineI2 + lineI2a;
    *contI = contI2 + contI2a;

  } else if ((zht > rout_blr) && ((rout_blr / zht) >= 1.0e-2)) {
    lincr = sqrt((SQR(zht) - SQR(rin_blr)));
    mustarincr = lincr / zht;

    loutcr = sqrt(SQR(zht) - SQR(rout_blr));
    mustaroutcr = loutcr / zht;

    if ((mustar > mustarincr) && (mustar > mustaroutcr)) {
      term1 = SQR(mustar) + SQR(rin_blr / zht) - 1.0;
      rootmin1 = (zht * mustar) - (zht * sqrt(term1));
      rootmin2 = (zht * mustar) + (zht * sqrt(term1));
      lmin_min = min(rootmin1, rootmin2);
      lmin_max = max(rootmin1, rootmin2);

      if (term1 < 0.0) {
	printf("\nSomething wrong in blrintens, Pos.3, integral 3a\n");
	printf("term1 = %e, zht = %e, mustar = %e, phcos = %e\n",
	       term1, zht, mustar, phcos);
	printf("sqmu = %e, sqrinz = %e\n", SQR(mustar), SQR(rin_blr/zht));
	exit(Exit_success);
      }

      if ((lmin_max <= lmin_min) || (lmin_min <= 0.0)) {
	printf("\nSomething wrong in blrintens, Pos.3, integral 3a\n");
	printf("rootmin1 = %e, rootmin2 = %e, lmin_max = %e, lmin_min = %e\n",
	       rootmin1, rootmin2, lmin_max, lmin_min);
	exit(Exit_success);
      }

      term2 = SQR(mustar) + SQR(rout_blr / zht) - 1.0;
      rootmax1 = (zht * mustar) - (zht * sqrt(term2));
      rootmax2 = (zht * mustar) + (zht * sqrt(term2));
      lmax_min = min(rootmax1, rootmax2);
      lmax_max = max(rootmax1, rootmax2);

      if (term2 < 0.0) {
	printf("\nSomething wrong in blrintens, Pos.3, integral 3b\n");
	printf("term1 = %e, term2 = %e\n", term1, term2);
	printf("sqmu = %e, sqrinz = %e, sqrouz = %e\n",
	       SQR(mustar), SQR(rin_blr/zht), SQR(rout_blr/zht));
	exit(Exit_success);
      }

      if ((lmax_max <= lmax_min) || (lmax_min <= 0.0)) {
	printf("\nSomething wrong in blrintens, Pos.3, integral 3b\n");
	printf("rootmin1 = %e, rootmin2 = %e, rootmax1 = %e, rootmax2 = %e\n",
	       rootmin1, rootmin2, rootmax1, rootmax2);
	printf("lmax_min = %e, lmax_max = %e\n", lmax_min, lmax_max);
	exit(Exit_success);
      }

      if (lmax_min > lmin_min) {
	logstep3a = pow((lmax_min / lmin_min), 1.0/((double) (GS - 1)));
	lstep3a = lmin_min;
      } else {
	logstep3a = pow((lmin_min / lmax_min), 1.0/((double) (GS - 1)));
	lstep3a = lmax_min;
      }

      if (lmax_max > lmin_max) {
	logstep3b = pow((lmax_max / lmin_max), 1.0/((double) (GS - 1)));
	lstep3b = lmin_max;
      } else {
	logstep3b = pow((lmin_max / lmax_max), 1.0/((double) (GS - 1)));
	lstep3b = lmax_max;
      }

      for (j = 0; j < GS; j++) {

	len[j] = lstep3a;
	r = sqrt(SQR(zht) - (2.0 * zht * len[j] * mustar) + SQR(len[j]));
	lineterm = LBLR * pow(r, linep1);
	jline[j] = lineterm / (JCOEFF * lineinteg);

	contterm = LBLR * optdep / (pow(r, contp1) * covfac);
	jcont[j] = contterm / (JCOEFF * continteg);

	lstep3a *= logstep3a;

      }

      for (j = 1; j < GS; j++) {

	lineI3a += (0.5 * (jline[j] + jline[j - 1]) * (len[j] - len[j - 1]));
	contI3a += (0.5 * (jcont[j] + jcont[j - 1]) * (len[j] - len[j - 1]));

      }

      for (j = 0; j < GS; j++) {

	len[j] = lstep3b;
	r = sqrt(SQR(zht) - (2.0 * zht * len[j] * mustar) + SQR(len[j]));
	lineterm = LBLR * pow(r, linep1);
	jline[j] = lineterm / (JCOEFF * lineinteg);

	contterm = LBLR * optdep / (pow(r, contp1) * covfac);
	jcont[j] = contterm / (JCOEFF * continteg);

	lstep3b *= logstep3b;

      }

      for (j = 1; j < GS; j++) {

	lineI3b += (0.5 * (jline[j] + jline[j - 1]) * (len[j] - len[j - 1]));
	contI3b += (0.5 * (jcont[j] + jcont[j - 1]) * (len[j] - len[j - 1]));

      }

      *lineI = lineI3a + lineI3b;
      *contI = contI3a + contI3b;

    } else if ((mustar <= mustarincr) && (mustar > mustaroutcr)){
      term2 = SQR(mustar) + SQR(rout_blr / zht) - 1.0;
      rootmax1 = (zht * mustar) - (zht * sqrt(term2));
      rootmax2 = (zht * mustar) + (zht * sqrt(term2));
      lmax_min = min(rootmax1, rootmax2);
      lmax_max = max(rootmax1, rootmax2);

      if (term2 < 0.0) {
	printf("\nSomething wrong in blrintens, Pos.3, integral 3d\n");
	printf("term2 = %e, sqmu = %e, sqrouz = %e\n", term2, SQR(mustar),
	       SQR(rout_blr/zht));
	exit(Exit_success);
      }

      if ((lmax_max <= lmax_min) || (lmax_min <= 0.0)) {
	printf("\nSomething wrong in blrintens, Pos.3, integral 3d\n");
	printf("rootmax1 = %e, rootmax2 = %e\n", rootmax1, rootmax2);
	printf("lmax_min = %e, lmax_max = %e\n", lmax_min, lmax_max);
	exit(Exit_success);
      }

      logstep3d = pow((lmax_max / lmax_min), 1.0/((double) (GS - 1)));
      lstep3d = lmax_min;

      for (j = 0; j < GS; j++) {

	len[j] = lstep3d;
	r = sqrt(SQR(zht) - (2.0 * zht * len[j] * mustar) + SQR(len[j]));
	lineterm = LBLR * pow(r, linep1);
	jline[j] = lineterm / (JCOEFF * lineinteg);

	contterm = LBLR * optdep / (pow(r, contp1) * covfac);
	jcont[j] = contterm / (JCOEFF * continteg);

	lstep3d *= logstep3d;

      }

      for (j = 1; j < GS; j++) {

	lineI3d += (0.5 * (jline[j] + jline[j - 1]) * (len[j] - len[j - 1]));
	contI3d += (0.5 * (jcont[j] + jcont[j - 1]) * (len[j] - len[j - 1]));

      }

      *lineI = lineI3d;
      *contI = contI3d;

    } else {
      *lineI = 0.0;
      *contI = 0.0;

    }

  } else {
    *lineI = 0.0;
    *contI = 0.0;

  }

  return;

}


// Test function to see the integration result of eta_ph
//
double integ_fun(double xval, double Gamma_bulk, double intenval)
{

  double expr, Beta_bulk, linedenom, bphterm;

  Beta_bulk = beta(Gamma_bulk);

  bphterm = 1.0 + (Beta_bulk * xval);
  linedenom = pow((Gamma_bulk * bphterm), 3.0);

  expr = intenval / linedenom;

  return(expr);

}


double etagaus_integ(double (*fungaus)(double, double, double), double etmin, double etmax, double *x, double *w, double Gamma_bulk, double LBLR, double optdep, double covfac, double zht, double rin_blr, double rout_blr, double lineinteg, double continteg, double *lineI, double *contI, double lineinten[], double continten[], char *Izfile, int index)
{

  int j, y, k = 0;
  double xr, xm, dx, sum = 0.0, phcos, phcosstar, term1, term2;

  xm = 0.5 * (etmax + etmin);
  xr = 0.5 * (etmax - etmin);

  for (j = 1; j <= 24; j++) {

    dx = xr * x[j];
    phcos = xm + dx;
    blrintens(LBLR, optdep, covfac, zht, rin_blr, rout_blr, phcos, Gamma_bulk, lineinteg, continteg, lineI, contI);
    lineinten[k] = *lineI;
    continten[k] = *contI;

    term1 = (*fungaus)((xm + dx), Gamma_bulk, lineinten[k]);
    k++;

    phcos = xm - dx;
    blrintens(LBLR, optdep, covfac, zht, rin_blr, rout_blr, phcos, Gamma_bulk, lineinteg, continteg, lineI, contI);
    lineinten[k] = *lineI;
    continten[k] = *contI;

    term2 = (*fungaus)((xm - dx), Gamma_bulk, lineinten[k]);
    k++;

    sum += (w[j] * (term1 + term2));

  }

  return(sum *= xr);

}

double etainteg_test(double Gamma_bulk, double *phcos, double *ablrinten, int n)
{
  double sum = 0.0, yofeta[n], bphterm, Beta_bulk, linedenom;
  int j;

  Beta_bulk = beta(Gamma_bulk);

  bphterm = 1.0 + (Beta_bulk * phcos[0]);
  linedenom = pow((Gamma_bulk * bphterm), 3.0);
  yofeta[0] = ablrinten[0] / linedenom;

  for (j = 1; j < n; j++) {

    bphterm = 1.0 + (Beta_bulk * phcos[j]);
    linedenom = pow((Gamma_bulk * bphterm), 3.0);
    yofeta[j] = ablrinten[j] / linedenom;

    sum += (0.5 * (yofeta[j] + yofeta[j - 1]) * (phcos[j] - phcos[j - 1]));

  }

  return(sum);

}


void nphblrtotuph(double *Nnuline, double BLF, double lineinteg, double continteg, double LBLR, double optdep, double covfac, double zht, double rin_blr, double rout_blr, double *lineI, double *contI, double **lineengy, double **contengy, double *lineta, double **nphcont, double **nphline, double contlabinteg, char *blrfile, double *lineinten_main, double *continten_main, int fileind)
{

  double denomsum, linedenom, contd1, expterm, etermfac, denom;
  double bphterm, lineinten[4000], continten[4000], Beta_BLF;
  int j, m;

  Beta_BLF = beta(BLF);

  for (j = 0; j < 4000; j++) {

    bphterm = 1.0 + (Beta_BLF * lineta[j]);
    etermfac = BLF * bphterm / AVTEMP;
    denom = SQR((BLF * bphterm)) * SQR((BLF * bphterm));

    blrintens(LBLR, optdep, covfac, zht, rin_blr, rout_blr, lineta[j], BLF, lineinteg, continteg, lineI, contI);
    lineinten[j] = *lineI;
    continten[j] = *contI;
    lineinten_main[j] = lineinten[j];
    continten_main[j] = continten[j];

    for (m = 0; m < LINEGRID; m++) {

      nphline[j][m] = (LINECOEFF * Nnuline[m] * lineinten[j]) /
	(lineengy[j][m] * denom);

      if (nphline[j][m] < 1.0e-100) {
        nphline[j][m] = 0.0;
      }

    }

    for (m = 0; m < CONTGRID; m++) {

      expterm = contengy[j][m] * etermfac;

      if (expterm >= 150.0) {
        nphcont[j][m] = 0.0;
        continue;

      } else if ((expterm >= 24.0) && (expterm < 150.0)) {
        contd1 = exp(expterm);
        nphcont[j][m] = (CONTCOEFF * continten[j] * SQR(contengy[j][m]))/
          (contd1 * contlabinteg);

      } else if ((expterm < 24.0) && (expterm > 1.0e-4)) {
        contd1 = exp(expterm) - 1.0;
        nphcont[j][m] = (CONTCOEFF * continten[j] * SQR(contengy[j][m]))/
          (contd1 * contlabinteg);

      } else {
        contd1 = expterm;
        nphcont[j][m] = (CONTCOEFF * continten[j] * SQR(contengy[j][m]))/
          (contd1 * contlabinteg);

      }

      if (nphcont[j][m] < 1.0e-100) {
        nphcont[j][m] = 0.0;
      }

    }

  }

  return;

}


// Start of nph calculation resulting from the BLR intensity profile,
// calculated in the plasma frame (PF). The formula has been taken from Liu &
// Bai 2006.
//
void nphblrtot(double *Nnuline, double BLF, double lineinteg, double continteg, double LBLR, double optdep, double covfac, double zht, double rin_blr, double rout_blr, double etagausmin, double etagausmax, double etalinmin, double etalinmax, double *lineI, double *contI, double lineengyp[][LINEGRID], double contengyp[][CONTGRID], double lineengym[][LINEGRID], double contengym[][CONTGRID], double lineengylin[][LINEGRID], double contengylin[][CONTGRID], double *lineta, double *gausetap, double *gausetam, double *x, double *w, double nphcontlin[][CONTGRID], double nphlinelin[][LINEGRID], double nphcontgausp[][CONTGRID], double nphlinegausp[][LINEGRID], double nphcontgausm[][CONTGRID], double nphlinegausm[][LINEGRID], double contlabinteg, char *blrfile, double *lineinten_linmain, double *continten_linmain, double *lineinten_gauspmain, double *continten_gauspmain, double *lineinten_gausmmain, double *continten_gausmmain, int fileind, int GAUSASZ)
{

  /*
  FILE *fp_gausinten, *fp_lininten, *fp_funnphgausp, *fp_funnphgausm;
  FILE *fp_funnphgaustot, *fp_funnphlin, *fp_funnphlinegausp;
  FILE *fp_funnphlinegausm, *fp_funnphcontgausp, *fp_funnphcontgausm;
  FILE *fp_funnphlinelin, *fp_funnphcontlin;
  */

  double denomsump, linedenomp, contd1p, exptermp, etermfacp, denomp;
  double bphtermp, lineinten_gausp[GAUSASZ], continten_gausp[GAUSASZ];
  double denomsumm, linedenomm, contd1m, exptermm, etermfacm, denomm;
  double bphtermm, lineinten_gausm[GAUSASZ], continten_gausm[GAUSASZ];
  double lineinten_lin[ETASTEP], continten_lin[ETASTEP], Beta_BLF, bphterm;
  double denomsum, linedenom, contd1, expterm, etermfac, denom;
  int j, m;
  /*
  char* filenamelin = (char*) malloc(sizeof(char) * 1024);
  char* filenamegaus = (char*) malloc(sizeof(char) * 1024);
  char* filenphgausp = (char*) malloc(sizeof(char) * 1024);
  char* filenphgausm = (char*) malloc(sizeof(char) * 1024);
  char* filenphlin = (char*) malloc(sizeof(char) * 1024);
  char* filenphgaustot = (char*) malloc(sizeof(char) * 1024);
  char* filenphlinegausp = (char*) malloc(sizeof(char) * 1024);
  char* filenphlinegausm = (char*) malloc(sizeof(char) * 1024);
  char* filenphcontgausp = (char*) malloc(sizeof(char) * 1024);
  char* filenphcontgausm = (char*) malloc(sizeof(char) * 1024);
  char* filenphlinelin = (char*) malloc(sizeof(char) * 1024);
  char* filenphcontlin = (char*) malloc(sizeof(char) * 1024);
  */
  Beta_BLF = beta(BLF);

  /*
  sprintf(filenamegaus, "%s%d%s%s", blrfile, fileind, "_gaus", ".dat");
  fp_gausinten = fopen(filenamegaus, "a");
  if(!fp_gausinten) {

    printf("\nCould not open the output file blr_Iz_gaus. \n");
    return;
  }

  sprintf(filenamelin, "%s%d%s%s", blrfile, fileind, "_lin", ".dat");
  fp_lininten = fopen(filenamelin, "a");
  if(!fp_lininten) {

    printf("\nCould not open the output file blr_Iz_lin. \n");
    return;
  }
  */

  for (j = 0; j < GAUSASZ; j++) {

    /*
    sprintf(filenphlinegausp, "%s%d%s%d%s", blrfile, fileind, "_nphlinegausp", j, ".dat");
    fp_funnphlinegausp = fopen(filenphlinegausp, "w");
    if(!fp_funnphlinegausp) {

      printf("\nCould not open the output file blr_nph_linegausp. \n");
      return;
    }

    sprintf(filenphlinegausm, "%s%d%s%d%s", blrfile, fileind, "_nphlinegausm", j, ".dat");
    fp_funnphlinegausm = fopen(filenphlinegausm, "w");
    if(!fp_funnphlinegausm) {

      printf("\nCould not open the output file blr_nph_linegausm. \n");
      return;
    }

    sprintf(filenphcontgausp, "%s%d%s%d%s", blrfile, fileind, "_nphcontgausp", j, ".dat");
    fp_funnphcontgausp = fopen(filenphcontgausp, "w");
    if(!fp_funnphcontgausp) {

      printf("\nCould not open the output file blr_nph_contgausp. \n");
      return;
    }

    sprintf(filenphcontgausm, "%s%d%s%d%s", blrfile, fileind, "_nphcontgausm", j, ".dat");
    fp_funnphcontgausm = fopen(filenphcontgausm, "w");
    if(!fp_funnphcontgausm) {

      printf("\nCould not open the output file blr_nph_contgausm. \n");
      return;
    }
    */

    bphtermp = 1.0 + (Beta_BLF * gausetap[j]);
    etermfacp = BLF * bphtermp / AVTEMP;
    denomp = SQR((BLF * bphtermp)) * SQR((BLF * bphtermp));

    blrintens(LBLR, optdep, covfac, zht, rin_blr, rout_blr, gausetap[j], BLF, lineinteg, continteg, lineI, contI);
    lineinten_gausp[j] = *lineI;
    continten_gausp[j] = *contI;
    lineinten_gauspmain[j] = lineinten_gausp[j];
    continten_gauspmain[j] = continten_gausp[j];

    /*
    fprintf(fp_gausinten, "%e %e %e\n", gausetap[j], lineinten_gausp[j],
	    continten_gausp[j]);
    */
    for (m = 0; m < LINEGRID; m++) {

      nphlinegausp[j][m] = (LINECOEFF * Nnuline[m] * lineinten_gausp[j]) /
        (lineengyp[j][m] * denomp);

      //fprintf(fp_funnphlinegausp, "%e %e\n", lineengyp[j][m], nphlinegausp[j][m]);

      if (nphlinegausp[j][m] < 1.0e-100) {
	nphlinegausp[j][m] = 0.0;
      }

    }

    //fclose(fp_funnphlinegausp);

    for (m = 0; m < CONTGRID; m++) {

      exptermp = contengyp[j][m] * etermfacp;

      if (exptermp >= 150.0) {
	nphcontgausp[j][m] = 0.0;
	continue;

      } else if ((exptermp >= 24.0) && (exptermp < 150.0)) {
	contd1p = exp(exptermp);
	nphcontgausp[j][m] = (CONTCOEFF * continten_gausp[j] * SQR(contengyp[j][m]))
	  / (contd1p * contlabinteg);

      } else if ((exptermp < 24.0) && (exptermp > 1.0e-4)) {
	contd1p = exp(exptermp) - 1.0;
	nphcontgausp[j][m] = (CONTCOEFF * continten_gausp[j] * SQR(contengyp[j][m]))
	  / (contd1p * contlabinteg);

      } else {
	contd1p = exptermp;
	nphcontgausp[j][m] = (CONTCOEFF * continten_gausp[j] * SQR(contengyp[j][m]))
	  / (contd1p * contlabinteg);

      }

      //fprintf(fp_funnphcontgausp, "%e %e\n", contengyp[j][m], nphcontgausp[j][m]);

      if (nphcontgausp[j][m] < 1.0e-100) {
	nphcontgausp[j][m] = 0.0;
      }

    }

    //fclose(fp_funnphcontgausp);

    // Calculating the corresponding gauseta_minus part.
    //
    bphtermm = 1.0 + (Beta_BLF * gausetam[j]);
    etermfacm = BLF * bphtermm / AVTEMP;
    denomm = SQR((BLF * bphtermm)) * SQR((BLF * bphtermm));

    blrintens(LBLR, optdep, covfac, zht, rin_blr, rout_blr, gausetam[j], BLF, lineinteg, continteg, lineI, contI);
    lineinten_gausm[j] = *lineI;
    continten_gausm[j] = *contI;
    lineinten_gausmmain[j] = lineinten_gausm[j];
    continten_gausmmain[j] = continten_gausm[j];

    /*
    fprintf(fp_gausinten, "%e %e %e\n", gausetam[j], lineinten_gausm[j],
	    continten_gausm[j]);
    */
    for (m = 0; m < LINEGRID; m++) {

      nphlinegausm[j][m] = (LINECOEFF * Nnuline[m] * lineinten_gausm[j]) /
        (lineengym[j][m] * denomm);

      //fprintf(fp_funnphlinegausm, "%e %e\n", lineengym[j][m], nphlinegausm[j][m]);

      if (nphlinegausm[j][m] < 1.0e-100) {
	nphlinegausm[j][m] = 0.0;
      }

    }

    //fclose(fp_funnphlinegausm);

    for (m = 0; m < CONTGRID; m++) {

      exptermm = contengym[j][m] * etermfacm;

      if (exptermm >= 150.0) {
	nphcontgausm[j][m] = 0.0;
	continue;

      } else if ((exptermm >= 24.0) && (exptermm < 150.0)) {
	contd1m = exp(exptermm);
	nphcontgausm[j][m] = (CONTCOEFF * continten_gausm[j] * SQR(contengym[j][m]))
	  / (contd1m * contlabinteg);

      } else if ((exptermm < 24.0) && (exptermm > 1.0e-4)) {
	contd1m = exp(exptermm) - 1.0;
	nphcontgausm[j][m] = (CONTCOEFF * continten_gausm[j] * SQR(contengym[j][m]))
	  / (contd1m * contlabinteg);

      } else {
	contd1m = exptermm;
	nphcontgausm[j][m] = (CONTCOEFF * continten_gausm[j] * SQR(contengym[j][m]))
	  / (contd1m * contlabinteg);

      }

      //fprintf(fp_funnphcontgausm, "%e %e\n", contengym[j][m], nphcontgausm[j][m]);

      if (nphcontgausm[j][m] < 1.0e-100) {
	nphcontgausm[j][m] = 0.0;
      }

    }

    //fclose(fp_funnphcontgausm);

  }

  //fclose(fp_gausinten);

  // Calculating the corresponding linear part.
  //
  for (j = 0; j < ETASTEP; j++) {

    bphterm = 1.0 + (Beta_BLF * lineta[j]);
    etermfac = BLF * bphterm / AVTEMP;
    denom = SQR((BLF * bphterm)) * SQR((BLF * bphterm));

    blrintens(LBLR, optdep, covfac, zht, rin_blr, rout_blr, lineta[j], BLF, lineinteg, continteg, lineI, contI);
    lineinten_lin[j] = *lineI;
    continten_lin[j] = *contI;
    lineinten_linmain[j] = lineinten_lin[j];
    continten_linmain[j] = continten_lin[j];

    /*
    fprintf(fp_lininten, "%e %e %e\n", lineta[j], lineinten_lin[j],
	    continten_lin[j]);

    sprintf(filenphlinelin, "%s%d%s%d%s", blrfile, fileind, "_nphlinelin", j, ".dat");
    fp_funnphlinelin = fopen(filenphlinelin, "w");
    if(!fp_funnphlinelin) {

      printf("\nCould not open the output file blr_nph_linelin. \n");
      return;
    }

    sprintf(filenphcontlin, "%s%d%s%d%s", blrfile, fileind, "_nphcontlin", j, ".dat");
    fp_funnphcontlin = fopen(filenphcontlin, "w");
    if(!fp_funnphcontlin) {

      printf("\nCould not open the output file blr_nph_contlin. \n");
      return;
    }
    */

    for (m = 0; m < LINEGRID; m++) {

      nphlinelin[j][m] = (LINECOEFF * Nnuline[m] * lineinten_lin[j]) /
        (lineengylin[j][m] * denom);

      //fprintf(fp_funnphlinelin, "%e %e\n", lineengylin[j][m], nphlinelin[j][m]);

      if (nphlinelin[j][m] < 1.0e-100) {
	nphlinelin[j][m] = 0.0;
      }

    }

    //fclose(fp_funnphlinelin);

    for (m = 0; m < CONTGRID; m++) {

      expterm = contengylin[j][m] * etermfac;

      if (expterm >= 150.0) {
	nphcontlin[j][m] = 0.0;
	continue;

      } else if ((expterm >= 24.0) && (expterm < 150.0)) {
	contd1 = exp(expterm);
	nphcontlin[j][m] = (CONTCOEFF * continten_lin[j] * SQR(contengylin[j][m]))/
	  (contd1 * contlabinteg);

      } else if ((expterm < 24.0) && (expterm > 1.0e-4)) {
	contd1 = exp(expterm) - 1.0;
	nphcontlin[j][m] = (CONTCOEFF * continten_lin[j] * SQR(contengylin[j][m]))/
	  (contd1 * contlabinteg);

      } else {
	contd1 = expterm;
	nphcontlin[j][m] = (CONTCOEFF * continten_lin[j] * SQR(contengylin[j][m]))/
	  (contd1 * contlabinteg);

      }

      //fprintf(fp_funnphcontlin, "%e %e\n", contengylin[j][m], nphcontlin[j][m]);

      if (nphcontlin[j][m] < 1.0e-100) {
	nphcontlin[j][m] = 0.0;
      }

    }

    //fclose(fp_funnphcontlin);

  }

  //fclose(fp_lininten);

  return;

}


// Start of gammadot calculation resulting from inverse Compton scattering of
// BLR photons in the emission region, calculated in the PF. The formula is
// based on head approx. given in Dermer & Menon book, 2007, Ch. 6.
//
double nphepinteg(double eengy, double phcos, double *photengy, double *photdensiy, int NGRID)
{
  double sum = 0.0, Mnph, Dnph, beta_enph, term1nph, term2nph, term3nph, num1nph;
  double den1nph, num2nph, den2nph, lnnph, D2nph, D3nph, gpenph;
  double num3nph, den3nph, funofep[NGRID], bphtermnph;
  double bpheengyconst;
  int j;

  beta_enph = beta(eengy);
  bphtermnph = 1.0 - (beta_enph * phcos);
  bpheengyconst = eengy * bphtermnph;

  gpenph = eengy  - photengy[0];
  Mnph = bpheengyconst * photengy[0];
  Dnph = 1.0 + (2.0 * Mnph);

  num1nph = (Mnph * (Mnph - 2.0)) - 2.0;
  den1nph = eengy * photengy[0] * SQR(Mnph);

  if ((2.0 * Mnph) < 1.0e-4) {

    lnnph = 2.0 * Mnph;
    D3nph = 1.0 + (6.0 * Mnph);
    D2nph = 1.0 + (4.0 * Mnph);

  } else {

    lnnph = log(Dnph);
    D3nph = pow(Dnph, 3.0);
    D2nph = SQR(Dnph);

  }

  term1nph = lnnph * ((gpenph * num1nph) - eengy) / den1nph;

  den2nph = 3.0 * photengy[0] * D3nph;
  num2nph = 6.0 * photengy[0] * Dnph * gpenph * (1.0 + Mnph) * bphtermnph;
  term2nph = 1.0 - D3nph + num2nph;

  num3nph = (2.0 * gpenph * Dnph) - (eengy * (Mnph * (Mnph - 1.0) - 1.0));
  den3nph = eengy * Mnph;
  term3nph = 6.0 * D2nph * num3nph / den3nph;

  funofep[0] = photdensiy[0] * (term1nph + ((term2nph + term3nph) / den2nph));

  if (funofep[0] < 0.0) {
    funofep[0] = 0.0;
  }

  for (j = 1; j < NGRID; j++) {

    gpenph = eengy  - photengy[j];
    Mnph = bpheengyconst * photengy[j];
    Dnph = 1.0 + (2.0 * Mnph);

    num1nph = (Mnph * (Mnph - 2.0)) - 2.0;
    den1nph = eengy * photengy[j] * SQR(Mnph);

    if ((2.0 * Mnph) < 1.0e-4) {

      lnnph = 2.0 * Mnph;
      D3nph = 1.0 + (6.0 * Mnph);
      D2nph = 1.0 + (4.0 * Mnph);

    } else {

      lnnph = log(Dnph);
      D3nph = pow(Dnph, 3.0);
      D2nph = SQR(Dnph);

    }

    term1nph = lnnph * ((gpenph * num1nph) - eengy) / den1nph;

    den2nph = 3.0 * photengy[j] * D3nph;
    num2nph = 6.0 * photengy[j] * Dnph * gpenph * (1.0 + Mnph) * bphtermnph;
    term2nph = 1.0 - D3nph + num2nph;

    num3nph = (2.0 * gpenph * Dnph) - (eengy * (Mnph * (Mnph - 1.0) - 1.0));
    den3nph = eengy * Mnph;
    term3nph = 6.0 * D2nph * num3nph / den3nph;

    funofep[j] = photdensiy[j] * (term1nph + ((term2nph + term3nph) / den2nph));

    if (funofep[j] < 0.0) {
      funofep[j] = 0.0;
      continue;
    }

    sum += (0.5 * (funofep[j] + funofep[j - 1]) *
	    (photengy[j] - photengy[j - 1]));

  }

  return(sum);

}


double nphepsumline(double eengy, double phcos, double *photengy, double *photdenblr)
{
  double sum = 0.0, term1, term2, M, D, num1, den1, den2, num2, beta_e;
  double den3, num3, term3, funofepline[LINEGRID], bphterm, gpe;
  double lnterm, D2term, D3term;
  int j;

  beta_e = beta(eengy);
  bphterm = 1.0 - (beta_e * phcos);

  for (j = 0; j < LINEGRID; j++) {

    gpe = eengy  - photengy[j];
    M = eengy * photengy[j] * bphterm;
    D = 1.0 + (2.0 * M);

    num1 = (M * (M - 2.0)) - 2.0;
    den1 = eengy * photengy[j] * SQR(M);

    if ((2.0 * M) < 1.0e-4) {

      lnterm = 2.0 * M;
      D3term = 1.0 + (6.0 * M);
      D2term = 1.0 + (4.0 * M);

    } else {

      lnterm = log(D);
      D3term = pow(D, 3.0);
      D2term = SQR(D);

    }

    term1 = lnterm * ((gpe * num1) - eengy) / den1;

    den2 = 3.0 * photengy[j] * D3term;
    num2 = 6.0 * photengy[j] * D * gpe * (1.0 + M) * bphterm;
    term2 = 1.0 - D3term + num2;

    num3 = (2.0 * gpe * D) - (eengy * (M * (M - 1.0) - 1.0));
    den3 = eengy * M;
    term3 = 6.0 * D2term * num3 / den3;

    funofepline[j] = photdenblr[j] * (term1 + ((term2 + term3) / den2));

    if (funofepline[j] < 0.0) {
      funofepline[j] = 0.0;
      continue;
    }

    sum += funofepline[j];

  }

  return(sum);

}


double gamdotblr_gaus(double etmin, double etmax, double eengy, double photengy_linep[][LINEGRID], double photengy_linem[][LINEGRID], double photengy_contp[][CONTGRID], double photengy_contm[][CONTGRID], double pdline_gausblrp[][LINEGRID], double pdline_gausblrm[][LINEGRID], double pdcont_gausblrp[][CONTGRID], double pdcont_gausblrm[][CONTGRID], double *w, double *gausetap, double *gausetam, int GAUSASZ)
{
  double sumcont = 0.0, sumline = 0.0, sum;
  double term1, term2, xr, term1line, term2line;
  int j, k = 0;

  xr = 0.5 * (etmax - etmin);

  for (j = 1; j <= GAUSASZ; j++) {

    term1 = nphepinteg(eengy, gausetap[k], photengy_contp[k], pdcont_gausblrp[k], CONTGRID);
    term1line = nphepsumline(eengy, gausetap[k], photengy_linep[k], pdline_gausblrp[k]);

    term2 = nphepinteg(eengy, gausetam[k], photengy_contm[k], pdcont_gausblrm[k], CONTGRID);
    term2line = nphepsumline(eengy, gausetam[k], photengy_linem[k], pdline_gausblrm[k]);

    sumcont += (w[j] * (term1 + term2));
    sumline += (w[j] * (term1line + term2line));

    k++;

  }

  sum = sumcont + sumline;

  return(sum *= xr);

}


double gamdotblr_lin(double eengy, double photengy_line[][LINEGRID], double photengy_cont[][CONTGRID], double pdline_linblr[][LINEGRID], double pdcont_linblr[][CONTGRID], double *lineta)
{
  double sumcont = 0.0, sumline = 0.0, yofline_lineta[ETASTEP];
  double sum, yofcont_lineta[ETASTEP], beta_e, bphterm, lineterm;
  int j;

  yofcont_lineta[0] = nphepinteg(eengy, lineta[0], photengy_cont[0], pdcont_linblr[0], CONTGRID);

  if (yofcont_lineta[0] < 0.0) {
    yofcont_lineta[0] = 0.0;
  }

  yofline_lineta[0] = nphepsumline(eengy, lineta[0], photengy_line[0], pdline_linblr[0]);

  if (yofline_lineta[0] < 0.0) {
    yofline_lineta[0] = 0.0;
  }

  for (j = 1; j < ETASTEP; j++) {

    yofcont_lineta[j] = nphepinteg(eengy, lineta[j], photengy_cont[j], pdcont_linblr[j], CONTGRID);

    if (yofcont_lineta[j] < 0.0) {
      yofcont_lineta[j] = 0.0;
      continue;
    }

    sumcont += (0.5 * (yofcont_lineta[j] + yofcont_lineta[j - 1]) *
		(lineta[j] - lineta[j - 1]));

    yofline_lineta[j] = nphepsumline(eengy, lineta[j], photengy_line[j], pdline_linblr[j]);

    if (yofline_lineta[j] < 0.0) {
      yofline_lineta[j] = 0.0;
      continue;
    }

    sumline += (0.5 * (yofline_lineta[j] + yofline_lineta[j - 1]) *
                (lineta[j] - lineta[j - 1]));

  }

  sum = sumline + sumcont;

  return(sum);

}


double gamdotblr(double eengy, double gausetmin, double gausetmax, double pe_linp[][LINEGRID], double pe_linm[][LINEGRID], double pe_contp[][CONTGRID], double pe_contm[][CONTGRID], double pe_linelin[][LINEGRID], double pe_contlin[][CONTGRID], double photdenline_linblr[][LINEGRID], double photdencont_linblr[][CONTGRID], double *lineta, double photdenline_gausblrp[][LINEGRID], double photdenline_gausblrm[][LINEGRID], double *w, double *gausetap, double *gausetam, double photdencont_gausblrp[][CONTGRID], double photdencont_gausblrm[][CONTGRID], int GAUSASZ)
{

  double gammadot, gam1, gam2;
  int j;

  gam1 = gamdotblr_gaus(gausetmin, gausetmax, eengy, pe_linp, pe_linm, pe_contp, pe_contm, photdenline_gausblrp, photdenline_gausblrm, photdencont_gausblrp, photdencont_gausblrm, w, gausetap, gausetam, GAUSASZ);
  gam2 = gamdotblr_lin(eengy, pe_linelin, pe_contlin, photdenline_linblr, photdencont_linblr, lineta);

  gammadot = GDOTBLRCOEFF * (gam1 + gam2);

  return(-gammadot);

}


// Start of ECBLR ndot, in PF, calculation using Compton cross-section under
// head-on approx. based on eqn. derived using Dermer & Menon book, Ch-6.
//
double Int_gamblr(double coski, double scatterphotengy, double photengy, double *e_engy, double *nep)
{
  double sum = 0.0, yofgam[EGRID], cross_sec, gamLL, gammid, term1;
  double gamepsterm, epepsterm, coskiterm, sqrtterm;
  int j, n = 0;

  coskiterm = 1.0 - coski;
  term1 = 2.0 / (photengy * scatterphotengy * coskiterm);

  if (term1 < 1.0e-4) {    //Some small number to enable binomial expansion.
    sqrtterm = 1.0 + (0.5 * term1);

  } else {
    sqrtterm = sqrt(1.0 + term1);
  }

  gamLL = 0.5 * scatterphotengy * (1.0 + sqrtterm);
  gamLL = max(gamLL, e_engy[0]);

  while ((n < EGRID) && (e_engy[n] < gamLL)) {
    n++;
  }

  if (((n < EGRID) && (nep[n] > 0.0)) && (gamLL > e_engy[0])) {

    gammid = 0.5 * (e_engy[n] + gamLL);

    if (gammid >= (0.1 / photengy)) {

      gamepsterm = (gammid - scatterphotengy) / gammid;
      epepsterm = scatterphotengy / (photengy * gammid * coskiterm *
				     (gammid - scatterphotengy));

      cross_sec = gamepsterm + (1.0 / gamepsterm) - (2.0 * epepsterm) +
	SQR(epepsterm);

    } else {

      epepsterm = scatterphotengy / (SQR(gammid) * photengy * coskiterm);
      cross_sec = 2.0 - (2.0 * epepsterm) + SQR(epepsterm);

    }

    yofgam[n] = (nep[n] * cross_sec) / SQR(gammid);

    sum += ((e_engy[n] - gamLL) * yofgam[n]);

    if (yofgam[n] < 0.0) {
      sum = 0.0;
    }

    n++;

  }

  if (n < EGRID) {

    if (nep[n] == 0.0) {

      yofgam[n] = 0.0;

    } else {

      if (e_engy[n] >= (0.1 / photengy)) {

	gamepsterm = (e_engy[n] - scatterphotengy) / e_engy[n];
	epepsterm = scatterphotengy / (photengy * e_engy[n] * coskiterm *
				       (e_engy[n] - scatterphotengy));

	cross_sec = gamepsterm + (1.0 / gamepsterm) - (2.0 * epepsterm) +
	  SQR(epepsterm);

      } else {

	epepsterm = scatterphotengy / (SQR(e_engy[n]) * photengy * coskiterm);
	cross_sec = 2.0 - (2.0 * epepsterm) + SQR(epepsterm);

      }

      yofgam[n] = (nep[n] * cross_sec) / SQR(e_engy[n]);

      if (yofgam[n] < 0.0) {
	yofgam[n] = 0.0;
      }

    }

    for (j = (n + 1); j < EGRID; j++) {

      if (nep[j] == 0.0) {

	yofgam[j] = 0.0;
	continue;
      }

      if (e_engy[j] >= (0.1 / photengy)) {

	gamepsterm = (e_engy[j] - scatterphotengy) / e_engy[j];
	epepsterm = scatterphotengy / (photengy * e_engy[j] * coskiterm *
				       (e_engy[j] - scatterphotengy));
	
	cross_sec = gamepsterm + (1.0 / gamepsterm) - (2.0 * epepsterm) +
	  SQR(epepsterm);

      } else {

        epepsterm = scatterphotengy / (SQR(e_engy[j]) * photengy * coskiterm);
	cross_sec = 2.0 - (2.0 * epepsterm) + SQR(epepsterm);

      }

      yofgam[j] = (nep[j] * cross_sec) / SQR(e_engy[j]);

      if (yofgam[j] < 0.0) {
	yofgam[j] = 0.0;
	continue;
      }

      sum += (0.5 * (yofgam[j] + yofgam[j - 1]) * (e_engy[j] - e_engy[j - 1]));

    }

  }

  return (sum);

}


double Int_epblr(double coski, double scatterp_engy, double epsmax, double *e_engy, double *nep, double *p_engy, double *np_blr)
{
  double sum = 0.0, yofep[CONTGRID], epUL, epmid;
  int j, n = 0;

  epUL = (2.0 * scatterp_engy) / (1.0 - coski);
  epUL = min(epUL, epsmax);

  if (epUL < epsmax) {

    while ((n < CONTGRID) && (p_engy[n] < epUL)) {
      n++;

    }

  } else {

    n = CONTGRID - 1;

  }

  yofep[0] = (np_blr[0] / p_engy[0]) *
    Int_gamblr(coski, scatterp_engy, p_engy[0], e_engy, nep);

  if (yofep[0] < 1.0e-100) yofep[0] = 0.0;

  for (j = 1; j < CONTGRID; j++) {

    if (p_engy[j] > epUL) {
      break;
    }

    yofep[j] = (np_blr[j] / p_engy[j]) *
      Int_gamblr(coski, scatterp_engy, p_engy[j], e_engy, nep);

    if (yofep[j] < 1.0e-100) {
      yofep[j] = 0.0;
      continue;
    }

    sum += (0.5 * (yofep[j] + yofep[j - 1]) * (p_engy[j] - p_engy[j - 1]));

  }

  return (sum);

}


double Int_epblrline(double coski, double scatterp_engy, double epsmax, double *e_engy, double *nep, double *p_engy, double *np_blr)
{
  double sum = 0.0, yofep[LINEGRID], epUL, epmid;
  int j, n = 0;

  epUL = (2.0 * scatterp_engy) / (1.0 - coski);
  epUL = min(epUL, epsmax);

  if (epUL < epsmax) {

    while ((n < LINEGRID) && (p_engy[n] < epUL)) {
      n++;

    }

  } else {

    n = LINEGRID - 1;

  }

  for (j = 0; j < LINEGRID; j++) {

    if (p_engy[j] > epUL) {
      break;
    }

    yofep[j] = (np_blr[j] / p_engy[j]) *
      Int_gamblr(coski, scatterp_engy, p_engy[j], e_engy, nep);

    if (yofep[j] < 1.0e-100) {
      yofep[j] = 0.0;
      continue;
    }

    sum += yofep[j];

  }

  return (sum);

}


double etablr_integ(double eps, double gausmin, double gausmax, double *e_engy, double *nep, double pe_linep[][LINEGRID], double pe_linem[][LINEGRID], double pe_contp[][CONTGRID], double pe_contm[][CONTGRID], double pe_linelin[][LINEGRID], double pe_contlin[][CONTGRID], double pdline_gausblrp[][LINEGRID], double pdline_gausblrm[][LINEGRID], double pdcont_gausblrp[][CONTGRID], double pdcont_gausblrm[][CONTGRID], double *w, double pdline_linblr[][LINEGRID], double pdcont_linblr[][CONTGRID], double *lineta, int GAUSASZ, double coskip[], double coskim[], double coskilin[])
{
  double sumgaus = 0.0, yofline_lineta[ETASTEP], yofcont_lineta[ETASTEP];
  double Beta_bulk, term1, term2, term1line, term2line;
  double sumlin = 0.0, sumgausline = 0.0, sumgauscont = 0.0, sumlinline = 0.0;
  double sumlincont = 0.0, sum = 0.0, xr;

  int j, k = 0;

  // For the gaus part of etaph
  //
  xr = 0.5 * (gausmax - gausmin);

  for (j = 1; j <= GAUSASZ; j++) {

    term1 = Int_epblr(coskip[k], eps, pe_contp[k][CONTGRID-1], e_engy, nep, pe_contp[k], pdcont_gausblrp[k]);
    term1line = Int_epblrline(coskip[k], eps, pe_linep[k][LINEGRID-1], e_engy, nep, pe_linep[k], pdline_gausblrp[k]);

    term2 = Int_epblr(coskim[k], eps, pe_contm[k][CONTGRID-1], e_engy, nep, pe_contm[k], pdcont_gausblrm[k]);
    term2line = Int_epblrline(coskim[k], eps, pe_linem[k][LINEGRID-1], e_engy, nep, pe_linem[k], pdline_gausblrm[k]);

    sumgauscont += (w[j] * (term1 + term2));
    sumgausline += (w[j] * (term1line + term2line));

    k++;

  }

  sumgaus = sumgauscont + sumgausline;
  sumgaus *= xr;

  // For the linear part of etaph
  //
  yofcont_lineta[0] = Int_epblr(coskilin[0], eps, pe_contlin[0][CONTGRID-1], e_engy, nep, pe_contlin[0], pdcont_linblr[0]);
  yofline_lineta[0] = Int_epblrline(coskilin[0], eps, pe_linelin[0][LINEGRID-1], e_engy, nep, pe_linelin[0], pdline_linblr[0]);

  for (j = 1; j < ETASTEP; j++) {

    yofcont_lineta[j] = Int_epblr(coskilin[j], eps, pe_contlin[j][CONTGRID-1], e_engy, nep, pe_contlin[j], pdcont_linblr[j]);
    yofline_lineta[j] = Int_epblrline(coskilin[j], eps, pe_linelin[j][LINEGRID-1], e_engy, nep, pe_linelin[j], pdline_linblr[j]);

    sumlincont += (0.5 * (yofcont_lineta[j] + yofcont_lineta[j - 1]) *
		   (lineta[j] - lineta[j - 1]));

    sumlinline += (0.5 * (yofline_lineta[j] + yofline_lineta[j - 1]) *
		   (lineta[j] - lineta[j - 1]));

  }

  sumlin = sumlincont + sumlinline;

  sum = sumgaus + sumlin;

  return(sum);

}


double ndotblr(double eps, double muobsprime, double *e_engy, double *nep, double gausetmin, double gausetmax, double pe_linp[][LINEGRID], double pe_linm[][LINEGRID], double pe_contp[][CONTGRID], double pe_contm[][CONTGRID], double pe_linelin[][LINEGRID], double pe_contlin[][CONTGRID], double photdenline_linblr[][LINEGRID], double photdencont_linblr[][CONTGRID], double *lineta, double photdenline_gausblrp[][LINEGRID], double photdenline_gausblrm[][LINEGRID], double *w, double *gausetap, double *gausetam, double photdencont_gausblrp[][CONTGRID], double photdencont_gausblrm[][CONTGRID], int GAUSASZ, double coskip[][GAUSASZ], double coskim[][GAUSASZ], double coskilin[][ETASTEP], double phistep)
{
  double finalterm, sum = 0.0, yofphi[PHIGRID];
  int i, j;

  if (!(fabs(muobsprime) < .99999999)) {
    printf ("\n muobsprime = %e !! ", muobsprime);
    return (0.0);
  }

  yofphi[0] = etablr_integ(eps, gausetmin, gausetmax, e_engy, nep, pe_linp, pe_linm, pe_contp, pe_contm, pe_linelin, pe_contlin, photdenline_gausblrp, photdenline_gausblrm, photdencont_gausblrp, photdencont_gausblrm, w, photdenline_linblr, photdencont_linblr, lineta, GAUSASZ, coskip[0], coskim[0], coskilin[0]);

  for (j = 1; j < PHIGRID; j++) {

    yofphi[j] = etablr_integ(eps, gausetmin, gausetmax, e_engy, nep, pe_linp, pe_linm, pe_contp, pe_contm, pe_linelin, pe_contlin, photdenline_gausblrp, photdenline_gausblrm, photdencont_gausblrp, photdencont_gausblrm, w, photdenline_linblr, photdencont_linblr, lineta, GAUSASZ, coskip[j], coskim[j], coskilin[j]);

    sum += (0.5 * (yofphi[j] + yofphi[j - 1]) * phistep);

  }


  finalterm = SSCOEFF * sum;

  return(2.0 * finalterm);

}


// Start of gammadot calculation resulting from inverse Compton scattering of
// BLR photons in the emission region, in the PF, in Thomson regime.
//
double gamdotblr_gausth(double etmin, double etmax, double eengy, double BLF, double contlabinteg, double *inline_gausblrp, double *inline_gausblrm, double *incont_gausblrp, double *incont_gausblrm, double *w, double *gausetap, double *gausetam, int GAUSASZ)
{
  double sum = 0.0, pdline_gausp[GAUSASZ], pdcont_gausp[GAUSASZ], xr;
  double pdline_gausm[GAUSASZ], yofgausetap[GAUSASZ], Beta_BLF, beta_e;
  double bphtermp, term1p, bphtermm, term1m, yofgausetam[GAUSASZ], BLFtermp;
  double pdcont_gausm[GAUSASZ], denomp, denomm, BLFtermm;
  int j, k = 0;

  xr = 0.5 * (etmax - etmin);

  Beta_BLF = beta(BLF);
  beta_e = beta(eengy);

  for (j = 1; j <= GAUSASZ; j++) {

    BLFtermp = 1.0 + (Beta_BLF * gausetap[k]);
    denomp = SQR(BLF * BLFtermp) * SQR(BLF * BLFtermp);

    bphtermp = 1.0 - (beta_e * gausetap[k]);
    term1p = bphtermp * ((SQR(eengy) * bphtermp) - 1.0);

    pdline_gausp[k] = inline_gausblrp[k];
    pdcont_gausp[k] = THCONTCOEFF * incont_gausblrp[k] / contlabinteg;

    BLFtermm = 1.0 + (Beta_BLF * gausetam[k]);
    denomm = SQR(BLF * BLFtermm) * SQR(BLF * BLFtermm);

    bphtermm = 1.0 - (beta_e * gausetam[k]);
    term1m = bphtermm * ((SQR(eengy) * bphtermm) - 1.0);

    pdline_gausm[k] = inline_gausblrm[k];
    pdcont_gausm[k] = THCONTCOEFF * incont_gausblrm[k] / contlabinteg;

    yofgausetap[k] = (pdline_gausp[k] + pdcont_gausp[k]) * term1p / denomp;

    if (yofgausetap[k] < 0.0) {
      yofgausetap[k] = 0.0;
      continue;
    }

    yofgausetam[k] = (pdline_gausm[k] + pdcont_gausm[k]) * term1m / denomm;

    if (yofgausetam[k] < 0.0) {
      yofgausetam[k] = 0.0;
      continue;
    }

    sum += (w[j] * (yofgausetap[k] + yofgausetam[k]));

    k++;

  }

  return(sum *= xr);

}


double gamdotblr_linth(double eengy, double BLF, double contlabinteg, double *inline_linblr, double *incont_linblr, double *lineta)
{

  double yoflineta[ETASTEP], pdcont[ETASTEP], sum = 0.0, beta_e, bphterm;
  double pdline[ETASTEP], BLFterm, Beta_BLF, denom, term1;
  int j;

  Beta_BLF = beta(BLF);
  beta_e = beta(eengy);

  BLFterm = 1.0 + (Beta_BLF * lineta[0]);
  denom = SQR(BLF * BLFterm) * SQR(BLF * BLFterm);

  bphterm = 1.0 - (beta_e * lineta[0]);
  term1 = bphterm * ((SQR(eengy) * bphterm) - 1.0);

  pdline[0] = inline_linblr[0];
  pdcont[0] = THCONTCOEFF * incont_linblr[0] / contlabinteg;

  yoflineta[0] = (pdline[0] + pdcont[0]) * term1 / denom;

  if (yoflineta[0] < 0.0) {
    yoflineta[0] = 0.0;
  }

  for (j = 1; j < ETASTEP; j++) {

    BLFterm = 1.0 + (Beta_BLF * lineta[j]);
    denom = SQR(BLF * BLFterm) * SQR(BLF * BLFterm);

    bphterm = 1.0 - (beta_e * lineta[j]);
    term1 = bphterm * ((SQR(eengy) * bphterm) - 1.0);

    pdline[j] = inline_linblr[j];
    pdcont[j] = THCONTCOEFF * incont_linblr[j] / contlabinteg;

    yoflineta[j] = (pdline[j] + pdcont[j]) * term1 / denom;

    if (yoflineta[j] < 0.0) {
      yoflineta[j] = 0.0;
      continue;
    }

    sum += (0.5 * (yoflineta[j] + yoflineta[j - 1]) *
		(lineta[j] - lineta[j - 1]));

  }

  return(sum);

}


double gamdotblrth(double eengy, double gausetmin, double gausetmax, double BLF, double contlabinteg, double *intenline_linblr, double *intencont_linblr, double *lineta, double *intenline_gausblrp, double *intenline_gausblrm, double *w, double *gausetap, double *gausetam, double *intencont_gausblrp, double *intencont_gausblrm, int GAUSASZ)
{

  double gamth1, gamth2, gammadotth;
  int i;

  gamth1 = gamdotblr_gausth(gausetmin, gausetmax, eengy, BLF, contlabinteg, intenline_gausblrp, intenline_gausblrm, intencont_gausblrp, intencont_gausblrm, w, gausetap, gausetam, GAUSASZ);

  gamth2 = gamdotblr_linth(eengy, BLF, contlabinteg, intenline_linblr, intencont_linblr, lineta);

  gammadotth = GDOTBLRTHCOEFF * CONTCOEFF * (gamth1 + gamth2);

  return(-gammadotth);

}


// Start of approximate gammadot calculation in Thomson regime for comparison 
// with the actual expression
//
double gamdotapprox_gausth(double etmin, double etmax, double BLF, double contlabinteg, double *inline_gausblrp, double *inline_gausblrm, double *incont_gausblrp, double *incont_gausblrm, double *w, double *gausetap, double *gausetam, int GAUSASZ)
{
  double sum = 0.0, pdline_gausp[GAUSASZ], pdcont_gausp[GAUSASZ], xr;
  double pdline_gausm[GAUSASZ], yofgausetap[GAUSASZ], Beta_BLF;
  double yofgausetam[GAUSASZ], BLFtermp;
  double pdcont_gausm[GAUSASZ], denomp, denomm, BLFtermm;
  int j, k = 0;

  xr = 0.5 * (etmax - etmin);

  Beta_BLF = beta(BLF);

  for (j = 1; j <= GAUSASZ; j++) {

    BLFtermp = 1.0 + (Beta_BLF * gausetap[k]);
    denomp = SQR(BLF * BLFtermp) * SQR(BLF * BLFtermp);

    pdline_gausp[k] = inline_gausblrp[k];
    pdcont_gausp[k] = THCONTCOEFF * incont_gausblrp[k] / contlabinteg;

    BLFtermm = 1.0 + (Beta_BLF * gausetam[k]);
    denomm = SQR(BLF * BLFtermm) * SQR(BLF * BLFtermm);

    pdline_gausm[k] = inline_gausblrm[k];
    pdcont_gausm[k] = THCONTCOEFF * incont_gausblrm[k] / contlabinteg;

    yofgausetap[k] = (pdline_gausp[k] + pdcont_gausp[k]) / denomp;

    if (yofgausetap[k] < 0.0) {
      yofgausetap[k] = 0.0;
      continue;
    }

    yofgausetam[k] = (pdline_gausm[k] + pdcont_gausm[k]) / denomm;

    if (yofgausetam[k] < 0.0) {
      yofgausetam[k] = 0.0;
      continue;
    }

    sum += (w[j] * (yofgausetap[k] + yofgausetam[k]));

    k++;

  }

  return(sum *= xr);

}


double gamdotapprox_linth(double BLF, double contlabinteg, double *inline_linblr, double *incont_linblr, double *lineta)
{
  double yoflineta[ETASTEP], pdcont[ETASTEP], sum = 0.0;
  double pdline[ETASTEP], BLFterm, Beta_BLF, denom;
  int j;

  Beta_BLF = beta(BLF);
  BLFterm = 1.0 + (Beta_BLF * lineta[0]);
  denom = SQR(BLF * BLFterm) * SQR(BLF * BLFterm);

  pdline[0] = inline_linblr[0];
  pdcont[0] = THCONTCOEFF * incont_linblr[0] / contlabinteg;

  yoflineta[0] = (pdline[0] + pdcont[0]) / denom;

  if (yoflineta[0] < 0.0) {
    yoflineta[0] = 0.0;
  }

  for (j = 1; j < ETASTEP; j++) {

    BLFterm = 1.0 + (Beta_BLF * lineta[j]);
    denom = SQR(BLF * BLFterm) * SQR(BLF * BLFterm);

    pdline[j] = inline_linblr[j];
    pdcont[j] = THCONTCOEFF * incont_linblr[j] / contlabinteg;

    yoflineta[j] = (pdline[j] + pdcont[j]) / denom;

    if (yoflineta[j] < 0.0) {
      yoflineta[j] = 0.0;
      continue;
    }

    sum += (0.5 * (yoflineta[j] + yoflineta[j - 1]) *
	    (lineta[j] - lineta[j - 1]));

  }

  return(sum);

}


double gamdotthapprox(double eengy, double gausetmin, double gausetmax, double BLF, double contlabinteg, double *intenline_linblr, double *intencont_linblr, double *lineta, double *intenline_gausblrp, double *intenline_gausblrm, double *w, double *gausetap, double *gausetam, double *intencont_gausblrp, double *intencont_gausblrm, int GAUSASZ)
{

  double gamth1, gamth2, gammadotth;

  gamth1 = gamdotapprox_gausth(gausetmin, gausetmax, BLF, contlabinteg, intenline_gausblrp, intenline_gausblrm, intencont_gausblrp, intencont_gausblrm, w, gausetap, gausetam, GAUSASZ);

  gamth2 = gamdotapprox_linth(BLF, contlabinteg, intenline_linblr, intencont_linblr, lineta);

  gammadotth = GDOTTHCOMPCOEFF * CONTCOEFF * SQR(eengy) * (gamth1 + gamth2);

  return(-gammadotth);

}



// Start of ndot calculation resulting from inverse Compton scattering of
// BLR photons in the emission region, in the PF, in Thomson regime.
//
double Int_gamblrth(double coski, double etermfactor, double *e_engy, double *nep)
{
  double sum = 0.0, yofgam[EGRID], beta_e, contd1, term1;
  double gamterm, epepsterm, coskiterm, expterm, epsterm;
  int m;

  beta_e = beta(e_engy[0]);
  coskiterm = 1.0 - (beta_e * coski);
  gamterm = SQR(e_engy[0]) * coskiterm;

  expterm = etermfactor / gamterm;

  if (expterm >= 150.0) {
    yofgam[0] = 0.0;

  } else if ((expterm >= 24.0) && (expterm < 150.0)) {
    contd1 = exp(expterm);
    yofgam[0] = nep[0] / (pow(e_engy[0], 6.0) * SQR(coskiterm) * contd1);

  } else if ((expterm < 24.0) && (expterm > 1.0e-4)) {
    contd1 = exp(expterm) - 1.0;
    yofgam[0] = nep[0] / (pow(e_engy[0], 6.0) * SQR(coskiterm) * contd1);

  } else {
    contd1 = etermfactor;
    yofgam[0] = nep[0] / (SQR(e_engy[0]) * SQR(e_engy[0]) * coskiterm * contd1);

  }

  if ((yofgam[0] < 0.0) || (yofgam[0] < 1.0e-100)) {
    yofgam[0] = 0.0;
  }

  for (m = 1; m < EGRID; m++) {

    if (nep[m] == 0.0) {

      yofgam[m] = 0.0;
      continue;
    }

    beta_e = beta(e_engy[m]);
    coskiterm = 1.0 - (beta_e * coski);
    gamterm = SQR(e_engy[m]) * coskiterm;

    expterm = etermfactor / gamterm;

    if (expterm >= 150.0) {
      yofgam[m] = 0.0;
      continue;

    } else if ((expterm >= 24.0) && (expterm < 150.0)) {
      contd1 = exp(expterm);
      yofgam[m] = nep[m] / (pow(e_engy[m], 6.0) * SQR(coskiterm) * contd1);

    } else if ((expterm < 24.0) && (expterm > 1.0e-4)) {
      contd1 = exp(expterm) - 1.0;
      yofgam[m] = nep[m] / (pow(e_engy[m], 6.0) * SQR(coskiterm) * contd1);

    } else {
      contd1 = etermfactor;
      yofgam[m] = nep[m] / (SQR(e_engy[m]) * SQR(e_engy[m]) * coskiterm * contd1);

    }

    if ((yofgam[m] < 0.0) || (yofgam[m] < 1.0e-100)) {
      yofgam[m] = 0.0;
      continue;
    }

    sum += (0.5 * (yofgam[m] + yofgam[m - 1]) * (e_engy[m] - e_engy[m - 1]));

  }

  return (sum);

}


double etablrth_integ(double epstemp, double BLF, double gausmin, double gausmax, double *w, double *gausetap, double *gausetam, double *lineta, double *e_engy, double *nep, double *incont_gausp, double *incont_gausm, double *incont_linblr, int GAUSASZ, double *coskip, double *coskim, double *coskilin)
{
  double yoflineta[ETASTEP], Beta_bulk, term1, term2, sum = 0.0;
  double bphtermp, etermfacp, bphtermm, etermfacm, bphterm, etermfac, xr;
  double sumgaus_cont = 0.0, sumlin_cont = 0.0;

  int j, k = 0, m;

  // For the gaus part of etaph
  //
  Beta_bulk = beta(BLF);

  xr = 0.5 * (gausmax - gausmin);

  for (j = 1; j <= GAUSASZ; j++) {

    bphtermp = 1.0 + (Beta_bulk * gausetap[k]);
    etermfacp = epstemp * BLF * bphtermp;

    bphtermm = 1.0 + (Beta_bulk * gausetam[k]);
    etermfacm = epstemp * BLF * bphtermm;

    term1 = incont_gausp[k] * Int_gamblrth(coskip[k], etermfacp, e_engy, nep);

    term2 = incont_gausm[k] * Int_gamblrth(coskim[k], etermfacm, e_engy, nep);

    sumgaus_cont += (w[j] * (term1 + term2));

    k++;

  }

  sumgaus_cont *= xr;

  // For the linear part of etaph
  //
  bphterm = 1.0 + (Beta_bulk * lineta[0]);
  etermfac = epstemp * BLF * bphterm;

  yoflineta[0] = incont_linblr[0] * Int_gamblrth(coskilin[0], etermfac, e_engy, nep);

  for (j = 1; j < ETASTEP; j++) {

    bphterm = 1.0 + (Beta_bulk * lineta[j]);
    etermfac = epstemp * BLF * bphterm;

    yoflineta[j] = incont_linblr[j] * Int_gamblrth(coskilin[j], etermfac, e_engy, nep);

    sumlin_cont += (0.5 * (yoflineta[j] + yoflineta[j - 1]) * (lineta[j] - lineta[j - 1]));

  }

  sum = sumgaus_cont + sumlin_cont;

  return(sum);

}


double Sum_gamblrthline(double eps, double coskiterm, double *Nnuline, double *peline, double *e_engy, double *nep) {

  double sum = 0.0, term1, term2, gammaline[LINEGRID], lnterm1, exptermline;
  double gamterm, negammaline[LINEGRID], rootval, lnterm2, lnterm3, ln1, ln2, ln3;
  int j = 0, m = 0;

  term2 = eps / coskiterm;

  for (m = 0; m < LINEGRID; m++) {

    rootval = term2 / peline[m];
    gammaline[m] = sqrt(rootval);
    negammaline[m] = 0.0;

  }

  for (m = 0; m < LINEGRID; m++) {

    if (gammaline[m] < 10.0) {
      continue;

    }

    for (j = 1; j < EGRID; j++) {
      
      if (nep[j] == 0.0) {
	continue;
	
      }

      if ((gammaline[m] > e_engy[j - 1]) && (gammaline[m] < e_engy[j])) {

	lnterm1 = gammaline[m] / e_engy[j - 1];
	ln1 = log(lnterm1);
	
	lnterm2 = nep[j] / nep[j - 1];
	ln2 = log(lnterm2);

	lnterm3 = e_engy[j] / e_engy[j - 1];
	ln3 = log(lnterm3);

	exptermline = (ln1 * ln2) / ln3;

	if (fabs(exptermline) > 1.0e-4) {
	  negammaline[m] = exp(exptermline) * nep[j - 1]; 

	} else {
	  negammaline[m] = (1.0 + exptermline) * nep[j - 1];
	  
	}

      }

    }

    term1 = Nnuline[m] / pow(peline[m], 1.5);
    
    sum += (term1 * negammaline[m]);

    if (sum < 0.0) {

      printf("eps = %e, coskiterm = %e\n", eps, coskiterm);
      printf("eps = %e, Nnuline[%d] = %e, peline[%d] = %e, gammaline[%d] = %e, negammaline[%d] = %e\n",
	     eps, m, Nnuline[m], m, peline[m], m, gammaline[m], m, negammaline[m]);
      printf("\nProblem in Sum_gamblrthline\n");
      exit(0);
    }

  }

  return(sum);

}



double etablrth_integline(double eps, double Gamma_bulk, double gausmin, double gausmax, double *w, double *gausetap, double *gausetam, double *lineta, double *Nnuline, double *e_engy, double *nep, double pe_linep[][LINEGRID], double pe_linem[][LINEGRID], double pe_linelin[][LINEGRID], double *inline_gausp, double *inline_gausm, double *inline_linblr, int GAUSASZ, double *coskip, double *coskim, double *coskilin) {

  double sumgaus_line = 0.0, yoflineta[ETASTEP], Beta_bulk, term1, term2, sum = 0.0;
  double sumlin_line = 0.0, xr, bphtermp, etermfacp, denomp, bphtermm, etermfacm;
  double etermfac, denom, pdline_linblr[ETASTEP], denomm, bphterm;
  double pdline_gausblrp[GAUSASZ], pdline_gausblrm[GAUSASZ];
  double coskipterm, coskimterm, coskilinterm;

  int j, k = 0, m;

  // For the gaus part of etaph
  //
  Beta_bulk = beta(Gamma_bulk);

  xr = 0.5 * (gausmax - gausmin);

  for (j = 1; j <= GAUSASZ; j++) {

    bphtermp = 1.0 + (Beta_bulk * gausetap[k]);
    denomp = SQR(bphtermp) * SQR(bphtermp);
    coskipterm = 1.0 - coskip[k];

    bphtermm = 1.0 + (Beta_bulk * gausetam[k]);
    denomm = SQR(bphtermm) * SQR(bphtermm);
    coskimterm = 1.0 - coskim[k];

    pdline_gausblrp[k] = (inline_gausp[k] * sqrt(coskipterm)) / denomp;
    term1 = pdline_gausblrp[k] * Sum_gamblrthline(eps, coskipterm, Nnuline, pe_linep[k], e_engy, nep);

    pdline_gausblrm[k] = (inline_gausm[k] * sqrt(coskimterm)) / denomm;
    term2 = pdline_gausblrm[k] * Sum_gamblrthline(eps, coskimterm, Nnuline, pe_linem[k], e_engy, nep);

    sumgaus_line += (w[j] * (term1 + term2));

    k++;

  }

  sumgaus_line *= xr;

  // For the linear part of eta_ph
  //
  bphterm = 1.0 + (Beta_bulk * lineta[0]);
  denom = SQR(bphterm) * SQR(bphterm);
  coskilinterm = 1.0 - coskilin[0];

  pdline_linblr[0] = (inline_linblr[0] * sqrt(coskilinterm)) / denom;
  yoflineta[0] = pdline_linblr[0] * Sum_gamblrthline(eps, coskilinterm, Nnuline, pe_linelin[0], e_engy, nep);

  for (j = 1; j < ETASTEP; j++) {

    bphterm = 1.0 + (Beta_bulk * lineta[j]);
    denom = SQR(bphterm) * SQR(bphterm);
    coskilinterm = 1.0 - coskilin[j];

    pdline_linblr[j] = (inline_linblr[j] * sqrt(coskilinterm)) / denom;
    yoflineta[j] = pdline_linblr[j] * Sum_gamblrthline(eps, coskilinterm, Nnuline, pe_linelin[j], e_engy, nep);

    sumlin_line += (0.5 * (yoflineta[j] + yoflineta[j - 1]) * (lineta[j] - lineta[j - 1]));

  }

  sum = sumgaus_line + sumlin_line;

  return(sum);

}


double ndotblrth(double eps, double Gamma_bulk, double contlabinteg, double gausetmin, double gausetmax, double *w, double *gausetap, double *gausetam, double *lineta, double *Nnuline, double *e_engy, double *nep, double pe_linep[][LINEGRID], double pe_linem[][LINEGRID], double pe_linelin[][LINEGRID], double *intencont_gausblrp, double *intencont_gausblrm, double *intencont_linblr, double *intenline_gausblrp, double *intenline_gausblrm, double *intenline_linblr, int GAUSASZ, double coskip[][GAUSASZ], double coskim[][GAUSASZ], double coskilin[][ETASTEP], double phistep) {

  double finalterm_cont, finalterm_line, finalterm, sum_cont = 0.0, epstemp;
  double yofphi_cont[PHIGRID], yofphi_line[PHIGRID], sum_line = 0.0;
  int i, j;

  epstemp = eps / AVTEMP;

  yofphi_cont[0] = etablrth_integ(epstemp, Gamma_bulk, gausetmin, gausetmax, w, gausetap, gausetam, lineta, e_engy, nep, intencont_gausblrp, intencont_gausblrm, intencont_linblr, GAUSASZ, coskip[0], coskim[0], coskilin[0]);

  yofphi_line[0] = etablrth_integline(eps, Gamma_bulk, gausetmin, gausetmax, w, gausetap, gausetam, lineta, Nnuline, e_engy, nep, pe_linep, pe_linem, pe_linelin, intenline_gausblrp, intenline_gausblrm, intenline_linblr, GAUSASZ, coskip[0], coskim[0], coskilin[0]);

  for (j = 1; j < PHIGRID; j++) {

    yofphi_cont[j] = etablrth_integ(epstemp, Gamma_bulk, gausetmin, gausetmax, w, gausetap, gausetam, lineta, e_engy, nep, intencont_gausblrp, intencont_gausblrm, intencont_linblr, GAUSASZ, coskip[j], coskim[j], coskilin[j]);
    sum_cont += (0.5 * (yofphi_cont[j] + yofphi_cont[j - 1]) * phistep);

    yofphi_line[j] = etablrth_integline(eps, Gamma_bulk, gausetmin, gausetmax, w, gausetap, gausetam, lineta, Nnuline, e_engy, nep, pe_linep, pe_linem, pe_linelin, intenline_gausblrp, intenline_gausblrm, intenline_linblr, GAUSASZ, coskip[j], coskim[j], coskilin[j]);
    sum_line += (0.5 * (yofphi_line[j] + yofphi_line[j - 1]) * phistep);

  }

  finalterm_cont = CONTCOEFF * SQR(eps) * sum_cont / contlabinteg;

  finalterm_line = LINECOEFF * sum_line / (2.0 * sqrt(eps) * SQR(Gamma_bulk) * SQR(Gamma_bulk));

  finalterm = NDOTBLRTHCOEFF * (finalterm_cont + finalterm_line);

  return(2.0 * finalterm);

}
