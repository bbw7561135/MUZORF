#define GGCOEFF 1.49625e-14  // Value of 3c*sigma_t/4

double sigmagg(double scatterphotengy, double *photengy, double *nph, double muval, int n)
{
  double sum = 0.0, rootx, x, sigxgg, t1, eLower, nphmid, epmid, yofsig[n];
  int j, k = 0;

  eLower = 2.0 / (scatterphotengy * (1.0 - muval));

  while ((k < n) && (photengy[k] < eLower)) {

    k++;

  }

  if ((k < n) && (nph[k] > 1.0e-15)) {

    epmid = 0.5 * (photengy[k] + eLower);
    x = eLower / epmid;

    if (x < 1) {

      rootx = sqrt((1.0 - x));
      t1 = (1.0 + rootx) / (1.0 - rootx);
      sigxgg = (1.0 - SQR(rootx)) * (((3.0 - SQR((rootx * rootx))) * log(t1))
				     - ((2.0 * rootx) * (2.0 - SQR(rootx))));

    } else sigxgg = 0.0;

    sum += ((photengy[k] - eLower) * sigxgg * nph[k]);
    k++;

  }

  if (k < n) {

    if (nph[k] <= 1.0e-15) {

      yofsig[k] = 0.0;

    } else {

      x = eLower / photengy[k];

      if (x < 1) {

	rootx = sqrt((1.0 - x));
	t1 = (1.0 + rootx) / (1.0 - rootx);
	sigxgg = (1.0 - SQR(rootx)) * (((3.0 - SQR((rootx * rootx))) * log(t1))
				       - ((2.0 * rootx) * (2.0 - SQR(rootx))));

      } else sigxgg = 0.0;

      yofsig[k] = sigxgg * nph[k];

    }

    for (j = (k + 1); j < n; j++) {

      if (nph[j] <= 1.0e-15) {

	yofsig[j] = 0.0;
	continue;

      }

      x = eLower / photengy[j];

      if (x < 1) {

	rootx = sqrt((1.0 - x));
	t1 = (1.0 + rootx) / (1.0 - rootx);
	sigxgg = (1.0 - SQR(rootx)) * (((3.0 - SQR((rootx * rootx))) * log(t1))
				       - ((2.0 * rootx) * (2.0 - SQR(rootx))));

      } else sigxgg = 0.0;

      yofsig[j] = sigxgg * nph[j];
      sum += (0.5 * (yofsig[j] + yofsig[j - 1]) * (photengy[j] - photengy[j - 1]));

    }

  }

  return(sum);

}



double taugg(double scatterphotengy, double *photengy, double *nph, int n)
{
  double sum = 0.0, mumin = -1.0, steps, mumax = 1.0, mugg[MUGRID], funofmugg[MUGRID];
  int j;

  steps = (mumax - mumin) / ((double) (MUGRID));

  mugg[0] = mumin;
  funofmugg[0] = (1.0 - mugg[0]) * sigmagg(scatterphotengy, photengy, nph, mugg[0], n);
  mumin += steps;

  for (j = 1; j < MUGRID; j++) {

    mugg[j] = mumin;
    funofmugg[j] = (1.0 - mugg[j]) * sigmagg(scatterphotengy, photengy, nph, mugg[j], n);
    sum += (0.5 * (funofmugg[j] + funofmugg[j - 1]) * steps);
    mumin += steps;

  }

  return(sum);

}


// Pair Injection.
//
double I(double photengy1, double photengy2, double ecmUL, double cpm)
{
  double Icmpm, ecpm;

  if (cpm > 0.0) {

    ecpm = (photengy1 * photengy2) + (cpm * SQR(ecmUL));
    Icmpm = (1.0 / sqrt(fabs(cpm))) * log((ecmUL * sqrt(fabs(cpm))) + sqrt(ecpm));

  } else {

    Icmpm = (1.0 / sqrt(fabs(cpm))) * asin(ecmUL * sqrt((fabs(cpm) / (photengy1 * photengy2))));

  }

  return(Icmpm);

}



double H(double photengy1, double photengy2, double ecmUL, double cpm, double dpm)
{
  double Hcmpm, ecpm;

  if (fabs(cpm) < 1.0e-5) {

    Hcmpm = (((0.083333 * pow(ecmUL, 3.0)) - (0.125 * ecmUL * dpm)) *
	     (1.0 / pow((photengy1 * photengy2), 1.5))) +
      (((0.166667 * pow(ecmUL, 3.0)) + (0.5 * ecmUL) + (0.25 / ecmUL)) *
       (1.0 / sqrt((photengy1 * photengy2))));

  } else {

    ecpm = (photengy1 * photengy2) + (cpm * SQR(ecmUL));
    Hcmpm = (0.25 * (2.0 - (((photengy1 * photengy2) - 1.0) / cpm)) *
	     I(photengy1, photengy2, ecmUL, cpm)) -
      (0.125 * ecmUL / sqrt(ecpm)) *
      ((dpm / (photengy1 * photengy2)) + (2.0 / cpm)) +
      (0.25 * sqrt(ecpm) * ((ecmUL / cpm) +
			    (1.0 / (ecmUL * photengy1 * photengy2))));

  }

  return(Hcmpm);

}



double ecm(double photengy1, double photengy2, double eengy)
{
  double sumengy, ecmLower, ecmstar, ecmUpper, ecmdag, dplus, dminus, cplus;
  double ecminus, engyterm, ecmplus, t1, tcmL, t2, tcmU, cminus, diff;

  sumengy = photengy1 + photengy2;
  if (sumengy < eengy) return(0.0);

  engyterm = SQR(((eengy * (sumengy - eengy)) + 1.0)) - SQR(sumengy);
  if (engyterm <= 0.0) return(0.0);

  ecminus = 0.5 * ((eengy * (sumengy - eengy)) + 1.0 - sqrt(engyterm));
  ecmdag = sqrt(ecminus);
  ecmLower = max(1.0, ecmdag);
  t1 = SQR(sumengy) - (4.0 * SQR(ecmLower));
  tcmL = 0.25 * sqrt(t1);

  ecmplus = 0.5 * ((eengy * (sumengy - eengy)) + 1.0 + sqrt(engyterm));
  ecmstar = sqrt(ecmplus);
  ecmUpper = min(sqrt((photengy1 * photengy2)), ecmstar);
  t2 = SQR(sumengy) - (4.0 * SQR(ecmUpper));
  tcmU = 0.25 * sqrt(t2);

  cplus = SQR((photengy1 - eengy)) - 1.0;
  cminus = SQR((photengy2 - eengy)) - 1.0;
  dplus = SQR(photengy1) + (photengy1 * photengy2) + (eengy * (photengy2 - photengy1));
  dminus = SQR(photengy2) + (photengy1 * photengy2) - (eengy * (photengy2 - photengy1));

  if (ecmLower >= ecmUpper) return(0.0);

  else {

    diff = tcmU + H(photengy1, photengy2, ecmUpper, cplus, dplus) +
      H(photengy1, photengy2, ecmUpper, cminus, dminus) - tcmL -
      H(photengy1, photengy2, ecmLower, cplus, dplus) - H(photengy1, photengy2, ecmLower, cminus, dminus);

    return(diff);

  }

}



double nph2_int(double photengy1, double eengy, double *photengy, double *nph, int n)
{
  double epmid2, phot2Lower, sum = 0.0, yofpe[n];
  int j, k = 0;

  phot2Lower = max((1.0 / photengy1), (eengy + 1.0 - photengy1));

  while ((k < n) && (photengy[k] < phot2Lower)) k++;

  if (k < n) {

    epmid2 = 0.5 * (photengy[k] + phot2Lower);
    sum += ((photengy[k] - phot2Lower) * nph[k] * ecm(photengy1, epmid2, eengy) / SQR(epmid2));
    k++;

  }

  if (k < n) {

    yofpe[k] = nph[k] * ecm(photengy1, photengy[k], eengy) / SQR(photengy[k]);

    for (j = (k + 1); j < n; j++) {

      yofpe[j] = nph[j] * ecm(photengy1, photengy[j], eengy) / SQR(photengy[j]);
      sum += (0.5 * (yofpe[j] + yofpe[j - 1]) * (photengy[j] - photengy[j - 1]));

    }

  }

  return(sum);

}



double nedotgg(double *photengy, double eengy, double *nph, int n)
{
  double sum = 0.0, funofe1[n];
  int j;

  if (eengy <= 0.0) {

    printf("\n\nIn function nedotgg : Value of eengy is incorrect -> %e\n\n", eengy);
    exit(Exit_success);
  }

  funofe1[0] = nph[0] * nph2_int(photengy[0], eengy, photengy, nph, n) /
    SQR(photengy[0]);

  for (j = 1; j < n; j++) {

    funofe1[j] = nph[j] * nph2_int(photengy[j], eengy, photengy, nph, n) /
      SQR(photengy[j]);
    sum += (0.5 * (funofe1[j] + funofe1[j - 1]) * (photengy[j] - photengy[j - 1]));

  }

  return(GGCOEFF * sum);

}
