#define SSCOEFF 7.48125e-15      // Value of 3c*sigma_t/8


double summ(double z)
{
  double sum = 0.0, nn, sd;
  int n = 1;

  do {

    nn = (double) (n);
    sd = pow(z, -nn) / SQR(nn);
    sum += sd;
    n++;

  } while (sd > 1.0e-10);

  return(sum);

}


double gammasscdot(double *ep, double *nph, double eengy, int n)
{
  double engyprod, gengyfun, sum = 0.0, gammadot, t1, yofep[ASZ];
  int j;

  engyprod = eengy * ep[0];
  t1 = 1.0 + (4.0 * engyprod);

  if (engyprod < 0.25) {
    gengyfun = SQR(engyprod) * (3.555556 - (22.4 * engyprod) +
				(125.44 * SQR(engyprod)));

  } else {
    gengyfun = (8.0 * engyprod * (1.0 + (5.0 * engyprod)) / (3.0 * SQR(t1)))
      - ((4.0 * engyprod / t1) * (0.666667 + (0.5 / engyprod) +
      (0.125 / SQR(engyprod)))) + (log(t1) * (1.0 + (3.0 / engyprod) +
      (0.75 / SQR(engyprod)) + (0.5 * log(t1) / engyprod) - (log(4.0 *
      engyprod) / engyprod))) - (2.5 / engyprod) + (summ(t1) / engyprod) -
      (1.644934 / engyprod) - 2.0;

  }

  yofep[0] = nph[0] * gengyfun / ep[0];

  for (j = 1; j < n; j++) {

    if (nph[j] < 1.0e-15) {
      yofep[j] = 0.0;
      continue;
    }

    engyprod = eengy * ep[j];
    t1 = 1.0 + (4.0 * engyprod);

    if (engyprod < 0.25) {
      gengyfun = SQR(engyprod) * (3.555556 - (22.4 * engyprod) +
				  (125.44 * SQR(engyprod)));

    } else {
      gengyfun = (8.0 * engyprod * (1.0 + (5.0 * engyprod)) / (3.0 * SQR(t1)))
	- ((4.0 * engyprod / t1) * (0.666667 + (0.5 / engyprod) +
        (0.125 / SQR(engyprod)))) + (log(t1) * (1.0 + (3.0 / engyprod) +
        (0.75 / SQR(engyprod)) + (0.5 * log(t1) / engyprod) - (log(4.0 *
	engyprod) / engyprod))) - (2.5 / engyprod) + (summ(t1) / engyprod) -
        (1.644934 / engyprod) - 2.0;

    }

    yofep[j] = nph[j] * gengyfun / ep[j];
    sum += (0.5 * (yofep[j] + yofep[j - 1]) * (ep[j] - ep[j - 1]));

  }

  gammadot = SSCOEFF * sum;

  return(-gammadot);

}



double epint(double scatterphotengy, double *photengy, double eengy, double *nph, int n)
{
  double sum = 0.0, engyratio, t1, eLower, eUpper, g_es_e_gam, funofep[ASZ], eengy_param;
  int j;
  double t1_param, eUpper_param, engyratio_param, scatter_param, gam_param, scatter_engy_param;

  t1 = 4.0 * SQR(eengy) * photengy[0];
  eLower = photengy[0] / (4.0 * SQR(eengy));
  eUpper = t1 / (1.0 + (4.0 * photengy[0] * eengy));
  engyratio = scatterphotengy / (t1 * (1.0 - (scatterphotengy / eengy)));

  if ((scatterphotengy >= eLower) && (scatterphotengy <= photengy[0])) {
    g_es_e_gam = (0.25 * scatterphotengy / SQR((eengy * photengy[0]))) -
      (6.25e-02 / (SQR((eengy * eengy)) * photengy[0]));
    funofep[0] = nph[0] * g_es_e_gam;

  } else if ((scatterphotengy > photengy[0]) && (scatterphotengy <= eUpper)) {
    g_es_e_gam = (1.0 / (photengy[0] * SQR(eengy))) *
      ((0.5 * engyratio * log(engyratio)) +
       (0.25 * (1.0 + (2.0 * engyratio)) * (1.0 - engyratio)) +
       (0.125 * SQR(scatterphotengy) * (1.0 - engyratio) / (eengy * (eengy - scatterphotengy))));
    funofep[0] = nph[0] * g_es_e_gam;

    } else funofep[0] = 0.0;

  t1_param = 4.0 * SQR(eengy);
  eUpper_param = 4.0 * eengy;
  engyratio_param = scatterphotengy / (4.0 * SQR(eengy) * (1.0 - (scatterphotengy / eengy)));
  scatter_param = 0.25 * scatterphotengy / SQR(eengy);
  gam_param = 6.25e-02 / SQR((eengy * eengy));
  eengy_param = 1.0 / SQR(eengy);
  scatter_engy_param = 0.125 * SQR(scatterphotengy) / (eengy * (eengy - scatterphotengy));

  for (j = 1; j < n; j++) {

    if (nph[j] < 1.0e-15) {
      funofep[j] = 0.0;
      continue;
    }

    t1 = t1_param * photengy[j];
    eLower = photengy[j] / t1_param;
    eUpper = t1 / (1.0 + (eUpper_param * photengy[j]));
    engyratio = engyratio_param / photengy[j];

    if ((scatterphotengy >= eLower) && (scatterphotengy <= photengy[j])) {
      g_es_e_gam = (scatter_param / SQR(photengy[j])) - (gam_param / photengy[j]);
      funofep[j] = nph[j] * g_es_e_gam;

    } else if ((scatterphotengy > photengy[j]) && (scatterphotengy <= eUpper)) {
      g_es_e_gam = (eengy_param / photengy[j]) * ((0.5 * engyratio * log(engyratio)) +
	           (0.25 * (1.0 + (2.0 * engyratio)) * (1.0 - engyratio)) +
                   (scatter_engy_param * (1.0 - engyratio)));
      funofep[j] = nph[j] * g_es_e_gam;

    } else funofep[j] = 0.0;

    sum += (0.5 * (funofep[j] + funofep[j-1]) * (photengy[j] - photengy[j-1]));

  }

  return(sum);

}



double ndotssc(double scatterphotengy, double *photengy, double *e_engy, double *nep, double *nph, int n)
{
  double sum = 0.0, inte, funofgam[EGRID];
  int j;

  if (e_engy[0] < scatterphotengy) funofgam[0] = 0.0;
  else funofgam[0] = nep[0] * epint(scatterphotengy, photengy, e_engy[0], nph, n);

  for (j = 1; j < EGRID; j++) {

    if (nep[j] == 0.0) {

      funofgam[j] = 0.0;
      continue;

    }

    if (e_engy[j] < scatterphotengy) funofgam[j] = 0.0;
    else funofgam[j] = nep[j] * epint(scatterphotengy, photengy, e_engy[j], nph, n);

    sum += (0.5 * (funofgam[j] + funofgam[j - 1]) * (e_engy[j] - e_engy[j - 1]));

  }

  return(SSCOEFF * 8.0 * sum);

}
