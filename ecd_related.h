#define ECDCOEFF 1.707299e15  // Value of (pi^5*r_e^2/30c^2)(m_e*c^2/h)^3
#define ERCOEFF 2.70130 // Value of Gamma(4)*RZeta(4)/Gamma(3)*RZeta(3)
#define ECDTHCOEFF 3.80716e16 // Value of (4*pi^5*sigma_T/45c^2)(m_e*c^2/h)^3
#define NDOTTHCOEFF 6.99802e14 // Value of (sigma_T/(2*c^2))*(m_e*c^2/h)^3
//#define NDOTCOEFF 4.17663e13 // Value of (3sigma_T/(32*pi*c^2))*(m_e*c^2/h)^3
#define NECDCOEFF 2.62426e14 // Value of (3sigma_T/(16*c^2))*(m_e*c^2/h)^3
#define ECDISOCOEFF 1.64887e15 // Value of (3*pi*sigma_T/8*c^2) * (m_e*c^2/h)^3
#define NECDISOCOEFF 2.204e29 // Value of (pi/c^3) * (m_e*c^2/h)^3


// The emission region height, zht, and disk radius array, raddist[GS], need to
// be in units of pc for this entire subroutine.
//

// Start of ECD gammadot calculation using Compton cross-section under
// Markus' eta_e = 0 approximation. Based on eqn.(13) of Boettcher et al., 1997
//
double Iegeta(double avegam, double engy, double etavalue)
{
  double Eterm, Aterm, Dterm, betaterm, Ivalue, beta_engy, term4;
  double s, p, sq, psqe, Nterm, lnterm, m, mad, ead, term1, term2, term3;

  beta_engy = beta(engy);
  betaterm = beta_engy * sqrt((1.0 - SQR(etavalue)));

  Eterm = engy * avegam;
  if (Eterm < 0.0) return(0.0);

  //if (Eterm < 0.25) {
  if (Eterm < 0.1) {
    Ivalue = (2.66667 * (1.0 + (0.5 * SQR(betaterm)))) -
      (11.2 * Eterm * (1.0 + (1.5 * SQR(betaterm)))) +
      (18.4 * SQR(Eterm) * (1.0 + (3.0 * SQR(betaterm)) +
			    (SQR(betaterm) * SQR(betaterm))));

  } else {
    s = 1.0 + (2.0 * Eterm);
    p = 2.0 * Eterm * betaterm;
    sq = sqrt((1.0 - SQR(betaterm)));
    psqe = Eterm * sq;
    Nterm = sqrt((SQR(s) - SQR(p)));
    Aterm = (s / p) - sqrt((SQR((s / p)) - 1.0));
    Dterm = (1.0 - sq) / betaterm;
    lnterm = log(Eterm * betaterm / Aterm);
    m = 1.0 - SQR(Dterm);
    mad = 1.0 - (Aterm * Dterm);
    ead = Eterm * betaterm * m;

    term1 = (6.0 / psqe) - 0.83333 + log(0.5 * (s + Nterm));
    term2 = (2.0 * lnterm / psqe) + (8.0 * Dterm * log(mad) / ead) +
      (3.0 * lnterm / (SQR(Eterm) * pow(sq, 3.0)));
    term3 = (24.0 * SQR(Dterm) / SQR(ead)) *
      ((((1.0 + SQR(Dterm)) / m) * log(mad)) - (Aterm * Dterm / mad));
    term4 = (((2.0 * SQR(s)) + SQR(p)) / (6.0 * pow(Nterm, 5.0))) -
      (s / (2.0 * pow(Nterm, 3.0))) - (1.0 / Nterm);

    Ivalue = (term1 - term2 - term3 + term4) / SQR(Eterm);

  }

  if (Ivalue < 0.0) Ivalue = 0.0;

  return(Ivalue * 2.0 * pi);

}


double gammaecdiskdot(double eengy, double Gamma_bulk, double zht, double *thetar, double *raddist)
{
  double sum = 0.0, etaph, x, rzterm, xterm;
  double yofr[RGRID], gammadot, Beta_bulk, eave;
  int j;

  Beta_bulk = beta(Gamma_bulk);

  rzterm = 0.5 * SQR((raddist[0] / zht));
  if (rzterm < 1.0e-6) {    //Some small number to enable binomial expansion.
    x = zht * (1.0 + rzterm);
    etaph = (1.0 - (Beta_bulk * rzterm) - Beta_bulk) /
      (1.0 + rzterm - Beta_bulk);

  } else {
    xterm = SQR(raddist[0]) + SQR(zht);
    x = sqrt(xterm);
    etaph = (zht - (Beta_bulk * x)) / (x - (Beta_bulk * zht));

  }

  if (fabs(etaph) < 0.99999999) {
    eave = ERCOEFF * thetar[0] / (Gamma_bulk * (1.0 + (Beta_bulk * etaph)));

    yofr[0] = raddist[0] * SQR(thetar[0]) * SQR(thetar[0]) *
      SQR((x - (Beta_bulk * zht))) * Iegeta(eave, eengy, etaph) /
      (SQR(x) * SQR(x));
  } else yofr[0] = 0.0;

  for (j = 1; j < RGRID; j++) {

    rzterm = 0.5 * SQR((raddist[j] / zht));
    if (rzterm < 1.0e-6) {
      x = zht * (1.0 + rzterm);
      etaph = (1.0 - (Beta_bulk * rzterm) - Beta_bulk) /
	(1.0 + rzterm - Beta_bulk);

    } else {
      xterm = SQR(raddist[j]) + SQR(zht);
      x = sqrt(xterm);
      etaph = (zht - (Beta_bulk * x)) / (x - (Beta_bulk * zht));

    }

    if (fabs(etaph) < .99999999) {
      eave = ERCOEFF * thetar[j] / (Gamma_bulk * (1.0 + (Beta_bulk * etaph)));

      yofr[j] = raddist[j] * SQR(thetar[j]) * SQR(thetar[j]) *
	SQR((x - (Beta_bulk * zht))) * Iegeta(eave, eengy, etaph) /
	(SQR(x) * SQR(x));
    } else yofr[j] = 0.0;

    sum += (0.5 * (yofr[j] + yofr[j - 1]) * (raddist[j] - raddist[j - 1]));

  }

  gammadot = ECDCOEFF * SQR((eengy * Gamma_bulk)) * sum;

  return(-gammadot);

}


// Start of ECD gammadot calculation using Thomson cross-section under
// Markus' eta_e = 0 approximation. Based on eqn.(12) of Boettcher et al., 1997
//
double gammaecdthdot(double eengy, double Gamma_bulk, double zht, double *thetar, double *raddist)
{
  double gammathdot, x, xterm, sum = 0.0, yofr[RGRID], rzterm, Beta_bulk;
  int j;

  Beta_bulk = beta(Gamma_bulk);

  rzterm = 0.5 * SQR((raddist[0] / zht));
  if (rzterm < 1.0e-6) x = zht * (1.0 + rzterm);
  else {
    xterm = SQR(raddist[0]) + SQR(zht);
    x = sqrt(xterm);
  }

  yofr[0] = raddist[0] * SQR(thetar[0]) * SQR(thetar[0]) *
    SQR((x - (Beta_bulk * zht))) / (SQR(x) * SQR(x));

  for (j = 1; j < RGRID; j++) {

    rzterm = 0.5 * SQR((raddist[j] / zht));
    if (rzterm < 1.0e-6) x = zht * (1.0 + rzterm);
    else {
      xterm = SQR(raddist[j]) + SQR(zht);
      x = sqrt(xterm);
    }

    yofr[j] = raddist[j] * SQR(thetar[j]) * SQR(thetar[j]) *
      SQR((x - (Beta_bulk * zht))) / (SQR(x) * SQR(x));

    sum += (0.5 * (yofr[j] + yofr[j - 1]) * (raddist[j] - raddist[j - 1]));

  }

  gammathdot = ECDTHCOEFF * SQR((eengy * Gamma_bulk)) * sum;

  return(-gammathdot);

}


// Start of ECD ndot calculation in the Thomson-regime under head-on approx.
// Based on G.I.(2.1.52) of Markus' thesis.
//
double Int_phi(double etaph, double etas, double beta_gam, double expterm)
{
  double sum = 0.0, phimin = 0.0, phimax = pi, steps, coski, bmterm, expindex;
  double yofphi[PHIGRID];
  int j;

  steps = (phimax - phimin) / ((double) (PHIGRID - 1));

  coski = (etaph * etas) + (sqrt(((1.0 - SQR(etaph)) * (1.0 - SQR(etas)))) *
			 cos(phimin));
  if (coski < -.99999999) coski = -.99999999;
  if (coski > .99999999) coski = .99999999;

  bmterm = 1.0 - (beta_gam * coski);
  expindex = expterm / bmterm;

  if (expindex >= 150.0) {
    yofphi[0] = 0.0;

  } else if ((expindex >= 24.0) && (expindex < 150.0)) {
    yofphi[0] = 1.0 / (SQR(bmterm) * exp(expindex));

  } else if ((expindex < 24.0) && (expindex > 1.0e-4)) {
    yofphi[0] = 1.0 / (SQR(bmterm) * (exp(expindex) - 1.0));

  } else {
    yofphi[0] = 1.0 / (bmterm * expterm);
  }

  if (yofphi[0] < 1.0e-100) yofphi[0] = 0.0;

  phimin += steps;

  for (j = 1; j < PHIGRID; j++) {

    coski = (etaph * etas) + (sqrt(((1.0 - SQR(etaph)) * (1.0 - SQR(etas)))) *
			 cos(phimin));
    if (coski < -.99999999) coski = -.99999999;
    if (coski > .99999999) coski = .99999999;

    bmterm = 1.0 - (beta_gam * coski);
    expindex = expterm / bmterm;

    if (expindex >= 150.0) {
      phimin += steps;
      yofphi[j] = 0.0;
      continue;

    } else if ((expindex >= 24.0) && (expindex < 150.0)) {
      yofphi[j] = 1.0 / (SQR(bmterm) * exp(expindex));

    } else if ((expindex < 24.0) && (expindex > 1.0e-4)) {
      yofphi[j] = 1.0 / (SQR(bmterm) * (exp(expindex) - 1.0));

    } else {
      yofphi[j] = 1.0 / (bmterm * expterm);
    }

    if (yofphi[j] < 1.0e-100) {
      phimin += steps;
      yofphi[j] = 0.0;
      continue;
    }

    sum += (0.5 * (yofphi[j] + yofphi[j - 1]) * steps);

    phimin += steps;

  }

  return(2.0 * sum);

}

double Rinteg(double gam, double eps, double muobsprime, double Gamma_bulk, double zht, double *thetar, double *raddist)
{
  double x, xterm, expterm, denom, sum = 0.0, yofr[RGRID];
  double etaph, xBzterm, beta_gam, Beta_bulk, rzterm;
  int j;

  beta_gam = beta(gam);
  Beta_bulk = beta(Gamma_bulk);

  rzterm = 0.5 * SQR((raddist[0] / zht));
  if (rzterm < 1.0e-6) {
    x = zht * (1.0 + rzterm);
    etaph = (1.0 - (Beta_bulk * rzterm) - Beta_bulk) /
      (1.0 + rzterm - Beta_bulk);

  } else {
    xterm = SQR(raddist[0]) + SQR(zht);
    x = sqrt(xterm);
    etaph = (zht - (Beta_bulk * x)) / (x - (Beta_bulk * zht));

  }

  xBzterm = x - (Beta_bulk * zht);
  denom = Gamma_bulk * SQR(gam) * thetar[0] * xBzterm;
  expterm = (eps * x) / denom;

  if (fabs(etaph) < .99999999) {
    yofr[0] = raddist[0] * Int_phi(etaph, muobsprime, beta_gam, expterm) /
      SQR(xBzterm);
  } else yofr[0] = 0.0;

  for (j = 1; j < RGRID; j++) {

    rzterm = 0.5 * SQR((raddist[j] / zht));
    if (rzterm < 1.0e-6) {
      x = zht * (1.0 + rzterm);
      etaph = (1.0 - (Beta_bulk * rzterm) - Beta_bulk) /
	(1.0 + rzterm - Beta_bulk);
    } else {
      xterm = SQR(raddist[j]) + SQR(zht);
      x = sqrt(xterm);
      etaph = (zht - (Beta_bulk * x)) / (x - (Beta_bulk * zht));
    }

    xBzterm = x - (Beta_bulk * zht);
    denom = Gamma_bulk * SQR(gam) * thetar[j] * xBzterm;
    expterm = (eps * x) / denom;

    if (fabs(etaph) < .99999999) {
      yofr[j] = raddist[j] *Int_phi(etaph, muobsprime, beta_gam, expterm) /
	SQR(xBzterm);
    } else yofr[j] = 0.0;

    sum += (0.5 * (yofr[j] + yofr[j - 1]) * (raddist[j] - raddist[j - 1]));

  }

  return(sum);

}

double ndotth(double eps, double Gamma_bulk, double zht, double muobsprime, double *thetar, double *raddist, double *e_engy, double *nep)
{
  double gamterm, finalterm, sum = 0.0, yofengy[EGRID];
  int j;

  if (!(fabs(muobsprime) < .99999999)) {
    printf ("\n muobsprime = %e !! ", muobsprime);
    return(0.0);
  }

  if (nep[0] == 0.0) yofengy[0] = 0.0;

  else {
    gamterm = nep[0] / (SQR(e_engy[0]) * SQR(e_engy[0]) * SQR(e_engy[0]));
    yofengy[0] = gamterm * Rinteg(e_engy[0], eps, muobsprime, Gamma_bulk, zht, thetar, raddist);
  }

  for (j = 1; j < EGRID; j++) {

    if (nep[j] == 0.0) {
      yofengy[j] = 0.0;
      continue;
    }

    gamterm = nep[j] / (SQR(e_engy[j]) * SQR(e_engy[j]) * SQR(e_engy[j]));
    yofengy[j] = gamterm * Rinteg(e_engy[j], eps, muobsprime, Gamma_bulk, zht, thetar, raddist);

    sum += (0.5 * (yofengy[j] + yofengy[j - 1]) * (e_engy[j] - e_engy[j - 1]));

  }

  finalterm = SQR((eps / Gamma_bulk)) * sum;

  return(NDOTTHCOEFF * finalterm);

}


// Start of ECD ndot calculation using Compton cross-section under
// head-on approx. Based on eqn. derived using Dermer & Menon book, Ch-6.
//
double Int_gamecd(double coski, double scatterphotengy, double photengy, double *e_engy, double *nep)
{
  double sum = 0.0, yofgam[EGRID], cross_sec, gamLL, gammid, term1;
  double gamepsterm, epepsterm, coskiterm;
  int j, k = 0;

  coskiterm = 1.0 - coski;
  term1 = 2.0 / (photengy * scatterphotengy * coskiterm);
  gamLL = 0.5 * scatterphotengy * (1.0 + sqrt(1.0 + term1));
  gamLL = max(gamLL, e_engy[0]);

  while ((k < EGRID) && (e_engy[k] < gamLL)) {
    k++;
  }

  if (((k < EGRID) && (nep[k] > 0.0)) && (gamLL > e_engy[0])) {

    gammid = 0.5 * (e_engy[k] + gamLL);

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

    yofgam[k] = (nep[k] * cross_sec) / SQR(gammid);

    sum += ((e_engy[k] - gamLL) * yofgam[k]);

    if (yofgam[k] < 0.0) {
      sum = 0.0;
    }

    k++;

  }

  if (k < EGRID) {

    if (nep[k] == 0.0) {

      yofgam[k] = 0.0;

    } else {

      if (e_engy[k] >= (0.1 / photengy)) {

	gamepsterm = (e_engy[k] - scatterphotengy) / e_engy[k];
	epepsterm = scatterphotengy / (photengy * e_engy[k] * coskiterm *
				       (e_engy[k] - scatterphotengy));
	
	cross_sec = gamepsterm + (1.0 / gamepsterm) - (2.0 * epepsterm) +
	  SQR(epepsterm);

      } else {

	epepsterm = scatterphotengy / (SQR(e_engy[k]) * photengy * coskiterm);
	cross_sec = 2.0 - (2.0 * epepsterm) + SQR(epepsterm);

      }

      yofgam[k] = (nep[k] * cross_sec) / SQR(e_engy[k]);

      if (yofgam[k] < 0.0) {
	yofgam[k] = 0.0;
      }

    }

    for (j = (k + 1); j < EGRID; j++) {

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


double Int_epecd(double coski, double scatterp_engy, double epsmax, double expterm, double *e_engy, double *nep, double *p_engy)
{
  double sum = 0.0, yofep[DISKGRID], epUL, epmid, expindex;
  int j, k = 0;

  epUL = (2.0 * scatterp_engy) / (1.0 - coski);
  epUL = min(epUL, epsmax);

  expindex = p_engy[0] * expterm;

  if (expindex >= 150.0) {
    yofep[0] = 0.0;

  } else if ((expindex >= 24.0) && (expindex < 150.0)) {
    yofep[0] = (p_engy[0] / exp(expindex)) *
      Int_gamecd(coski, scatterp_engy, p_engy[0], e_engy, nep);

  } else if ((expindex < 24.0) && (expindex > 1.0e-4)) {
    yofep[0] = (p_engy[0] / (exp(expindex) - 1.0)) *
      Int_gamecd(coski, scatterp_engy, p_engy[0], e_engy, nep);

  } else {
    yofep[0] = (1.0 / expterm) *
      Int_gamecd(coski, scatterp_engy, p_engy[0], e_engy, nep);
  }

  if (yofep[0] < 1.0e-100) yofep[0] = 0.0;

  for (j = 1; j < DISKGRID; j++) {

    if (p_engy[j] > epUL) {
      break;
    }

    expindex = p_engy[j] * expterm;

    if (expindex >= 150.0) {
      yofep[j] = 0.0;
      continue;

    } else if ((expindex >= 24.0) && (expindex < 150.0)) {
      yofep[j] = (p_engy[j] / exp(expindex)) *
	Int_gamecd(coski, scatterp_engy, p_engy[j], e_engy, nep);

    } else if ((expindex < 24.0) && (expindex > 1.0e-4)) {
      yofep[j] = (p_engy[j] / (exp(expindex) - 1.0)) *
	Int_gamecd(coski, scatterp_engy, p_engy[j], e_engy, nep);

    } else {
      yofep[j] = (1.0 / expterm) *
	Int_gamecd(coski, scatterp_engy, p_engy[j], e_engy, nep);
    }

    if (yofep[j] < 1.0e-100) {
      yofep[j] = 0.0;
      continue;
    }

    sum += (0.5 * (yofep[j] + yofep[j - 1]) * (p_engy[j] - p_engy[j - 1]));

  }

return (sum);

}


double Recd_integ(double phival, double eps, double etas, double Gamma_bulk, double zht, double *thetar, double *raddist, double *e_engy, double *nep, double *xecd, double *etaphd, double p_engy[][DISKGRID])
{
  double expterm, denom, sum = 0.0, yofr[RGRID];
  double xBzterm, Beta_bulk, coski; 
  int j, i;

  Beta_bulk = beta(Gamma_bulk);

  xBzterm = xecd[0] - (Beta_bulk * zht);
  denom = Gamma_bulk * thetar[0] * xBzterm;
  expterm = xecd[0] / denom;

  coski = (etaphd[0] * etas) + (sqrt(((1.0 - SQR(etaphd[0])) * (1.0 - SQR(etas)))) *
			    cos(phival));

  if (coski < -0.99999999) coski = -0.99999999;
  if (coski > 0.99999999) coski = 0.99999999;

  yofr[0] = raddist[0] * Int_epecd(coski, eps, p_engy[0][DISKGRID - 1], expterm, e_engy, nep, p_engy[0]) /
    SQR(xBzterm);

  for (j = 1; j < RGRID; j++) {

    xBzterm = xecd[j] - (Beta_bulk * zht);
    denom = Gamma_bulk * thetar[j] * xBzterm;
    expterm = xecd[j] / denom;

    coski = (etaphd[j] * etas) + (sqrt(((1.0 - SQR(etaphd[j])) * (1.0 - SQR(etas)))) *
			      cos(phival));
    
    if (coski < -0.99999999) coski = -0.99999999;
    if (coski > 0.99999999) coski = 0.99999999;

    yofr[j] = raddist[j] * Int_epecd(coski, eps, p_engy[j][DISKGRID - 1], expterm, e_engy, nep, p_engy[j]) /
      SQR(xBzterm);

    sum += (0.5 * (yofr[j] + yofr[j - 1]) * (raddist[j] - raddist[j - 1]));

  }

  return(sum);

}


double ndotecd(double eps, double Gamma_bulk, double zht, double muobsprime, double *thetar, double *raddist, double *e_engy, double *nep, double *xecd, double *etaphd, double p_engy[][DISKGRID], double *phiecd, double phistep)
{
  double finalterm, sum = 0.0, yofphi[PHIGRID];
  int j;

  if (!(fabs(muobsprime) < .99999999)) {
    printf ("\n muobsprime = %e !! ", muobsprime);
    return (0.0);
  }

  yofphi[0] = Recd_integ(phiecd[0], eps, muobsprime, Gamma_bulk, zht, thetar, raddist, e_engy, nep, xecd, etaphd, p_engy);

  for (j = 1; j < PHIGRID; j++) {

    yofphi[j] = Recd_integ(phiecd[j], eps, muobsprime, Gamma_bulk, zht, thetar, raddist, e_engy, nep, xecd, etaphd, p_engy);

    sum += (0.5 * (yofphi[j] + yofphi[j - 1]) * phistep);

  }

  finalterm = (NECDCOEFF * sum) / SQR(Gamma_bulk);

  return(2.0 * finalterm);

}



// Start of ECD gammadot calculation using Boettcher et al. 1997 formula
//
double crossseciso(double prodiso)
{
  double gengyiso, t1iso;

  t1iso = 1.0 + (4.0 * prodiso);

  if (prodiso < 0.25) {
    gengyiso = SQR(prodiso) * (3.555556 - (22.4 * prodiso) +
				(125.44 * SQR(prodiso)));

  } else {

    gengyiso = (8.0 * prodiso * (1.0 + (5.0 * prodiso)) / (3.0 * SQR(t1iso)))
      - ((4.0 * prodiso / t1iso) * (0.666667 + (0.5 / prodiso) +
	(0.125 / SQR(prodiso)))) + (log(t1iso) * (1.0 + (3.0 / prodiso) +
	(0.75 / SQR(prodiso)) + (0.5 * log(t1iso) / prodiso) - (log(4.0 *
	 prodiso) / prodiso))) - (2.5 / prodiso) + (summ(t1iso) / prodiso) -
        (1.644934 / prodiso) - 2.0;

  }

  return(gengyiso);

}


double Int_epecdiso(double eengy, double expterm, double *p_engy)
{
  double sum = 0.0, yofep[DISKGRID], expindex, prodecd;
  int j;

  expindex = p_engy[0] * expterm;
  prodecd = eengy * p_engy[0];

  if (expindex >= 150.0) {
    yofep[0] = 0.0;

  } else if ((expindex >= 24.0) && (expindex < 150.0)) {
    yofep[0] = p_engy[0] * crossseciso(prodecd) / exp(expindex);

  } else if ((expindex < 24.0) && (expindex > 1.0e-4)) {
    yofep[0] = p_engy[0] * crossseciso(prodecd) / (exp(expindex) - 1.0);

  } else {
    yofep[0] = crossseciso(prodecd) / expterm;

  }

  if (yofep[0] < 1.0e-100) yofep[0] = 0.0;

  for (j = 1; j < DISKGRID; j++) {

    expindex = p_engy[j] * expterm;
    prodecd = eengy * p_engy[j];

    if (expindex >= 150.0) {
      yofep[j] = 0.0;
      continue;

    } else if ((expindex >= 24.0) && (expindex < 150.0)) {
      yofep[j] = p_engy[j] * crossseciso(prodecd) / exp(expindex);

    } else if ((expindex < 24.0) && (expindex > 1.0e-4)) {
      yofep[j] = p_engy[j] * crossseciso(prodecd) / (exp(expindex) - 1.0);

    } else {
      yofep[j] = crossseciso(prodecd) / expterm;

    }

    if (yofep[j] < 1.0e-100) {
      yofep[j] = 0.0;
      continue;
    }

    sum += (0.5 * (yofep[j] + yofep[j - 1]) * (p_engy[j] - p_engy[j - 1]));

  }

  return(sum);

}


double gammaecdiso(double eengy, double Gamma_bulk, double zht, double *thetar, double *raddist, double *plab_engy)
{
  double sum = 0.0, x, rzterm, xterm, etaph, xBzterm, denom, expterm;
  double yofr[RGRID], gammadot, Beta_bulk, eave, p_engy[DISKGRID];
  int j, i;

  Beta_bulk = beta(Gamma_bulk);

  rzterm = 0.5 * SQR((raddist[0] / zht));
  if (rzterm < 1.0e-6) {    //Some small number to enable binomial expansion.
    x = zht * (1.0 + rzterm);
    etaph = (1.0 - (Beta_bulk * rzterm) - Beta_bulk) /
      (1.0 + rzterm - Beta_bulk);

  } else {
    xterm = SQR(raddist[0]) + SQR(zht);
    x = sqrt(xterm);
    etaph = (zht - (Beta_bulk * x)) / (x - (Beta_bulk * zht));

  }

  xBzterm = x - (Beta_bulk * zht);
  denom = Gamma_bulk * thetar[0] * xBzterm;
  expterm = x / denom;

  if (fabs(etaph) < 0.99999999) {

    for (i = 0; i < DISKGRID; i++) {

      p_engy[i] = plab_engy[i] / (1.0 + (Beta_bulk * etaph));

    }

    yofr[0] = raddist[0] * Int_epecdiso(eengy, expterm, p_engy)/ SQR(xBzterm);

  } else yofr[0] = 0.0;

  for (j = 1; j < RGRID; j++) {

    rzterm = 0.5 * SQR((raddist[j] / zht));
    if (rzterm < 1.0e-6) {
      x = zht * (1.0 + rzterm);
      etaph = (1.0 - (Beta_bulk * rzterm) - Beta_bulk) /
	(1.0 + rzterm - Beta_bulk);

    } else {
      xterm = SQR(raddist[j]) + SQR(zht);
      x = sqrt(xterm);
      etaph = (zht - (Beta_bulk * x)) / (x - (Beta_bulk * zht));

    }

    xBzterm = x - (Beta_bulk * zht);
    denom = Gamma_bulk * thetar[j] * xBzterm;
    expterm = x / denom;

    if (fabs(etaph) < 0.99999999) {

      for (i = 0; i < DISKGRID; i++) {

	p_engy[i] = plab_engy[i] / (1.0 + (Beta_bulk * etaph));

      }

      yofr[j] = raddist[j] * Int_epecdiso(eengy, expterm, p_engy)/ SQR(xBzterm);

    } else yofr[j] = 0.0;

    sum += (0.5 * (yofr[j] + yofr[j - 1]) * (raddist[j] - raddist[j - 1]));

  }

  gammadot = ECDISOCOEFF * sum / SQR(Gamma_bulk);

  return(-gammadot);

}


// Start of ECD ndot calculation using Jones formula
//
double gesgiso(double scatterphotengy, double pengy, double eengy)
{
  double engyratio, t1, eLower, eUpper, gesegam_iso;

  t1 = 4.0 * SQR(eengy) * pengy;
  eLower = pengy / (4.0 * SQR(eengy));
  eUpper = t1 / (1.0 + (4.0 * pengy * eengy));
  engyratio = scatterphotengy / (t1 * (1.0 - (scatterphotengy / eengy)));

  if ((scatterphotengy >= eLower) && (scatterphotengy <= pengy)) {

    gesegam_iso = (0.25 * scatterphotengy / SQR((eengy * pengy))) -
      (6.25e-02 / (SQR((eengy * eengy)) * pengy));

  } else if ((scatterphotengy > pengy) && (scatterphotengy <= eUpper)) {

    gesegam_iso = (1.0 / (pengy * SQR(eengy))) *
      ((0.5 * engyratio * log(engyratio)) +
       (0.25 * (1.0 + (2.0 * engyratio)) * (1.0 - engyratio)) +
       (0.125 * SQR(scatterphotengy) * (1.0 - engyratio) / (eengy * (eengy - scatterphotengy))));

  } else gesegam_iso = 0.0;

  return(gesegam_iso);

}


double epintiso(double eps, double eengy, double expterm, double *p_engy)
{
  double sum = 0.0, yofep[DISKGRID], expindex;
  int j;

  expindex = p_engy[0] * expterm;

  if (expindex >= 150.0) {
    yofep[0] = 0.0;

  } else if ((expindex >= 24.0) && (expindex < 150.0)) {
    yofep[0] = SQR(p_engy[0]) * gesgiso(eps, p_engy[0], eengy) / exp(expindex);

  } else if ((expindex < 24.0) && (expindex > 1.0e-4)) {
    yofep[0] = SQR(p_engy[0]) * gesgiso(eps, p_engy[0], eengy) /
      (exp(expindex) - 1.0);

  } else {
    yofep[0] = p_engy[0] * gesgiso(eps, p_engy[0], eengy) / expterm;

  }

  if (yofep[0] < 1.0e-100) yofep[0] = 0.0;

  for (j = 1; j < DISKGRID; j++) {

    expindex = p_engy[j] * expterm;

    if (expindex >= 150.0) {
      yofep[j] = 0.0;
      continue;

    } else if ((expindex >= 24.0) && (expindex < 150.0)) {
      yofep[j] = SQR(p_engy[j]) * gesgiso(eps, p_engy[j], eengy) / exp(expindex);

    } else if ((expindex < 24.0) && (expindex > 1.0e-4)) {
      yofep[j] = SQR(p_engy[j]) * gesgiso(eps, p_engy[j], eengy) /
	(exp(expindex) - 1.0);

    } else {
      yofep[j] = p_engy[j] * gesgiso(eps, p_engy[j], eengy) / expterm;

    }

    if (yofep[j] < 1.0e-100) {
      yofep[j] = 0.0;
      continue;
    }

    sum += (0.5 * (yofep[j] + yofep[j - 1]) * (p_engy[j] - p_engy[j - 1]));

  }

  return(sum);

}


double Int_gamecdiso(double scatterphotengy, double expterm, double *e_engy, double *nep, double *photengy)
{
  double sum = 0.0, funofgam[EGRID];
  int j;

  if (e_engy[0] < scatterphotengy) funofgam[0] = 0.0;
  else funofgam[0] = nep[0] * epintiso(scatterphotengy, e_engy[0], expterm, photengy);

  for (j = 1; j < EGRID; j++) {

    if (nep[j] == 0.0) {

      funofgam[j] = 0.0;
      continue;

    }

    if (e_engy[j] < scatterphotengy) funofgam[j] = 0.0;
    else funofgam[j] = nep[j] * epintiso(scatterphotengy, e_engy[j], expterm, photengy);

    sum += (0.5 * (funofgam[j] + funofgam[j - 1]) * (e_engy[j] - e_engy[j - 1]));

  }

  return(sum);

}


double ndotecdiso(double eps, double Gamma_bulk, double zht, double *thetar, double *raddist, double *e_engy, double *nep, double *plab_engy)
{
  double sum = 0.0, x, rzterm, xterm, etaph, finaltermiso;
  double yofr[RGRID], gammadot, Beta_bulk, eave, p_engy[DISKGRID];
  double xBzterm, denom, expterm;
  int j, i;

  Beta_bulk = beta(Gamma_bulk);

  rzterm = 0.5 * SQR((raddist[0] / zht));
  if (rzterm < 1.0e-6) {    //Some small number to enable binomial expansion.
    x = zht * (1.0 + rzterm);
    etaph = (1.0 - (Beta_bulk * rzterm) - Beta_bulk) /
      (1.0 + rzterm - Beta_bulk);

  } else {
    xterm = SQR(raddist[0]) + SQR(zht);
    x = sqrt(xterm);
    etaph = (zht - (Beta_bulk * x)) / (x - (Beta_bulk * zht));

  }

  xBzterm = x - (Beta_bulk * zht);
  denom = Gamma_bulk * thetar[0] * xBzterm;
  expterm = x / denom;

  if (fabs(etaph) < 0.99999999) {

    for (i = 0; i < DISKGRID; i++) {

      p_engy[i] = plab_engy[i] / (1.0 + (Beta_bulk * etaph));

    }

    yofr[0] = raddist[0] * Int_gamecdiso(eps, expterm, e_engy, nep, p_engy)/ SQR(xBzterm);

  } else yofr[0] = 0.0;

  for (j = 1; j < RGRID; j++) {

    rzterm = 0.5 * SQR((raddist[j] / zht));
    if (rzterm < 1.0e-6) {
      x = zht * (1.0 + rzterm);
      etaph = (1.0 - (Beta_bulk * rzterm) - Beta_bulk) /
	(1.0 + rzterm - Beta_bulk);

    } else {
      xterm = SQR(raddist[j]) + SQR(zht);
      x = sqrt(xterm);
      etaph = (zht - (Beta_bulk * x)) / (x - (Beta_bulk * zht));

    }

    xBzterm = x - (Beta_bulk * zht);
    denom = Gamma_bulk * thetar[j] * xBzterm;
    expterm = x / denom;

    if (fabs(etaph) < 0.99999999) {

      for (i = 0; i < DISKGRID; i++) {

	p_engy[i] = plab_engy[i] / (1.0 + (Beta_bulk * etaph));

      }

      yofr[j] = raddist[j] * Int_gamecdiso(eps, expterm, e_engy, nep, p_engy)/ SQR(xBzterm);

    } else yofr[j] = 0.0;

    sum += (0.5 * (yofr[j] + yofr[j - 1]) * (raddist[j] - raddist[j - 1]));

  }

  finaltermiso = (NECDISOCOEFF * sum) / SQR(Gamma_bulk);

  return(finaltermiso * SSCOEFF * 8.0);

}


 /////////////////////////////////////////////////
/*
double Int_phiecd(double gam, double eps, double etaph, double etas, double beta_gam, double expterm, int n, int rgrid)
{
  double sum = 0.0, phimin = 0.0, phimax = pi, steps, coski, Dterm, Mterm;
  double bmterm, bmobsterm, expindex, kappaprime, term1, bmepsterm, yterm;
  double yofphi[n];
  int j;

  steps = (phimax - phimin) / ((double) (n));
  bmobsterm = 1.0 - (beta_gam * etas);
  kappaprime = (etas - beta_gam) / bmobsterm;

  //printf("gam = %e, eps = %e, kappaprime = %e\n", gam, eps, kappaprime);

  coski = (etaph * etas) + (sqrt(((1.0 - SQR(etaph)) * (1.0 - SQR(etas)))) *
			 cos(phimin));
  if (coski < -.99999999) coski = -.99999999;
  if (coski > .99999999) coski = .99999999;

  //if (!(fabs(kappaprime) < .99999999)) yofphi[0] = 0.0;
  if (kappaprime < -.99999999) kappaprime = -.99999999;
  if (kappaprime > .99999999) kappaprime = .99999999;

  //else {
    bmterm = 1.0 - (beta_gam * coski);
    Mterm = gam * bmterm * (1.0 - beta_gam) * (1.0 + etas);
    bmepsterm = bmterm - (Mterm * eps);
    Dterm = bmobsterm / bmepsterm;
    expindex = Dterm * expterm;
    term1 = (bmepsterm / bmterm) + (bmterm / bmepsterm) - 1.0;
    yterm = term1 + SQR(kappaprime);

    if (bmterm < 1.0e-20) {
      yofphi[0] = 0.0;
      printf("For j = 0, eps = %e, gam = %e, bmterm = %e, rgrid = %d\n", eps, gam, bmterm, rgrid);
    }

    if (expindex > 150.0 || expindex < 1.0e-20) {
      yofphi[0] = 0.0;
      printf("For j = 0, eps = %e, gam = %e, expindex = %e, rgrid = %d\n", eps, gam, expindex, rgrid);
    }

    if (expindex > 1.0e-4) {
      yofphi[0] = (yterm * bmobsterm) / (SQR(bmepsterm) * (exp(expindex) - 1.0));

    } else yofphi[0] = yterm / (bmepsterm * expterm);

    //printf("yofphi[0] = %e\n", yofphi[0]);

    if (yofphi[0] < 1.0e-100) yofphi[0] = 0.0;

    //}

  phimin += steps;

  //printf("etaph = %e, expterm = %e, beta_gam = %e\n", etaph, expterm, beta_gam);

  //printf("In Int_phi, j = 0, coski = %e, bmterm = %e, Dterm = %e, expindex = %e, yofphi[0] = %e, phimin = %e\n", coski, bmterm, Dterm, expindex, yofphi[0], phimin);

  for (j = 1; j < n; j++) {

    coski = (etaph * etas) + (sqrt(((1.0 - SQR(etaph)) * (1.0 - SQR(etas)))) *
			 cos(phimin));
    if (coski < -.99999999) coski = -.99999999;
    if (coski > .99999999) coski = .99999999;

    //if (!(fabs(kappaprime) < .99999999)) yofphi[j] = 0.0;
    //if (kappaprime < -.99999999) kappaprime = -.99999999;
    //if (kappaprime > .99999999) kappaprime = .99999999;

    //else {
    bmterm = 1.0 - (beta_gam * coski);
    if (bmterm < 1.0e-20) {
      phimin += steps;
      yofphi[j] = 0.0;
      printf("1. Coming out of for loop in Int_phi at j = %d, bmterm = %e, gam = %e, eps = %e, rgrid = %d\n", j, bmterm, gam, eps, rgrid);
      continue;
    }

    Mterm = gam * bmterm * (1.0 - beta_gam) * (1.0 + etas);
    bmepsterm = bmterm - (Mterm * eps);
    Dterm = bmobsterm / bmepsterm;
    expindex = Dterm * expterm;
    if (expindex > 150.0 || expindex < 1.0e-20) {
      phimin += steps;
      yofphi[j] = 0.0;
      printf("2. Coming out of for loop in Int_phi at j = %d, expindex = %e, gam = %e, eps = %e, rgrid = %d\n", j, expindex, gam, eps, rgrid);
      continue;
    }

    term1 = (bmepsterm / bmterm) + (bmterm / bmepsterm) - 1.0;
    yterm = term1 + SQR(kappaprime);

    if (expindex > 1.0e-4) {
      yofphi[j] = (yterm * bmobsterm) / (SQR(bmepsterm) * (exp(expindex) - 1.0));
    } else yofphi[j] = yterm / (bmepsterm * expterm);

    if (yofphi[j] < 1.0e-100) {
      phimin += steps;
      //printf("yofphi[%d] = %e\n", j, yofphi[j]);
      yofphi[j] = 0.0;
      printf("3. Coming out of for loop in Int_phi at j = %d, gam = %e, eps = %e, rgrid = %d\n", j, gam, eps, rgrid);
      continue;
    }

    //}

    sum += (0.5 * (yofphi[j] + yofphi[j - 1]) * steps);
    phimin += steps;

    if (sum == 0.0) {
      printf("gam = %e, eps = %e, kappaprime = %e, yofphi[%d] = %e, yofphi[%d] = %e\n", gam, eps, kappaprime, j, yofphi[j], (j-1), yofphi[j-1]);
    }
    //printf("In Int_phi, j = %d, coski = %e, bmterm = %e, Dterm = %e, expindex = %e, yofphi[j] = %e, phimin = %e\n", j, coski, bmterm, Dterm, expindex, yofphi[j], phimin);

  }

  if (sum == 0.0) printf("j = %d, sum = %e, gam = %e, eps = %e, kappaprime = %e\n", j, sum, gam, eps, kappaprime);

  return(2.0 * sum);

}


double Recd_integ(double gam, double eps, double muobsprime, double Gamma_bulk, double zht, double *thetar, double *raddist, int n)
{

  double x, xterm, expterm, denom, sum = 0.0, yofr[GS];
  double etaph, xBzterm, beta_gam, Beta_bulk, rzterm;
  int j;

  beta_gam = beta(gam);
  Beta_bulk = beta(Gamma_bulk);

  rzterm = 0.5 * SQR((raddist[0] / zht));
  if (rzterm < 1.0e-6) {
    x = zht * (1.0 + rzterm);
    etaph = (1.0 - (Beta_bulk * rzterm) - Beta_bulk) /
      (1.0 + rzterm - Beta_bulk);

  } else {
    xterm = SQR(raddist[0]) + SQR(zht);
    x = sqrt(xterm);
    etaph = (zht - (Beta_bulk * x)) / (x - (Beta_bulk * zht));

  }

  xBzterm = x - (Beta_bulk * zht);
  denom = Gamma_bulk * thetar[0] * xBzterm;
  expterm = (eps * x) / denom;

  if (fabs(etaph) < .99999999) {
    yofr[0] = raddist[0]*Int_phiecd(gam, eps, etaph, muobsprime, beta_gam, expterm, n, 0) /
      SQR(xBzterm);
  } else yofr[0] = 0.0;

  //printf("eps = %e, gam = %e\n", eps, gam);
  //printf("Int_phifun = %e\n", Int_phiecd(gam, eps, etaph, muobsprime, beta_gam, expterm, n));
  //printf("In Recd_integ, for j = 0, rzterm = %e, x = %e, etaph = %e, expterm = %e, xBzterm = %e, yofr[0] = %e\n", rzterm, x, etaph, expterm, xBzterm, yofr[0]);

  for (j = 1; j < GS; j++) {

    rzterm = 0.5 * SQR((raddist[j] / zht));
    if (rzterm < 1.0e-6) {
      x = zht * (1.0 + rzterm);
      etaph = (1.0 - (Beta_bulk * rzterm) - Beta_bulk) /
	(1.0 + rzterm - Beta_bulk);
    } else {
      xterm = SQR(raddist[j]) + SQR(zht);
      x = sqrt(xterm);
      etaph = (zht - (Beta_bulk * x)) / (x - (Beta_bulk * zht));
    }

    xBzterm = x - (Beta_bulk * zht);
    denom = Gamma_bulk * thetar[j] * xBzterm;
    expterm = (eps * x) / denom;

    if (fabs(etaph) < .99999999) {
      yofr[j] = raddist[j]*Int_phiecd(gam, eps, etaph, muobsprime, beta_gam, expterm, n, j) /
	SQR(xBzterm);
  } else yofr[j] = 0.0;

    //printf("Int_phifun = %e\n", Int_phiecd(gam, eps, etaph, muobsprime, beta_gam, expterm, n));
    //printf("In Recd_integ, j = %d, rzterm = %e, x = %e, etaph = %e, expterm = %e, xBzterm = %e, yofr[j] = %e\n", j, rzterm, x, etaph, expterm, xBzterm, yofr[j]);

    sum += (0.5 * (yofr[j] + yofr[j - 1]) * (raddist[j] - raddist[j - 1]));

  }

  return(sum);

}

double ndotecd(double eps, double Gamma_bulk, double zht, double muobsprime, double *thetar, double *raddist, double *e_engy, double *nep, int n)
{

  double finalterm, M, D, gamterm, sum = 0.0, yofgam[n];
  int j;

  if (!(fabs(muobsprime) < .99999999)) {
    printf ("\n muobsprime = %e !! ", muobsprime);
    return (0.0);
  }

  if (nep[0] == 0.0) yofgam[0] = 0.0;
  else {
    gamterm = nep[0] / SQR(e_engy[0]);
    yofgam[0] = gamterm * Recd_integ(e_engy[0], eps, muobsprime, Gamma_bulk, zht, thetar, raddist, n);
  }

  for (j = 1; j < n; j++) {

    if (nep[j] == 0.0) {
      yofgam[j] = 0.0;
      continue;
    }

    gamterm = nep[j] / SQR(e_engy[j]);
    yofgam[j] = gamterm * Recd_integ(e_engy[j], eps, muobsprime, Gamma_bulk, zht, thetar, raddist, n);

    sum += (0.5 * (yofgam[j] + yofgam[j - 1]) * (e_engy[j] - e_engy[j - 1]));

  }

  finalterm = SQR((eps / Gamma_bulk)) * sum;
  printf("eps = %e, ndotecd = %e\n", eps, (NDOTCOEFF * finalterm));

  return(NDOTCOEFF * finalterm);

}
*/
