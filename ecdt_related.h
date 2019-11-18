#define DTPHOTCOEFF 1.40311e29 // Value of 2*(m_e*c/h)^3   
#define GAMDTTHCOEFF 1.905e-10 // Value of 2*pi*c*sigma_T*pi^4*THETA_DT^4*DTPHOT/15
#define NDOTDTTHCOEFF 2.79932e15 // Value of c*sigma_T*DTPHOT
#define NDOTDTISOCOEFF 6.28318 // Value of 2*pi
#define THETADT 2.02089e-7 // Value of THETA_DT for T = 1200 K      
#define RAVGDT 2.954e18  // Value of lab frame R_DT = 0.957 pc in cm  
#define DTDIST 6.040e18  // Value of lab frame DT_dist = 1.957 pc in cm
#define ETADT 30
#define DTGRID 60
//#define FCOVDT 1.0


// Start of calculation for differential photon number density from Dusty 
// Torus. Calculating the angle at which the dusty torus photons enter the jet
// in the comoving frame of the jet. 
//                                                            
void etadet(double BLF, double emitdist, double *mincos, double *maxcos, double *etadtarr)
{

  double rtterm, Rdtterm, x, Beta_BLF, zterm, termin, termax, term1, etastep;
  double etadtmin, etadtmax, dtetmin, etadtmaxenter, etadtminenter;
  double etaminstar, etamaxstar, etamidminusmax, etamidplusmin;
  double etastepminus, etastepplus, dtetminplus;
  int i, count = 0, negval = 0, negplus = 0;

  Beta_BLF = beta(BLF);
  Rdtterm = SQR(DTDIST) - SQR(RAVGDT);
  rtterm = DTDIST * RAVGDT;
  zterm = SQR(emitdist) + SQR(DTDIST);

  x = Rdtterm / SQR(emitdist);

  if (x < 2.0e-6) {

    term1 = 0.5 * ((2.0 * SQR(emitdist)) + Rdtterm);

  } else {

    term1 = SQR(emitdist) * sqrt(1.0 + x);

  }

  termin = term1 - rtterm - (Beta_BLF * zterm);
  etaminstar = (term1 + rtterm) / zterm;

  *mincos = termin / (zterm - (Beta_BLF * (term1 - rtterm)));

  if (*mincos < -.99999999) *mincos = -.99999999;
  if (*mincos > .99999999) *mincos = .99999999;  

  termax = term1 + rtterm - (Beta_BLF * zterm);
  etamaxstar = (term1 - rtterm) / zterm;

  *maxcos = termax / (zterm - (Beta_BLF * (term1 + rtterm)));

  if (*maxcos < -.99999999) *maxcos = -.99999999;
  if (*maxcos > .99999999) *maxcos = .99999999;

  etadtmin = *mincos;
  etadtmax = *maxcos;

  if ((etadtmin < 0.0) && (etadtmax < 0.0)) {

    etastep = pow((fabs(etadtmin) / fabs(etadtmax)), 1.0/((double) (ETADT - 1)));
    //printf("etastep = %e\n", etastep);

    dtetmin = etadtmin;

    negval = 1;
    negplus = 0;

  } else if ((etadtmin > 0.0) && (etadtmax > 0.0)) {

    etastep = pow((fabs(etadtmax) / fabs(etadtmin)), 1.0/((double) (ETADT - 1)));
    //printf("etastep = %e\n", etastep);

    dtetmin = etadtmin;

    negval = 0;
    negplus = 0;

  } else if ((etadtmin < 0.0) && (etadtmax > 0.0)) {

    etamidminusmax = -1.0e-3;
    etamidplusmin = 1.0e-3;

    etastepminus = pow((fabs(etadtmin) / fabs(etamidminusmax)), 1.0/((double) (ETADT - (ETADT/2) - 1)));
    etastepplus = pow((fabs(etadtmax) / fabs(etamidplusmin)), 1.0/((double) (ETADT - (ETADT/2) - 1)));

    dtetmin = etadtmin;
    dtetminplus = etamidplusmin;

    negplus = 1;
    negval = 0;
  }

  /*
  printf("etaminstar = %e, etamaxstar = %e, etadtmin = %e, etadtmax = %e\n",
         etaminstar, etamaxstar, etadtmin, etadtmax);
  */

  if ((negval == 1) && (negplus == 0)) {
    
    for (i = 0; i < ETADT; i++) {
      
      etadtarr[i] = dtetmin;
      
      //To ensure that for negative values the ETA grid proceeds properly.

      dtetmin = dtetmin / etastep;  

    }
    
  } else if ((negval == 0) && (negplus == 1)) {

    for (i = 0; (i < ETADT/2); i++) {

      etadtarr[i] = dtetmin;

      dtetmin = dtetmin / etastepminus;

    }

    for (i = (ETADT/2); i < ETADT; i++) {

      etadtarr[i] = dtetminplus;

      dtetminplus *= etastepplus;

    }

  } else {

    for (i = 0; i < ETADT; i++) {

      etadtarr[i] = dtetmin;

      dtetmin *= etastep;

    }

  }

  /*
  for (i = 0; i < ETADT; i++) {

    printf("etadtarr[%d] = %e\n", i, etadtarr[i]);

  }
  */

  return;

}


void nphdt(double BLF, double *etadtarr, double *dtlabpe, double dtep[][DTGRID], double photden_dt[][DTGRID], double FCOVDT) 
{

  double Beta_BLF, bphterm, expterm, contdenom;
  int i, k;

  Beta_BLF = beta(BLF);

  for (i = 0; i < ETADT; i++) {

    bphterm = 1.0 + (Beta_BLF * etadtarr[i]);

    for (k = 0; k < DTGRID; k++) {
      
      dtep[i][k] = dtlabpe[k] / (BLF * bphterm);

      expterm = dtlabpe[k] / THETADT;

      if (expterm >= 150.0) {
	photden_dt[i][k] = 0.0;
	continue;

      } else if ((expterm >= 24.0) && (expterm < 150.0)) {
	contdenom = exp(expterm);
	photden_dt[i][k] = DTPHOTCOEFF * FCOVDT * SQR(dtep[i][k]) / contdenom;

      } else if ((expterm < 24.0) && (expterm > 1.0e-4)) {
	contdenom = exp(expterm) - 1.0;
	photden_dt[i][k] = DTPHOTCOEFF * FCOVDT * SQR(dtep[i][k]) / contdenom;

      } else {
	contdenom = expterm;
	photden_dt[i][k] = DTPHOTCOEFF * FCOVDT * SQR(dtep[i][k]) / contdenom;

      }

      if (photden_dt[i][k] < 1.0e-100) {
	
	photden_dt[i][k] = 0.0;

      }

    }

  }

  return;
  
}


// Start of gammadot calculation resulting from inverse Compton scattering of 
// Dusty Torus photons in the emission region, calculated in the PF. 
// The formula is based on head approx. given in Dermer & Menon book, 2007, 
// Ch. 6.                                                           
//         
double nphepsumdt(double eengy, double phcosdt, double *photengydt, double *photdendt)
{
  double sum = 0.0, M, D, beta_e, term1, term2, term3, num1, den1, num2, den2;
  double num3, den3, funofep[DTGRID], bphterm, lnterm, D2term, D3term, gpe;
  int j;

  beta_e = beta(eengy);
  bphterm = 1.0 - (beta_e * phcosdt);

  gpe = eengy  - photengydt[0];
  M = eengy * photengydt[0] * bphterm;
  D = 1.0 + (2.0 * M);

  num1 = (M * (M - 2.0)) - 2.0;
  den1 = eengy * photengydt[0] * SQR(M);

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

  den2 = 3.0 * photengydt[0] * D3term;
  num2 = 6.0 * photengydt[0] * D * gpe * (1.0 + M) * bphterm;
  term2 = 1.0 - D3term + num2;

  num3 = (2.0 * gpe * D) - (eengy * (M * (M - 1.0) - 1.0));
  den3 = eengy * M;
  term3 = 6.0 * D2term * num3 / den3;

  funofep[0] = photdendt[0] * (term1 + ((term2 + term3) / den2));

  if (funofep[0] < 0.0) {
    funofep[0] = 0.0;
  }

  for (j = 1; j < DTGRID; j++) {

    gpe = eengy  - photengydt[j];
    M = eengy * photengydt[j] * bphterm;
    D = 1.0 + (2.0 * M);

    num1 = (M * (M - 2.0)) - 2.0;
    den1 = eengy * photengydt[j] * SQR(M);

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

    den2 = 3.0 * photengydt[j] * D3term;
    num2 = 6.0 * photengydt[j] * D * gpe * (1.0 + M) * bphterm;
    term2 = 1.0 - D3term + num2;

    num3 = (2.0 * gpe * D) - (eengy * (M * (M - 1.0) - 1.0));
    den3 = eengy * M;
    term3 = 6.0 * D2term * num3 / den3;

    funofep[j] = photdendt[j] * (term1 + ((term2 + term3) / den2));

    if (funofep[j] < 0.0) {
      funofep[j] = 0.0;
      continue;
    }

    sum += (0.5 * (funofep[j] + funofep[j - 1]) *
            (photengydt[j] - photengydt[j - 1]));

  }

  return(sum);

}


double gamdotdt(double eengy, double *etadtarr, double pedt[][DTGRID], double pddt[][DTGRID])
{

  double sum = 0.0, gammadotdt, yofdteta[ETADT];
  int j;

  yofdteta[0] = nphepinteg(eengy, etadtarr[0], pedt[0], pddt[0], DTGRID);

  for (j = 1; j < ETADT; j++) {

    yofdteta[j] = nphepinteg(eengy, etadtarr[j], pedt[j], pddt[j], DTGRID);

    sum += (0.5 * (yofdteta[j] + yofdteta[j - 1]) * 
	    (etadtarr[j] - etadtarr[j - 1]));

  }

  gammadotdt = GDOTBLRCOEFF * sum;

  return(-gammadotdt);

}


// Start of ECDT ndot, in PF, calculation using Compton cross-section under 
// head-on approx. based on eqn. derived using Dermer & Menon book, Ch-6.  
// Int_gamblr function taken from ecblr_related.h file as the gamma integration
// was exactly the same.
//
double Int_epdt(double coskidt, double epsdt, double pedtmax, double *e_engy, double *nep, double *pengydt, double *np_dt)
{
  double sum = 0.0, yofepdt[DTGRID], epULdt, epmiddt;
  int j, n = 0;

  epULdt = (2.0 * epsdt) / (1.0 - coskidt);

  yofepdt[0] = (np_dt[0] / pengydt[0]) *
    Int_gamblr(coskidt, epsdt, pengydt[0], e_engy, nep);

  if (yofepdt[0] < 1.0e-100) yofepdt[0] = 0.0;

  for (j = 1; j < DTGRID; j++) {

    if (pengydt[j] > epULdt) {

      break;
    }

    yofepdt[j] = (np_dt[j] / pengydt[j]) *
      Int_gamblr(coskidt, epsdt, pengydt[j], e_engy, nep);

    if (yofepdt[j] < 1.0e-100) {
      yofepdt[j] = 0.0;
      continue;
    }

    sum += (0.5 * (yofepdt[j] + yofepdt[j - 1]) * 
	    (pengydt[j] - pengydt[j - 1]));

  }

  return(sum);

}


double etadt_integ(double eps, double etas, double phivaldt, double *e_engy, double *nep, double *etadtarr, double photengydt[][DTGRID], double photdendt[][DTGRID])
{
  double sum = 0.0, yofetadt[ETADT], coskidt;
  int j, k = 0;

  coskidt = (etadtarr[0] * etas) + (sqrt(((1.0 - SQR(etadtarr[0])) * 
					  (1.0 - SQR(etas)))) * cos(phivaldt));

  if (coskidt < -0.99999999) coskidt = -0.99999999;
  if (coskidt > 0.99999999) coskidt = 0.99999999;

  yofetadt[0] = Int_epdt(coskidt, eps, photengydt[0][DTGRID - 1], e_engy, nep, photengydt[0], photdendt[0]);

  for (j = 1; j < ETADT; j++) {

    coskidt = (etadtarr[j] * etas) + 
      (sqrt(((1.0 - SQR(etadtarr[j])) * (1.0 - SQR(etas)))) * cos(phivaldt));
    
    if (coskidt < -0.99999999) coskidt = -0.99999999;
    if (coskidt > 0.99999999) coskidt = 0.99999999;

    yofetadt[j] = Int_epdt(coskidt, eps, photengydt[j][DTGRID - 1], e_engy, nep, photengydt[j], photdendt[j]);

    sum += (0.5 * (yofetadt[j] + yofetadt[j - 1]) * 
	    (etadtarr[j] - etadtarr[j - 1]));

  }

  return(sum);

}


double ndotdt(double eps, double muobsprime, double *e_engy, double *nep, double *etadtarr, double *phidt, double pedt[][DTGRID], double pddt[][DTGRID], double phistep)
{
  double finaltermdt, sum = 0.0, yofphi[PHIGRID];
  int j;

  if (!(fabs(muobsprime) < .99999999)) {
    printf ("\n muobsprime = %e !! ", muobsprime);
    return (0.0);
  }

  yofphi[0] = etadt_integ(eps, muobsprime, phidt[0], e_engy, nep, etadtarr, pedt, pddt);

  for (j = 1; j < PHIGRID; j++) {

    yofphi[j] = etadt_integ(eps, muobsprime, phidt[j], e_engy, nep, etadtarr, pedt, pddt);

    sum += (0.5 * (yofphi[j] + yofphi[j - 1]) * phistep);

  }

  finaltermdt = SSCOEFF * sum;

  return(2.0 * finaltermdt);

}


// Start of gammadot calculation resulting from inverse Compton scattering of 
// DT photons in the emission region, in the PF, in Thomson regime.       
//    
double gamdotdtth(double eengy, double BLF, double *etadtarr, double FCOVDT)
{
  double yofetath[ETADT], sum = 0.0, beta_e, bphterm, BLFterm, Beta_BLF;
  double denom, term1, gammadotth;
  int j;

  Beta_BLF = beta(BLF);
  beta_e = beta(eengy);

  BLFterm = 1.0 + (Beta_BLF * etadtarr[0]);
  denom = SQR(BLFterm) * SQR(BLFterm);

  bphterm = 1.0 - (beta_e * etadtarr[0]);
  term1 = bphterm * ((SQR(eengy) * bphterm) - 1.0);

  yofetath[0] = term1 / denom;

  if (yofetath[0] < 0.0) {
    yofetath[0] = 0.0;
  }

  for (j = 1; j < ETADT; j++) {

    BLFterm = 1.0 + (Beta_BLF * etadtarr[j]);
    denom = SQR(BLFterm) * SQR(BLFterm);

    bphterm = 1.0 - (beta_e * etadtarr[j]);
    term1 = bphterm * ((SQR(eengy) * bphterm) - 1.0);

    yofetath[j] = term1 / denom;

    if (yofetath[j] < 0.0) {
      yofetath[j] = 0.0;
      continue;
    }

    sum += (0.5 * (yofetath[j] + yofetath[j - 1]) *
	    (etadtarr[j] - etadtarr[j - 1]));

  }

  gammadotth = GAMDTTHCOEFF * FCOVDT * sum / (SQR(BLF) * SQR(BLF));

  return(-gammadotth);

}


// Start of ndot calculation resulting from inverse Compton scattering of 
// DT photons in the emission region, in the PF, in Thomson regime.      
// Int_gamblrth function taken from ecblr_related.h file.
//                                                      
double etadtth_integ(double phivaldt, double epstempdt, double BLF, double etas, double *etadtarr, double *e_engy, double *nep)
{
  double sum = 0.0, yofetadtth[ETADT], Beta_bulk, bphtermdt, etermfacdt, coskidt;
  int j;

  Beta_bulk = beta(BLF);

  bphtermdt = 1.0 + (Beta_bulk * etadtarr[0]);
  etermfacdt = epstempdt * BLF * bphtermdt;

  coskidt = (etadtarr[0] * etas) + (sqrt(((1.0 - SQR(etadtarr[0])) *
                                          (1.0 - SQR(etas)))) * cos(phivaldt));

  if (coskidt < -0.99999999) coskidt = -0.99999999;
  if (coskidt > 0.99999999) coskidt = 0.99999999;

  yofetadtth[0] = Int_gamblrth(coskidt, etermfacdt, e_engy, nep);

  if ((yofetadtth[0] < 0.0) || (yofetadtth[0] < 1.0e-100)) {
    yofetadtth[0] = 0.0;
  }

  for (j = 1; j < ETADT; j++) {

    bphtermdt = 1.0 + (Beta_bulk * etadtarr[j]);
    etermfacdt = epstempdt * BLF * bphtermdt;

    coskidt = (etadtarr[j] * etas) + (sqrt(((1.0 - SQR(etadtarr[j])) *
					    (1.0 - SQR(etas)))) * cos(phivaldt));

    if (coskidt < -0.99999999) coskidt = -0.99999999;
    if (coskidt > 0.99999999) coskidt = 0.99999999;

    yofetadtth[j] = Int_gamblrth(coskidt, etermfacdt, e_engy, nep);

    if ((yofetadtth[j] < 0.0) || (yofetadtth[j] < 1.0e-100)) {
      yofetadtth[j] = 0.0;
      continue;
    }

    sum += (0.5 * (yofetadtth[j] + yofetadtth[j - 1]) * (etadtarr[j] - etadtarr[j - 1]));

  }

  return(sum);

}


double ndotdtth(double eps, double Gamma_bulk, double muobsprime, double *phidt, double *etadtarr, double *e_engy, double *nep, double phistep, double FCOVDT)
{
  double finaltermdtth, sum = 0.0, yofphidt[PHIGRID], epstempdt;
  int i, j;

  if (!(fabs(muobsprime) < .99999999)) {
    printf ("\n muobsprime = %e !! ", muobsprime);
    return (0.0);
  }

  epstempdt = eps / THETADT;

  yofphidt[0] = etadtth_integ(phidt[0], epstempdt, Gamma_bulk, muobsprime, etadtarr, e_engy, nep);

  for (j = 1; j < PHIGRID; j++) {

    yofphidt[j] = etadtth_integ(phidt[j], epstempdt, Gamma_bulk, muobsprime, etadtarr, e_engy, nep);

    sum += (0.5 * (yofphidt[j] + yofphidt[j - 1]) * phistep);

  }

  finaltermdtth = NDOTDTTHCOEFF * FCOVDT * SQR(eps) * sum;

  return(2.0 * finaltermdtth);

}


// Start of ECDT gammadot calculation using Boettcher et al. 1997 formula
//
double Int_epdtiso(double eengy, double *dtcom_engy, double *pddtcom)
{
  double sum = 0.0, funofepdt[DTGRID], proddt;
  int j;

  proddt = eengy * dtcom_engy[0];

  funofepdt[0] = crossseciso(proddt) * pddtcom[0] / dtcom_engy[0];

  if (funofepdt[0] < 1.0e-100) {
    funofepdt[0] = 0.0;
  }

  for (j = 1; j < DTGRID; j++) {

    proddt = eengy * dtcom_engy[j];

    funofepdt[j] = crossseciso(proddt) * pddtcom[j] / dtcom_engy[j];

    if (funofepdt[j] < 1.0e-100) {
      funofepdt[j] = 0.0;
      continue;
    }

    sum += (0.5 * (funofepdt[j] + funofepdt[j - 1]) *
            (dtcom_engy[j] - dtcom_engy[j - 1]));

  }

  return(sum);
    
}


double gamdotdtiso(double eengy, double *etadtarr, double pedt[][DTGRID], double pddt[][DTGRID])
{
  
  double sum = 0.0, gammadotdtiso, yofetadt[ETADT];
  int j;
  
  yofetadt[0] = Int_epdtiso(eengy, pedt[0], pddt[0]);

  for (j = 1; j < ETADT; j++) {

    yofetadt[j] = Int_epdtiso(eengy, pedt[j], pddt[j]);

    sum += (0.5 * (yofetadt[j] + yofetadt[j - 1]) *
            (etadtarr[j] - etadtarr[j - 1]));

  }

  gammadotdtiso = GDOTBLRCOEFF * sum;

  return(-gammadotdtiso);

}


// Start of ECD ndot calculation using Jones formula                                                                                                                           
//
double epdtintiso(double eps, double eengy, double *dtcomengy, double *pddtcom)
{
  double sum = 0.0, yofepdtiso[DTGRID];
  int j;

  yofepdtiso[0] = pddtcom[0] * gesgiso(eps, dtcomengy[0], eengy);

  if (yofepdtiso[0] < 1.0e-100) yofepdtiso[0] = 0.0;

  for (j = 1; j < DTGRID; j++) {

    yofepdtiso[j] = pddtcom[j] * gesgiso(eps, dtcomengy[j], eengy);

    if (yofepdtiso[j] < 1.0e-100) {
      yofepdtiso[j] = 0.0;
      continue;
    }

    sum += (0.5 * (yofepdtiso[j] + yofepdtiso[j - 1]) * 
	    (dtcomengy[j] - dtcomengy[j - 1]));

  }

  return(sum);

}


double etadtintiso(double eps, double eengy, double *etadtarr, double dtcom_engy[][DTGRID], double pddt_com[][DTGRID])
{
  double sum = 0.0, yofetadtiso[ETADT];
  int j;

  yofetadtiso[0] = epdtintiso(eps, eengy, dtcom_engy[0], pddt_com[0]);

  for (j = 1; j < ETADT; j++) {

    yofetadtiso[j] = epdtintiso(eps, eengy, dtcom_engy[j], pddt_com[j]);

    sum += (0.5 * (yofetadtiso[j] + yofetadtiso[j - 1]) *
            (etadtarr[j] - etadtarr[j - 1]));

  }

  return(sum);
  
}


double ndotdtiso(double eps, double *e_engy, double *nep, double *etadtarr, double dtcom_engy[][DTGRID], double pddt_com[][DTGRID])
{
  double sum = 0.0, nphdotdtiso, funofgam[EGRID];
  int j;
  
  if ((e_engy[0] < eps) || (nep[0] == 0.0)) {

    funofgam[0] = 0.0;

  } else {

    funofgam[0] = nep[0] * etadtintiso(eps, e_engy[0], etadtarr, dtcom_engy, pddt_com);

  }

  for (j = 1; j < EGRID; j++) {

    if ((e_engy[j] < eps) || (nep[j] == 0.0)) {

      funofgam[j] = 0.0;
      continue;

    } else {

      funofgam[j] = nep[j] * etadtintiso(eps, e_engy[j], etadtarr, dtcom_engy, pddt_com);

    }

    sum += (0.5 * (funofgam[j] + funofgam[j - 1]) * (e_engy[j] - e_engy[j - 1]));

  }

  nphdotdtiso = NDOTDTISOCOEFF * SSCOEFF * 8.0 * sum;

  return(nphdotdtiso);

}
