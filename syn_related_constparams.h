#define C1 1.088946                           // 1.2042 * 0.90429
#define C2 2.358614e-3                        // 2.60825e-3 * 0.90429
#define P1 0.20949
#define P2 0.79051
#define NU_C 4.19597e+06		      // Coefficient for 
                                              // characteristic syn.
					      // freq.(3*e/(4*pi*m_e*c))
#define SYNCOEFF 1.29095e-09                  // Value of 4*c*sigma_t/(24*pi*m_e*c^2)
#define NSYNCOEFF 3.47882e+23                 // Value of sqrt(3)*e^3/(4*pi*h^2)
#define XCOEFF 2.38324e-07                    // Value of 4*pi*m_e*c/(3*e)
#define TAUSYNCOEFF 1.02257e+04               // Value of sqrt(3)*e^3/(8*pi*(m_e*c)^2) 


double ndotsyn(double nu, double nu_1, double mf, double *e_engy, double *nep, int n)
{
  double x1, Rx_app, sum = 0.0, t1, t2, factor, nu_norm, x1_param, x1nu_param, yofx[EGRID];
  int j;

  nu_norm = 0.01 * nu_1;
  factor = NSYNCOEFF * mf / nu;  

  if (nu < nu_norm) {

      x1 = XCOEFF * nu_norm / (mf * SQR(e_engy[0]));

    } else {
      
      x1 = XCOEFF * nu / (mf * SQR(e_engy[0]));

    }

    t1 = C1 * pow(x1, P1);
    t2 = C2 / pow(x1, P2);

    Rx_app = (t1 - t2) * exp(-x1); 
    if (nu < nu_norm) {
      
      Rx_app = pow((nu/nu_norm), (1.0/3.0)) * Rx_app;
      
    }

    if (Rx_app > 1.0e-52) {

      yofx[0] = Rx_app * nep[0];
    
    } else yofx[0] = 0.0;

  x1_param = XCOEFF * nu_norm / mf;
  x1nu_param = XCOEFF * nu / mf;

  for (j = 1; j < n; j++) {

    if (nep[j] == 0.0) {

      yofx[j] = 0.0;
      continue;

    }

    if (nu < nu_norm) {

      x1 = x1_param / SQR(e_engy[j]);

    } else {
      
      x1 = x1nu_param / SQR(e_engy[j]);

    }

    t1 = C1 * pow(x1, P1);
    t2 = C2 / pow(x1, P2);

    Rx_app = (t1 - t2) * exp(-x1); 
    if (nu < nu_norm) {
      
      Rx_app = pow((nu/nu_norm), (1.0/3.0)) * Rx_app;
      
    }

    if (Rx_app > 1.0e-52) {

      yofx[j] = Rx_app * nep[j];
    
    } else yofx[j] = 0.0;

    sum += (0.5 * (yofx[j] + yofx[j - 1]) * (e_engy[j] - e_engy[j - 1]));
    
  }

  return (factor * sum);

}



double alpha_ssa(double nu, double nu_1, double mf, double *e_engy, double *nep) {

  double sum = 0.0, dnepge, t1, t2, x1, Rx_app, factor, nu_norm, nepgarr[EGRID], funofx[EGRID];
  double x1_param, x1nu_param;
  int j;

  nu_norm = 0.01 * nu_1;

  factor = TAUSYNCOEFF * mf / SQR(nu);

  nepgarr[0] = nep[0] / SQR(e_engy[0]);
  funofx[0] = 0.0;

  x1_param = XCOEFF * nu_norm / mf;
  x1nu_param = XCOEFF * nu / mf;

  for (j = 1; j < EGRID; j++) {

    if (nep[j] == 0.0) {

      nepgarr[j] = 0.0;
      funofx[j] = 0.0;
      continue;

    }

    if (nu < nu_norm) {
      
      x1 = x1_param / SQR(e_engy[j]);

    } else {

      x1 = x1nu_param / SQR(e_engy[j]);

    }

    t1 = C1 * pow(x1, P1);
    t2 = C2 / pow(x1, P2);

    Rx_app = (t1 - t2) * exp(-x1); 
    if (nu < nu_norm) {

      Rx_app = pow((nu/nu_norm), (1.0/3.0)) * Rx_app;

    }

    nepgarr[j] = nep[j] / SQR(e_engy[j]);

    if (Rx_app > 1.0e-52) {

      funofx[j] = Rx_app * SQR(e_engy[j]) * (nepgarr[j] - nepgarr[j - 1]) / 
	(e_engy[j] - e_engy[j - 1]);

    } else funofx[j] = 0.0;

    sum += (0.5 * (funofx[j] + funofx[j - 1]) * (e_engy[j] - e_engy[j - 1]));

  }
  
  sum *= factor;

  return (-sum);

}



double gammasyndot(double gammax, double mf)
{
  double gammadot;

  gammadot = SYNCOEFF * SQR(mf) * (SQR(gammax) - 1.0);

  return(-gammadot);

}
