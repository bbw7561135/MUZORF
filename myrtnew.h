void funcal(double x, double* funl, double* dfun, double term4, double term3, double term2, double term1, double term0)
{
  
  double t1, t2, t3, t4, t5;
  double difft1, difft2, difft3, difft4;
  
  t1 = term4 * pow(x, 4.0);
  t2 = term3 * pow(x, 3.0);
  t3 = term2 * pow(x, 2.0);
  t4 = term1 * x;
  t5 = term0;
  
  *funl = t1 + t2 + t3 + t4 + t5;
  
  difft1 = 4.0 * term4 * pow(x, 3.0);
  difft2 = 3.0 * term3 * pow(x, 2.0);
  difft3 = 2.0 * term2 * x;
  difft4 = term1;
  
  *dfun = difft1 + difft2 + difft3 + difft4;
}



double myrtnewt(void (*funcd)(double, double *, double *, double, double, double, double, double), double po, double xacc, double term4, double term3, double term2, double term1, double term0)
{
      int j;
      double p, fpo, dfpo;

      for (j = 1; j <= MAXIT; j++) {
	(*funcd)(po, &fpo, &dfpo, term4, term3, term2, term1, term0);
	p = po - (fpo/dfpo);
	if (fabs(p - po) < xacc) return p;
	else po = p;
      }
      printf("\nMaximum number of iterations exceeded in myrtnewt\n");
      return 0.0;
}



void funcalq(double gmin, double *fung, double *dfung, double index, double coeff, double coeffq, double dmax, double bmax, double amax)
{
  
  if (index == 1) {
    
    *fung = (coeff * log(gmin)) - gmin + dmax;
    *dfung = (coeff / gmin) - 1.0;

  } else if (index == 2) {

    *fung = log(gmin) + (coeff / gmin) - bmax;
    *dfung = (1.0 / gmin) - (coeff / pow(gmin, 2.0));

  } else {

    *fung = pow(gmin, (2.0 - index)) - (coeffq * pow(gmin, (1.0 - index))) - amax;
    *dfung = ((2.0 - index) * pow(gmin, (1.0 - index))) - ((1.0 - index) * pow(gmin, (-index)) * coeffq);

  }

}



double myrtnewtq(void (*funcd)(double, double *, double *, double, double, double, double, double, double), double po, double xacc, double index, double coeff, double coeffq, double dmax, double bmax, double amax)
{
      int j;
      double p, fpo, dfpo;

      for (j = 1; j <= MAXIT; j++) {

           (*funcd)(po, &fpo, &dfpo, index, coeff, coeffq, dmax, bmax, amax);
           p = po - (fpo/dfpo);
           if (((fabs(p - po)) / p) < xacc) return p;
           else po = p;
      }
      printf("\nMaximum number of iterations exceeded in myrtnewtq\n");
      return 0.0;
}
