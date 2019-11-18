/* VEGAS code from numerical recipes in C, 3rd Edition */

#define ALPH 1.5
#define NDMX  50
#define MXDIM 10
#define TINY 1.0e-30
//#define SQR(x) (x*x)

long idum;                    //For random number initialization in main.


int IMAX(int x, int y)
{

if (x > y) return (x);
else return (y);

}



int IMIN(int x, int y)
{

if (x > y) return (y);
else return (x);

}


void vegas(double regn[], int ndim, double (*fxn)(double [], double, double), int init, unsigned long ncall, int itmx, int nprn, double *tgral, double *sd, double *chi2a, double arh)
{
  double ran2(long *idum);

  void rebin(double rc, int nd, double r[], double xin[], double xi[]);

  static int i, it, j, k, mds, nd, ndo, ng, npg, ia[MXDIM+1], kg[MXDIM+1];
  static double calls, dv2g, dxg, f, f2, f2b, fb, rc, ti, tsi, wgt, xjac,
    xn, xnd, xo;
  static double d[NDMX+1][MXDIM+1], di[NDMX+1][MXDIM+1], dt[MXDIM+1],
    dx[MXDIM+1], r[NDMX+1], x[MXDIM+1], xi[MXDIM+1][NDMX+1], xin[NDMX+1];
  static double schi, si, swgt;
  // Best make everything static, allowing restarts

  if(init <= 0)          // Normal entry, enter here on a cold start.
  {                      // Change to mds=0 to disable stratified sampling,
    mds = ndo = 1;       // i.e, use importance sampling only

    idum = 0;            // Initialize random number generator to 0 for every new call of vegas

    for(j = 1; j <= ndim; j++) xi[j][1] = 1.0;
  }
  if(init <= 1) si = swgt = schi = 0.0;

  //Enter here to use the previous grid and make it better.
  if(init <= 2) {

    nd = NDMX;
    ng = 1;

    if(mds) {

      ng = (int)pow(((ncall/2.0)+0.25), (1.0/ndim));
      mds = 1;

      if((2*ng-NDMX) >= 0) {

        mds = -1;
        npg = ng/NDMX + 1;
        nd = ng/npg;
        ng = npg*nd;

      }

    }

    for(k = 1, i = 1; i <= ndim; i++) k *= ng;
    npg = IMAX(ncall/k, 2);
    calls = (double)npg * (double)k;
    dxg = 1.0/ng;
    for(dv2g = 1, i = 1; i <= ndim; i++) dv2g *= dxg;
    dv2g = SQR((calls*dv2g))/npg/npg/(npg-1.0);
    xnd = nd;
    dxg *= xnd;
    xjac = 1.0/calls;

    for(j = 1; j <= ndim; j++) {

      dx[j] = regn[j+ndim] - regn[j];
      xjac *= dx[j];

    }

    if(nd != ndo) {              //Do binning if necessary

      for(i = 1; i <= IMAX(nd, ndo); i++) r[i] = 1.0;
      for(j = 1; j <= ndim; j++) rebin(ndo/xnd, nd, r, xin, xi[j]);
      ndo = nd;

    }

    if(nprn >= 0) {
      printf("%s: ndim=%3d ncall=%8.0f\n",
             " Input parameters for vegas", ndim, calls);
      printf("%28s it=%5d itmx=%5d\n", " ", it, itmx);
      printf("%28s nprn=%5d ALPH=%5.2f\n", " ", nprn, ALPH);
      printf("%28s mds=%3d nd=%4d\n", " ", mds, nd);

      for(j = 1; j <= ndim; j++) {
        printf("%30s xl[%2d]= %11.4g xu[%2d]=%11.4g\n",
               " ", j, regn[j], j, regn[j+ndim]);

      }

    }

  }


  for(it = 1; it <= itmx; it++) {
    // Main iteration loop. Can enter here (init>=3) to do an 
    // additional itmx iterations with all other parameters unchanged.

    ti = tsi = 0.0;

    for(j = 1; j <= ndim; j++) {

      kg[j] = 1;

      for(i = 1; i <= nd; i++) d[i][j] = di[1][j] = 0.0;

    }

    for(;;) {

      fb = f2b = 0.0;

      for(k = 1; k <= npg; k++) {

        wgt = xjac;

        for(j = 1; j <= ndim; j++) {

          xn = (kg[j] - ran2(&idum))*dxg + 1.0;
          ia[j] = IMAX(IMIN((int)(xn), NDMX), 1);

          if(ia[j] > 1) {

            xo = xi[j][ia[j]] - xi[j][ia[j]-1];
            rc = xi[j][ia[j]-1] + (xn-ia[j])*xo;

          } else {

            xo = xi[j][ia[j]];
            rc = (xn - ia[j])*xo;

          }

          x[j] = regn[j] + rc*dx[j];
          wgt *= xo*xnd;

        }

        f = wgt*(*fxn)(x, wgt, arh);
        f2 = f*f;
        fb += f;
	f2b += f2;

        for(j = 1; j <= ndim; j++) {

          di[ia[j]][j] += f;
          if(mds >= 0) d[ia[j]][j] += f2;

        }

      }

      f2b = sqrt((f2b*npg));
      f2b = (f2b-fb)*(f2b+fb);
      if(f2b <= 0.0) f2b = TINY;
      ti += fb;
      tsi += f2b;

      if(mds < 0) {

        for(j = 1; j <= ndim; j++) d[ia[j]][j] += f2b;

      }

      for(k = ndim; k >= 1; k--) {

        kg[k] %= ng;
        if(++kg[k] != 1)  break;

      }

      if(k < 1) break;

    }

    // Compute final result for this iteration.
    tsi *= dv2g;
    wgt = 1.0/tsi;
    si += wgt*ti;
    schi += wgt*ti*ti;
    swgt += wgt;
    *tgral = si/swgt;
    *chi2a = (schi - si*(*tgral)) / (it - 0.9999);

    if(*chi2a < 0.0) *chi2a = 0.0;
    *sd = sqrt((1.0/swgt));
    tsi = sqrt(tsi);

    if (nprn >= 0) {

      printf("%s %3d : integral = %14.7g+/-%9.2g\n",
             " iteration no.", it, ti, tsi);
      printf("%s integral = %14.7g+/-%9.2g chi**2/IT n = %9.2g\n",
             " all iterations: ", *tgral, *sd, *chi2a);

      if(nprn) {

        for(j = 1; j <= ndim; j++) {

          printf("DATA FOR axis %2d\n", j);
          printf("%6s%13s%11s%13s%11s%13s\n",
                 "X", "delta i", "X", "delta i", "X", "delta i");

          for(i = 1+nprn/2; i <= nd; i += nprn+2) {

            printf("%8.5f%12.4g%12.5f%12.4g%12.5f%12.4g\n",
                   xi[j][i], di[i][j], xi[j][i+1],
                   di[i+1][j], xi[j][i+2], di[i+2][j]);

          }

        }

      }

    }

    // Refine grid
    for(j = 1; j <= ndim; j++) {

      xo = d[1][j];
      xn = d[2][j];
      d[1][j] = (xo + xn)/2.0;
      dt[j] = d[1][j];

      for(i = 2; i < nd; i++) {

        rc = xo + xn;
        xo = xn;
        xn = d[i+1][j];
        d[i][j] = (rc + xn)/3.0;
        dt[j] += d[i][j];

      }

      d[nd][j] = (xo + xn)/2.0;
      dt[j] += d[nd][j];

    }

    for(j = 1; j <= ndim; j++) {

      rc = 0.0;

      for(i = 1; i <= nd; i++) {

        if(d[i][j] < TINY) d[i][j] = TINY;
        r[i] = pow(((1.0-(d[i][j]/dt[j])) / (log(dt[j])-log(d[i][j]))), ALPH);
        rc += r[i];

      }

      rebin(rc/xnd, nd, r, xin, xi[j]);

    }

  } // end of main iteration loop


} // end of vegas



/********************************************************************/



void rebin(double rc, int nd, double r[], double xin[], double xi[])
{

  // Utility routine used by vegas, to rebin a vector of 
  // densities xi into new bins defined by a vector r.

  int i, k = 0;
  double dr = 0.0, xn = 0.0, xo = 0.0;

  for(i = 1; i < nd; i++) {

    while (rc > dr)

      dr += r[++k];

    if(k > 1) xo = xi[k-1];

    xn = xi[k];
    dr -= rc;
    xin[i] = xn - (xn - xo)*dr/r[k];

  }

  for(i = 1; i < nd; i++) xi[i] = xin[i];
  xi[nd] = 1.0;

} // end of rebin



/********************************************************************/
#define IM1 2147483563 
#define IM2 2147483399 
#define AM (1.0/IM1) 
#define IMM1 (IM1-1) 
#define IA1 40014 
#define IA2 40692 
#define IQ1 53668 
#define IQ2 52774 
#define IR1 12211 
#define IR2 3791 
#define NTAB 32 
#define NDIV (1+IMM1/NTAB) 
#define EPS 1.2e-7 
#define RNMX (1.0-EPS) 


double ran2(long *idum)
{
  /* Random number generator from Numerical Recipes */

  int j;
  long k;
  static long idum2 = 123456789;
  static long iy = 0;
  static long iv[NTAB];
  double temp;

  if(*idum <= 0) {

    if(-(*idum) < 1) *idum = 1;
    else *idum = -(*idum);
    idum2 = (*idum);

    for (j = NTAB+7; j >= 0; j--) {

      k = (*idum)/IQ1;
      *idum = IA1*(*idum-k*IQ1) - k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;

    }

    iy = iv[0];

  }

  k = (*idum)/IQ1;
  *idum = IA1 * (*idum - k*IQ1) - k*IR1;
  if (*idum < 0) *idum += IM1;
  k = idum2/IQ2;
  idum2 = IA2*(idum2 - k*IQ2) - k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j = iy/NDIV;
  iy = iv[j] - idum2;
  iv[j] = *idum;
  if (iy < 1) iy += IMM1;
  if ((temp = AM*iy) > RNMX) return RNMX;
  else return temp;

} // End of ran2.
