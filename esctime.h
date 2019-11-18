#include "vegas_arh.h"


double I1(double par_in[], double dummy, double a)
{

  double t1, t2, ans;

  t1 = 1.0 - SQR((par_in[1] * sin(par_in[2])));
  t2 = sqrt(t1) - (par_in[1] * cos(par_in[2]));
  ans = par_in[1]*par_in[3]*log(SQR(t2)+SQR((a*par_in[3])));

  return(ans);

} // end of I1



double I2(double inpar[], double secjunk, double a)
{

  double t1 = 0.0, t2 = 0.0, ans = 0.0;

  t1 = 1.0 - SQR((inpar[1] * sin(inpar[2]))); 
  t2 = sqrt(t1) - (inpar[1] * cos(inpar[2]));
  ans = inpar[1] * t2 * asin(a / sqrt((SQR(t2) + SQR(a))));

  return(ans);

} // end of I2



double I3(double par_in[], double dummy, double a)
{

  double t1, t2, ans;

  t1 = 1.0 - SQR((par_in[1] * sin(par_in[2]))); 
  t2 = sqrt(t1) - (par_in[1] * cos(par_in[2]));
  ans = par_in[1] * SQR(t2) * log(1.0 + (SQR(a) / SQR(t2)));

  return(ans);

} // end of I3




/********************************************************************/



double esctimecalc(double radius, double height)
{

  double region[7], regn[5];     // regions being integrated over
  const int num_dim = 2, num_dim_f1 = 3;  // number of dimensions of integral
  int init = 0;                      // Number of previous calls to vegas
  int num_calls = 1000;              // Number of calls to integration function
  int num_iter = 50;                 // number of grid refinements
  int print_par = -1;                 // controls printing out by vegas
  double arh, tesc;             // Radius, height, their ratio & esc. time

  // outputs of vegas integration routine
  double integ1, integ1_sigma, integ1_chi2;            // result of integral,
  double integ2, integ2_sigma, integ2_chi2;            // standard deviation of integral,
  double integ3, integ3_sigma, integ3_chi2;            // & chi squared of integration

  regn[1] = 0.0;
  regn[2] = 0.0;
  regn[3] = 1.0;
  regn[4] = 2.0 * pi;

  region[1] = 0.0;
  region[2] = 0.0;
  region[3] = 0.0;
  region[4] = 1.0;
  region[5] = 2.0 * pi;
  region[6] = 1.0;

  arh = height / radius;

  vegas(region, num_dim_f1, I1, init, num_calls, num_iter, print_par, &integ1, &integ1_sigma, &integ1_chi2, arh);
  init = 1;
  num_calls = 60000;
  num_iter = 1;
  vegas(region, num_dim_f1, I1, init, num_calls, num_iter, print_par, &integ1, &integ1_sigma, &integ1_chi2, arh);
  //printf("\n I1 = %e, sigma = %e, chi2 = %e\n", integ1, integ1_sigma, integ1_chi2);

  init = 0;
  num_calls = 1000;
  num_iter = 50;
  vegas(regn, num_dim, I2, init, num_calls, num_iter, print_par, &integ2, &integ2_sigma, &integ2_chi2, arh);
  init = 1;
  num_calls = 60000;
  num_iter = 1;
  vegas(regn, num_dim, I2, init, num_calls, num_iter, print_par, &integ2, &integ2_sigma, &integ2_chi2, arh);
  //printf("\n I2 = %e, sigma = %e, chi2 = %e\n", integ2, integ2_sigma, integ2_chi2);

  init = 0;
  num_calls = 1000;
  num_iter = 50;
  vegas(regn, num_dim, I3, init, num_calls, num_iter, print_par, &integ3, &integ3_sigma, &integ3_chi2, arh);
  init = 1;
  num_calls = 60000;
  num_iter = 1;
  vegas(regn, num_dim, I3, init, num_calls, num_iter, print_par, &integ3, &integ3_sigma, &integ3_chi2, arh);
  //printf("\n I3 = %e, sigma = %e, chi2 = %e\n", integ3, integ3_sigma, integ3_chi2);

  tesc = (height/(4.0*c)) - (height*log(arh)/(2.0*c)) + (height*integ1/(2.0*pi*c)) + 
    (radius*integ2/(pi*c)) - (SQR(radius)*integ3/(2.0*pi*height*c));
  //printf("\n arh = %e, tesc = %e\n", arh, tesc);

  return(tesc);

} // end of esctimecalc.
