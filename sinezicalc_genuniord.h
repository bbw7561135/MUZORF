/*

To calculate the angle between the magnetic field and the aberrated
line of sight in the plasma frame depending on the various
orientations of the magnetic field in the comoving frame of the
jet. The resultant quantity would be used in calculating the modified
synchrotron and subsequently the SSC emission.

*/

#define STEPSIZE 45

 
void parallelgeom(double dopfac, double thetaobs, double blfsh, double betash, double *sineangle, double *cosineangle)
{

  *cosineangle = dopfac * blfsh * (cos(thetaobs) - betash); 
  *sineangle = dopfac * sin(thetaobs);

  //printf("In function, cosineangle = %e, sineangle = %e\n", *cosineangle, *sineangle);

  return;

}


void transgeom(double dopfac, double thetaobs, double theta_xy, double *sineangle, double *cosineangle)
{

  *cosineangle = dopfac * sin(thetaobs) * cos(theta_xy);
  *sineangle = sqrt(1.0 - SQR(*cosineangle));

  //printf("In function, cosineangle = %e, sineangle = %e\n", *cosineangle, *sineangle);

  return;

}


void obliquegeom(double dopfac, double thetaobs, double theta_xy, double theta_z, double blfsh, double betash, double *sineangle, double *cosineangle)
{

  double term1, term2;

  term1 = sin(thetaobs) * sin(theta_z) * cos(theta_xy);
  term2 = blfsh * cos(theta_z) * (cos(thetaobs) - betash);
  
  *cosineangle = dopfac * (term1 + term2);
  *sineangle = sqrt(1.0 - SQR(*cosineangle));

  //printf("In function, cosineangle = %e, sineangle = %e\n", *cosineangle, *sineangle);

  return;

}


void toroidgeom(double dopfac, double thetaobs, double *toroidsine, double *cosineangle)
{

  double coeff = -1.0, angle;
  int i, phiangle = 0;

  for (i = 0; i < QUANTIZE; i++) {

    angle = (pi * phiangle) / 180.0;
    cosineangle[i] = coeff * dopfac * sin(angle) * sin(thetaobs);
    toroidsine[i] = sqrt(1.0 - SQR(cosineangle[i]));

    phiangle += STEPSIZE;
    /*
    printf("In function, cosineangle[%d] = %e, toroidsine[%d] = %e, phiangle = %d, angle = %e\n",
	   i, cosineangle[i], i, toroidsine[i], phiangle, angle);
    */
  }

  return;

}


void helixgeom(double dopfac, double blfsh, double betash, double thetaobs, double theta_z, double *helixsine, double *cosineangle)
{

  double angle, term1, term2;
  int i, phiangle = 0;

  for (i = 0; i < QUANTIZE; i++) {

    angle = (pi * phiangle) / 180.0;
    term1 = blfsh * cos(theta_z) * (cos(thetaobs) - betash);
    term2 = sin(thetaobs) * sin(theta_z) * sin(angle);    
    cosineangle[i] = dopfac * (term1 - term2);
    helixsine[i] = sqrt(1.0 - SQR(cosineangle[i]));

    phiangle += STEPSIZE;
    /*
    printf("In function, cosineangle[%d] = %e, helixsine[%d] = %e, phiangle = %d, angle = %e\n",
           i, cosineangle[i], i, helixsine[i], phiangle, angle);
    */
  }

  return;

}
