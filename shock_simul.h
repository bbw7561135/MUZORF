/*   This code is based on the work of Maddalena Spada et al. from MNRAS, 2001,
      325, 1559. The code calculates the parameters such as the magnetic field 
     and minimum energy of the electrons in the presence of a shock wave, 
     which will be used as input to the code shock_rad.c 
     Only the collision of two consecutive shells has been considered unlike 
     the paper where they take into account the collision of multiple 
     consecutive shells at regular intervals.
     */ 

/* Function to calculate the position of the collision/CD in the AGN frame. */

double colrad(gami, gamo, tv, outerrad, innerrad, outerwidth)
double gami, gamo, tv, outerrad, innerrad, outerwidth;
{
    double R_c, innerbeta, outerbeta, a, b, deltat, ro;

       a = 1. - (1./SQR(gami));
       b = 1. - (1./SQR(gamo));
       innerbeta = sqrt(a);
       outerbeta = sqrt(b);
       ro = outerrad - innerrad;
       deltat = (ro - outerwidth)/(c * (innerbeta - outerbeta));
       R_c = innerrad + (c * innerbeta * deltat); 

       return (R_c); 
}


/* Parameters of the merged shell calculated here. */

double blfm(double mi, double mo, double gami, double gamo)
{
    double gamm, a, b;

        a = (mi * gami) + (mo * gamo);
        b = (mi/gami) + (mo/gamo);
        gamm = sqrt((a/b));

        return(gamm);
}

double energyinm(double mi, double mo, double gami, double gamo)
{
  double gammam, ein, a, b;

  gammam = blfm(mi, mo, gami, gamo);
  a = gami - gammam;
  b = gamo - gammam;
  ein = ((mi * a) + (mo * b)) * 9.0e20;   // E = gamma*m*c^2
  printf("\n BLF and energy of merged shell = %e, %e\n", gammam, ein);

  return(ein);
}


/* Shock parameters calculated here. */

double widthfs(double gamcap, double gamfs, double outerwidth, double gamo)
{
     double delfs, a, roratio;
 
     //a = gamcap * gamfs;
     //roratio = (gamcap - 1.)/(a + 1.);
     //delfs = outerwidth * roratio;
     //delfs = gamo * outerwidth * roratio;
     delfs = gamo * outerwidth / ((4.0 * gamfs) + 3.0);
     
     return (delfs);
}

double widthrs(double gamcap, double gamrs, double innerwidth, double gami)
{
     double delrs, a, roratio;

     //a = gamcap * gamrs;
     //roratio = (gamcap - 1.)/(a + 1.);
     //delrs = innerwidth * roratio;
     //delrs = gami * innerwidth * roratio;
     delrs = gami * innerwidth / ((4.0 * gamrs) + 3.0);

        return (delrs);
}

void energies(double rhobar_o, double rhobar_i, double gamcap, double gamfs, double gamrs, double *ufs, double *urs)  
//energies(mi, mo, gami, gamo, gamfs, gamrs, gamcap, deli, delo, r, efs, ers)
//double mi, mo, gami, gamo, gamfs, gamrs, gamcap, deli, delo, *efs, *ers, *r;
{
     double delshi, delsho, Ein;
        
     //delshi = widthrs(gamcap, gamrs, deli);
     //delsho = widthfs(gamcap, gamfs, delo); 
     //printf("inner/outer shocked shell widths = %e, %e\n", delshi, delsho);
     //*r = delshi / delsho;
     //printf("Shocked shell widths' ratio = %e\n", *r);
     //Ein = energyinm(mi, mo, gami, gamo);
     //*efs = Ein/(1. + *r);
     //*ers = (Ein * *r)/(1. + *r);
     //*ufs = rhobar_o * SQR(c) * ((gamcap * gamfs) + 1.) * 
     //(gamfs - 1.) / (gamcap - 1.);
     //*urs = rhobar_i * SQR(c) * ((gamcap * gamrs) + 1.) * 
     //(gamrs - 1.) / (gamcap - 1.);
     *ufs = rhobar_o * SQR(c) * ((4.0 * gamfs) + 3.0) * (gamfs - 1.);
     *urs = rhobar_i * SQR(c) * ((4.0 * gamrs) + 3.0) * (gamrs - 1.);
     return;
}
