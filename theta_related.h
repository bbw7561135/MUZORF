#define RGCOEFF 1.48e5 // Value of G*M_solar/c^2
#define THETACOEFF 3.58816e2 // Value of 1.44*sqrt(1.48e5)*(0.176)^(1/4)
#define JCOEFF 157.9137 // Value of 16 * pi^2

void disktemp(double BHM, double disclumin46, double acceff, double *raddist, double *thetar)
{
  double Rg, Rindisk, Routdisk, logstepR, thetafactor;
  double rgfactor, rterm;
  int j;

  Rg = RGCOEFF * BHM;
  Rindisk = 6.002 * Rg;
  Routdisk = 6.0e3 * Rg;
  logstepR = pow((Routdisk / Rindisk), 1.0/((double) (GS - 1)));
  thetafactor = THETACOEFF * pow((disclumin46 * Rg / acceff), 0.25);
  rgfactor = 6.0 * Rg;

  for (j = 0; j < GS; j++) {

    raddist[j] = Rindisk;
    rterm = 1.0 - sqrt((rgfactor / raddist[j]));
    thetar[j] = thetafactor * pow(rterm, 0.25) / pow(raddist[j], 0.75);

    Rindisk *= logstepR;

  }

  return;

}


void blremissiv(double LBLR, double optdep, double covfac, double rinblr, double routblr, double *blrrad, double *jline, double *jcont, double *lineinteg, double *continteg)
{
  double lineterm, contterm, logstepr, lineyofr[GS], contyofr[GS];
  double linep1, linep2, contp1, q, p = 1.5, s = 1.0, linesum = 0.0;
  double contsum = 0.0;
  int j;

  logstepr = pow((routblr / rinblr), 1.0/((double) (GS - 1)));

  q = 1.0 / 3.0;
  linep1 = (2.0 * q) - p - 2.0;
  linep2 = (2.0 * q) - p;
  contp1 = s + 2.0;

  if ((q == 1.0/3.0) && (p == 1.5) && (s == 1.0)) {
    linep1 = -2.83333;
    linep2 = -0.83333;
    contp1 = 3.0;

    linesum = 6.0 * (pow(routblr, (1.0/6.0)) - pow(rinblr, (1.0/6.0)));
    contsum = log(routblr / rinblr);

    for (j = 0; j < GS; j++) {
      blrrad[j] = rinblr;
      rinblr *= logstepr;
    }

  } else {
    lineyofr[0] = pow(rinblr, linep2);

    if (s == 1.0) {
      contsum = log(routblr / rinblr);
    }

    if (s != 1.0) {
      contyofr[0] = 1.0 / pow(rinblr, s);
    }

    blrrad[0] = rinblr;

    rinblr *= logstepr;

    for (j = 1; j < GS; j++) {

      blrrad[j] = rinblr;

      lineyofr[j] = pow(blrrad[j], linep2);
      linesum += (0.5 * (lineyofr[j] + lineyofr[j - 1]) *
		  (blrrad[j] - blrrad[j - 1]));

      if (s != 1.0) {
	contyofr[j] = 1.0 / pow(blrrad[j], s);
	contsum += (0.5 * (contyofr[j] + contyofr[j - 1]) *
		    (blrrad[j] - blrrad[j - 1]));
      }

      rinblr *= logstepr;
    }
  }

  *lineinteg = linesum;
  *continteg = contsum;

  for (j = 0; j < GS; j++) {

    lineterm = LBLR * pow(blrrad[j], linep1);
    jline[j] = lineterm / (JCOEFF * linesum);

    contterm = LBLR * optdep / (pow(blrrad[j], contp1) * covfac);
    jcont[j] = contterm / (JCOEFF * contsum);

  }

  return;

}
