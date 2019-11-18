/*
 Calculation of shock propagation and the subsequent radiation
 processes using synchrotron power with a new empirical expression for
 a power law electron distribution with independent synchrotron
 frequency array. The blob changes its position in time and the
 subsequent evolution of the electron/photon population has been taken
 into account. The evolution of the populations has been calculated
 according to the evolution of radiative mechanisms, such as
 synchrotron, SSC, and EC. The EC mechanism includes the contribution
 of 3 external seed photon fields, namely the accretion disk, BLR, and
 DT, entering the jet in an anisotropic manner. All the quantities
 used and calculated are in the comoving & observer's frame.

 The model also includes the impact of the orientation of the magnetic
 field on the resulting spectral energy distribution and lightcurves
 of a blazar.

 Use of command line arguments, so while compiling this program, use the
 following command:
 output file name file#
 e.g. syn_power_fitting_newsph 1

 Output: nuFnu{nu} [Jy Hz] vs. Photon frequency [Hz]
 */

//#define BLRDT_UPH

#include <stdio.h>
#include <omp.h>

#ifdef __APPLE__
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <time.h>
#include <signal.h>

//#define NUM_THREADS 16

#define H0 70.0
#define VARNUM 44                               // Declares the number of parameters being used in the eve file
#define files 0
#define ZON 100
#define TRUE 1
#define FALSE 0
#define Energy 8.198451e-07           // Value of m_e*c^2
#define MAXPOINT 24
#define EPDMIN 3.23280e-06
#define EPDMAX 1.95422e-02
#define CONTMIN 2.42460e-08
#define CONTMAX 6.97928e-03
#define EPDTMIN 1.6164e-08
#define EPDTMAX 2.04474e-06
#define DTTEMP  1.2e3        // In degree Kelvin
#define TAULIMIT 1.0e-3
#define MAXIT 100
#define c 3.0e10
#define pi 3.141592654
#define SIGMASB 5.6704e-05        // Stefan-Boltzmann Constant in gs^(-3)K^(-4)
#define SIGMAT 6.65e-25	          // Value of Thomson scattering cross section
#define EMASS 9.1093897e-28	  // CGS Mass of an electron
#define PLANCKCONST 6.6260755e-27 // CGS Value of Planck's constant 
#define ECHARGE 4.803206e-10	  // CGS Value of charge

#define FREQ_GRID 150                          // Grid points for photon related arrays
#define ASZ 150             // Grid points for photon/e- related arrays
#define EGRID 150           // Grid points for e- related arrays
#define GS 30              // Grid points for distance related arrays
#define RGRID 30          // Grid points for the disc radius
#define MUGRID 60
#define PHIGRID 30
#define DISKGRID 30
#define Exit_success 0
#define ED 100          // Grid points for emission distance array
#define QUANTIZE  8     // Grid points for the phi array of a toroidal B-field 

#define GAMMAXLIM 3.0e2
#define GAMMAXDTLIM 3.0e2
#define NUFNUMAXLIM 1.0e10

#define SQR(x) (x*x)
void handleInterrupt(int);
int exitSimulation (clock_t, char*, int, char*, char*);


#include "util.h"
#include "tridag.h"                  // Borrowed from Numerical Recipes in C
#include "myrtnew.h"
#include "shock_simul.h"             // Self-made-library programs in C
#include "syn_related_constparams.h"
#include "syn_related_uniform.h"
#include "ssc_related_constparams.h"
#include "theta_related.h"
#include "pair_prod_constparams.h"
#include "ecd_related.h"
#include "ecblr_related.h"
#include "ecdt_related.h"
#include "esctime.h"
#include "telescope/telescope.h"
#include "sinezicalc_genuniord.h"

// the hash map for storing freq and flux that is referenced by observer time
//
shockrad *map = NULL;

//
// GLOBALS
//

double *syn_nu, D_Z = 0.0;
int Ntot = 0, f = 1, r = 1;
time_t start;
char* logfile;
char* sOutputDirectory;
char* nfnfile;
char startbuffer[25];

//
// END GLOBALS
//


double epsummain(double BLF, double betaphterm, double *phengy) {

  double sum = 0.0, expterm, denomterm, blfterm, yofep[CONTGRID];
  int j;

  if (betaphterm == 0.0) {
    blfterm = 1.0;
  } else {
    blfterm = BLF * betaphterm;
  }

  expterm = phengy[0] / AVTEMP;
  denomterm = exp(expterm) - 1.0;

  if (expterm >= 150.0) {
    yofep[0] = 0.0;

  } else if ((expterm >= 24.0) && (expterm < 150.0)) {
    yofep[0] = pow(phengy[0], 3.0) / exp(expterm);

  } else if ((expterm < 24.0) && (expterm > 1.0e-4)) {
    yofep[0] = pow(phengy[0], 3.0) / denomterm;

  } else {
    yofep[0] = pow(phengy[0], 3.0) / expterm;

  }

  if (yofep[0] < 1.0e-100)
    yofep[0] = 0.0;

  for (j = 1; j < CONTGRID; j++) {

    expterm = phengy[j] / AVTEMP;
    denomterm = exp(expterm) - 1.0;

    if (expterm >= 150.0) {
      yofep[j] = 0.0;

    } else if ((expterm >= 24.0) && (expterm < 150.0)) {
      yofep[j] = pow(phengy[j], 3.0) / exp(expterm);

    } else if ((expterm < 24.0) && (expterm > 1.0e-4)) {
      yofep[j] = pow(phengy[j], 3.0) / denomterm;

    } else {
      yofep[j] = pow(phengy[j], 3.0) / expterm;

    }

    if (yofep[j] < 1.0e-100) {
      yofep[j] = 0.0;
      continue;
    }

    sum += (0.5 * (yofep[j] + yofep[j - 1]) * (phengy[j] - phengy[j - 1]));

  }

  return (sum);

}


int main(int argc, char *argv[]) {

  char *sRunId = argv[1];
  sOutputDirectory = argv[2];

  time_t timer;
  struct tm* tm_info;

  int restrict_wall_time = 0;                
// If set to 0 the program will stop when it completes. If it is
// set to 1 it will stop when run time reaches MAX_WALL_TIME - FILE_WRITE_TIME

  double max_wall_time = 0;     // max wall time the program can run in seconds
  double file_write_time = 0;  
// the time, in seconds, to set aside for the files to be generated from the 
// data structure

	// Vars for openmp environment information
	//
	int nthreads, tid, procs, maxt, inpar, dynamic, nested;

        FILE *fp_eve, *fp_log, *fp_freq, *fp_nfn, *fp_theta, *fp_blrj;
	FILE *fp_2401, *fp_eta, *fp_blrin, *fp_nfnblr, *fp_nfnecd;
	FILE *fp_blrIz, *fp_blrIz2, *fp_blrIz3, *fp_dtepnph, *fp_fsnew;
	FILE *fp_dtupnph, *fp_blrupnph;

	FILE *fp_fsynnfn, *fp_fstgg, *fp_fssa, *fp_fseloss, *fp_fsinj;
	FILE *fp_fsscnfn, *fp_fsinit, *fp_fsecdnfn, *fp_fsblrnfn, *fp_fsdtnfn;

	FILE *fp_rsynnfn, *fp_rsinit, *fp_rseloss, *fp_rsinj, *fp_rsscnfn;
	FILE *fp_rsecdnfn, *fp_rsblrnfn, *fp_rsdtnfn, *fp_rsnew;

	time_t walltime, end;
	double elapsed;

	double steps = 0.0, ge = 1.01, gemax = 1.0e7, nui = 0.0, stepsnu = 0.0;
	double D = 0.0, dl = 0.0, dt = 0.0, te_esc = 0.0, mu = 0.0;
	double simtime = 0.0, Beta_sh = 0.0, M_i = 0.0, beta_o = 0.0;
	double contsteps = 0.0, V_i = 0.0, V_o = 0.0, Y = 0.0, xacc = 0.0;
	double x0 = 0.0, Gamma_sh = 0.0, nuFnu_max = 1.0, beta_i = 0.0;
	double z_c = 0.0, freq = 0.0, mu_prime = 0.0, mu_prime_abs = 0.0;
	double dt_param = 0.0, steps_param = 0.0, contepmin = 0.0;
	double dtstepsize = 0.0, dt_epmin = 0.0, dt_area = 0.0, Ldt = 0.0;
	double xterm_main, rzterm_main;
	double emitdist_min, logstepz, stepetaph, etaph_linemitmin;
        double uphblrster_contgausp = 0.0, uphblrster_contgausm = 0.0;
        double uphblrster_linegausp = 0.0, uphblrster_linegausm = 0.0;
	double sinezi = 0.0, thetaobs_rad = 0.0, thetaxy_rad = 0.0;
	double thetaz_rad = 0.0, etassczival = 0.0, sinezifeed = 0.0;

	double tph_esc_fs = 0.0, dtcool_fs = 0.0, Volcyl_fs = 0.0;
	double lmean_fs = 0.0, del_fs = 0.0, znht_fs = 0.0, volfs_param = 0.0;
	double timeobs_up_fs = 0.0, timeobs_si_fs = 0.0, gamma_break_fs = 0.0;
	double dtinj_fs = 0.0, Q0_fs = 0.0, dt_j_fs = 1.0, tcr_fs_zn = 0.0;
	double B_fs = 0.0, time_fs = 1.0, Voldl_fs = 0.0, gmax_upper_fs = 0.0;
	double nu_gamma1_fs = 0.0, nu_B_fs = 0.0, Gamma_fs = 0.0, dt_fs = 0.0;
	double gamma_b_fs = 0.0, beta_fs = 0.0, tcr_fs_com = 0.0, num_fs = 0.0;
 	double gamma_max_fs = 0.0, coeff_fs = 0.0, coeffq_fs = 0.0;
	double amax_fs = 0.0, distph_fs = 0.0, bmax_fs = 0.0, dmax_fs = 0.0;
	double rL_fs = 0.0, voldlfs_param = 0.0, P_up_fs = 0.0;
	double P_side_fs = 0.0, lmeanfs_param = 0.0, e_numden_fs = 0.0;
	double pupfs_param = 0.0, psidefs_param = 0.0;
	double gammadotsynmax_fs = 0.0, totgammadotmin_fs = 0.0;
	double totgammadotmin_zn_fs = 0.0, dtcool_zn_fs = 0.0;
	double dsh_fs = 0.0, Volcylsh_fs = 0.0, volshfs_param = 0.0;
	double Voldlsh_fs = 0.0, voldlshfs_param = 0.0, tph_esc_shfs = 1.0e20;
	double lmean_shfs = 0.0, lmeanshfs_param = 0.0, P_up_shfs = 0.0;
	double pupshfs_param = 0.0, P_side_shfs = 0.0, psideshfs_param = 0.0;
	double dshfs = 0.0, Voldlshfs = 0.0, voldlshfsparam = 0.0;
	double tdelsifs = 0.0, nuFnuecd_fsmax = 1.0, nuFnublr_fsmax = 1.0;
	double nuFnudt_fsmax = 1.0;
	double nu_gamma1_fsrand = 0.0, nu_gamma1_rsrand = 0.0;
	double nu_gamma1_fsfeed = 0.0, nu_gamma1_rsfeed = 0.0;
	//double nu_Bdisord_rs = 0.0, nu_Bdisord_fs = 0.0;
	//double Bord_rs = 0.0, Bord_fs = 0.0, Bdisord_rs = 0.0, Bdisord_fs = 0.0;

	double tph_esc_rs = 0.0, dtcool_rs = 0.0, Volcyl_rs = 0.0;
	double lmean_rs = 0.0, del_rs = 0.0, znht_rs = 0.0, volrs_param = 0.0;
	double timeobs_si_rs = 0.0, e_numden_rs = 0.0, dtinj_rs = 0.0;
	double Q0_rs = 0.0, dt_j_rs = 1.0, tcr_rs_zn = 0.0, B_rs = 0.0;
	double time_rs = 1.0, gmax_upper_rs = 0.0, Voldl_rs = 0.0;
	double distph_rs = 0.0, rL_rs = 0.0, nu_gamma1_rs = 0.0, dmax_rs = 0.0;
	double Gamma_rs = 0.0, gamma_b_rs = 0.0, beta_rs = 0.0, bmax_rs = 0.0;
	double tcr_rs_com = 0.0, gamma_max_rs = 0.0, num_rs = 0.0;
	double coeff_rs = 0.0, coeffq_rs = 0.0, amax_rs = 0.0, dt_rs = 0.0;
	double P_up_rs = 0.0, P_side_rs = 0.0, timeobs_do_rs = 0.0;
	double puprs_param, voldlrs_param = 0.0, lmeanrs_param = 0.0;
	double psiders_param = 0.0, gammadotsynmax_rs = 0.0;
	double gamma_break_rs = 0.0, nu_B_rs = 0.0;
	double totgammadotmin_zn_rs = 0.0, dtcool_zn_rs = 0.0;
	double dsh_rs = 0.0, Volcylsh_rs = 0.0, volshrs_param = 0.0;
	double Voldlsh_rs = 0.0, voldlshrs_param = 0.0, tph_esc_shrs = 1.0e20;
	double lmean_shrs = 0.0, lmeanshrs_param = 0.0, P_up_shrs = 0.0;
	double pupshrs_param = 0.0, P_side_shrs = 0.0, psideshrs_param = 0.0;
	double dshrs = 0.0, Voldlshrs = 0.0, voldlshrsparam = 0.0;
	double tdelsirs = 0.0, nuFnuecd_rsmax = 1.0, nuFnublr_rsmax = 1.0;
        double nuFnudt_rsmax = 1.0;

	double robar_i, robar_o, U_fs, U_rs, eta_IS, Gamma_m;
	double Eint_m, gmt1, gmt2, eint1, eint2, etat1, lineinteg, continteg;
	double lineI, contI, expterm_main, contd1_main, integ_sample = 0.0;
	double diskepmin, eplabinteg_main, stepepdisk, bphtermp, bphtermm;
	double bphterm_lin, temp1, temp2, phimin, phistep, mincos, maxcos;
	double bphterm_main;

	char c1;

	char* evefile = (char*) malloc(sizeof(char) * 1024);
	logfile = (char*) malloc(sizeof(char) * 1024);
	char* freqfile = (char*) malloc(sizeof(char) * 1024);
	nfnfile = (char*) malloc(sizeof(char) * 1024);
	char* thetafile = (char*) malloc(sizeof(char) * 1024);
	char* ecdthlossfile = (char*) malloc(sizeof(char) * 1024);
	char* ecdlossfile = (char*) malloc(sizeof(char) * 1024);
	char* blrjfile = (char*) malloc(sizeof(char) * 1024);
	char* blrIzfile = (char*) malloc(sizeof(char) * 1024);
	char* blrIzfile2 = (char*) malloc(sizeof(char) * 1024);
	char* blrIzfile3 = (char*) malloc(sizeof(char) * 1024);
	char* blrIznphfile = (char*) malloc(sizeof(char) * 1024);
	char* blrnphfile_lin = (char*) malloc(sizeof(char) * 1024);
	char* blrnphfile_gausp = (char*) malloc(sizeof(char) * 1024);
	char* blrnphfile_gausm = (char*) malloc(sizeof(char) * 1024);
	char* blrgamdotfile = (char*) malloc(sizeof(char) * 1024);
	char* ndotblr_file = (char*) malloc(sizeof(char) * 1024);
	char* ndotecd_file = (char*) malloc(sizeof(char) * 1024);
	char* ndotecdth_file = (char*) malloc(sizeof(char) * 1024);
	char* nfnecd_file = (char*) malloc(sizeof(char) * 1024);
	char* nfnblr_file = (char*) malloc(sizeof(char) * 1024);
	char* blrupnph_file = (char*) malloc(sizeof(char) * 1024);
	char* dtepnph_file = (char*) malloc(sizeof(char) * 1024);
	char* dtupnph_file = (char*) malloc(sizeof(char) * 1024);

	char* filename_nfnecd = (char*) malloc(sizeof(char) * 1024);
	char* filename_blrIz = (char*) malloc(sizeof(char) * 1024);
	char* filename_nfnblr = (char*) malloc(sizeof(char) * 1024);
	char* filename_blrupnph = (char*) malloc(sizeof(char) * 1024);
	char* filename_dtepnph = (char*) malloc(sizeof(char) * 1024);
	char* filename_dtupnph = (char*) malloc(sizeof(char) * 1024);

	char* fsinitfile = (char*) malloc(sizeof(char) * 1024);
	char* fssafile = (char*) malloc(sizeof(char) * 1024);
	char* fsynnfnfile = (char*) malloc(sizeof(char) * 1024);
	char* fselossfile = (char*) malloc(sizeof(char) * 1024);
	char* fsnewedenfile = (char*) malloc(sizeof(char) * 1024);
	char* fsinjfile = (char*) malloc(sizeof(char) * 1024);
	char* fsphotfile = (char*) malloc(sizeof(char) * 1024);
	char* fsscnfnfile = (char*) malloc(sizeof(char) * 1024);
	char* fstauggfile = (char*) malloc(sizeof(char) * 1024);
	char* fsecdnfnfile = (char*) malloc(sizeof(char) * 1024);
	char* fsblrnfnfile = (char*) malloc(sizeof(char) * 1024);
	char* fsdtnfnfile = (char*) malloc(sizeof(char) * 1024);
	char* fscomphotfile = (char*) malloc(sizeof(char) * 1024);

	char* rsinitfile = (char*) malloc(sizeof(char) * 1024);
	char* rsynnfnfile = (char*) malloc(sizeof(char) * 1024);
	char* rselossfile = (char*) malloc(sizeof(char) * 1024);
	char* rsnewedenfile = (char*) malloc(sizeof(char) * 1024);
	char* rsinjfile = (char*) malloc(sizeof(char) * 1024);
	char* rsphotfile = (char*) malloc(sizeof(char) * 1024);
	char* rsscnfnfile = (char*) malloc(sizeof(char) * 1024);
	char* rsecdnfnfile = (char*) malloc(sizeof(char) * 1024);
	char* rsblrnfnfile = (char*) malloc(sizeof(char) * 1024);
	char* rsdtnfnfile = (char*) malloc(sizeof(char) * 1024);
	char* rscomphotfile = (char*) malloc(sizeof(char) * 1024);

	// To keep track of files generated at each time step.
	char* filename0_fs = (char*) malloc(sizeof(char) * 1024); 	
	char* filename1_fs = (char*) malloc(sizeof(char) * 1024);
	char* filename2_fs = (char*) malloc(sizeof(char) * 1024);
	char* filename3_fs = (char*) malloc(sizeof(char) * 1024);
	char* filename4_fs = (char*) malloc(sizeof(char) * 1024);
	char* filename5_fs = (char*) malloc(sizeof(char) * 1024);
	char* filename9_fs = (char*) malloc(sizeof(char) * 1024);
	char* filename12_fs = (char*) malloc(sizeof(char) * 1024);
	char* filename8 = (char*) malloc(sizeof(char) * 1024);

	char* filename1_rs = (char*) malloc(sizeof(char) * 1024);
	char* filename2_rs = (char*) malloc(sizeof(char) * 1024);
	char* filename3_rs = (char*) malloc(sizeof(char) * 1024);
	char* filename5_rs = (char*) malloc(sizeof(char) * 1024);
	char* filename9_rs = (char*) malloc(sizeof(char) * 1024);
	char* filename12_rs = (char*) malloc(sizeof(char) * 1024);

	// arrays used in the main program are declared here
	//
	double blrrad[GS], jlineemissiv[GS], jcontemissiv[GS], raddist[GS];
	double thetar[GS], linelabfreq[LINEGRID], nphcontlab[CONTGRID];
	double nphlinelab[LINEGRID], linelablam[LINEGRID], absarr[MAXPOINT];
	double weights[MAXPOINT], etaph_lin[ETASTEP], etaphstar_lin[ETASTEP];
	double etaphstar_gausp[MAXPOINT], etaph_gausp[MAXPOINT];
	double etaphstar_gausm[MAXPOINT], etaph_gausm[MAXPOINT];
	double linelabep[LINEGRID], lineNV[LINEGRID], contlabep[CONTGRID];
	double nph_contlin[ETASTEP][CONTGRID], nph_linelin[ETASTEP][LINEGRID];
	double lineepm[MAXPOINT][LINEGRID], contepm[MAXPOINT][CONTGRID];
	double lineepp[MAXPOINT][LINEGRID], lineeplin[ETASTEP][LINEGRID];
	double contepp[MAXPOINT][CONTGRID], nph_contgausp[MAXPOINT][CONTGRID];
	double conteplin[ETASTEP][CONTGRID], nph_linegausp[MAXPOINT][LINEGRID];
	double coskip_main[PHIGRID][MAXPOINT], disklabep[DISKGRID];
	double nph_contgausm[MAXPOINT][CONTGRID], phi_main[PHIGRID];
	double nph_linegausm[MAXPOINT][LINEGRID];
	double coskim_main[PHIGRID][MAXPOINT], coskilin_main[PHIGRID][ETASTEP];
	double lineinten_gauspmain[MAXPOINT], lineinten_gausmmain[MAXPOINT];
	double continten_gauspmain[MAXPOINT], continten_gausmmain[MAXPOINT];
	double lineinten_linmain[ETASTEP], continten_linmain[ETASTEP];
	double etadtarr[ETADT], mudtstar[ETADT], dtlabpe[DTGRID];
	double dtep[ETADT][DTGRID], nph_dt[ETADT][DTGRID];
	double hypt[RGRID], etaphd[RGRID], p_engy[RGRID][DISKGRID];
	double emit_dist[ED], Izero[ED]; 
	double sineord[QUANTIZE], nu_gam1fs_ord[QUANTIZE], nu_gam1rs_ord[QUANTIZE];
	double nu_gminord[ZON][QUANTIZE], etassczipval[QUANTIZE];

#ifdef BLRDT_UPH

	double uph_dtcom[ED], uph_dtcomster[ED][ETADT], uph_dtcomtot[ED];
	double uphblrcom_gauscont[ED], uphblrcom_gausline[ED];
	double uph_blrcomgaus[ED], uph_blrcometa[ED], uph_blrcomtot[ED], inten_dtcom[ETADT]; 
	
	// ARRAY[4000]
	//
	double **uphblrster_eta, **nph_cont, **nph_line, **uphblrster_conteta, **uphblrster_lineeta;
	double **contep, **lineep, *lineinten_main, *continten_main;
	double *etaph_linemit, *etaphstar_linemit;

#endif

	// ARRAY[ZON]
	//
	double *e_numden, *gam_min, *gam_max, *nu_gmin, *gammadot_synmin;
	double *gammadot_synmax, *gammadot_sscmax, *gammadot_blrmax;
	double *gammadot_blrmin, *tdelsi_rs, *gammadot_sscmin, *totgammadotmax;
	double *gammadot_ecdmax, *gammadot_ecdmin, *tdelsi_fs;
	double *gammadot_dtmax, *gammadot_dtmin;
	double *nu_gmin_rand, *nu_gminfeed; 
	//double *nu_gmin_disord; 

	// ARRAY[FREQ_GRID]
	//
	double *epsil, *nuFnu, *epsilscatter, *nuFnu_fs_up;
	double *nuFnu_fs_si, *nuFnu_rs_do, *nuFnu_rs_si;
	double *nuFnusyn_fs, *nuFnuecd_fs, *nuFnuecd_rs, *nuFnusyn_rs;
	double *nuFnublr_fs, *nuFnublr_rs, *nuFnussc_fs, *nuFnussc_rs;
	double *nuFnudt_fs, *nuFnudt_rs;

	// ARRAY[EGRID]
	//
	double *ecdloss, *ecblrloss, *ecdthloss, *gamma_e, *a_i, *Q_inj_fs;
	double *ne_curnt_fs, *n_e0_fs, *n_e_p0_fs, *esynloss_fs, *Q_inj_rs;
	double *ne_curnt_rs, *n_e0_rs, *n_e_p0_rs, *esynloss_rs, *ecblrthloss;
	double *ec_dtthloss, *ec_dtloss, *ecblrthapprox;

	// ARRAY[ZON][EGRID]
	//
	double **eloss, **esscloss, **b_i, **c_i, **u_i, **n_e_p, **n_e, **r_i;
	double **ndot_gg;

	// ARRAY[ZON][ASZ]
	//
	double **power_syn, **tau_ssa;
	double **n_ph, **tau_gg, **ndot_ph, **ndot_ssc;
	double **power_ssc, **nuLnu_ecd;
	double **ndot_syn_ep, **nuLnu_syn, **nuLnu_ssc;
	double **ndot_ph_up, **ndot_ph_do, **ndot_ph_si, **ndot_ecd;
	double **power_ecd, **nuLnu_blr, **ndot_blr, **power_blr;
	double **ndot_dt, **nuLnu_dt, **power_dt;
	double **n_phsyn, **n_phssc, **n_phecd, **n_phblr, **n_phdt;
	double **ndot_phsyn_up, **ndot_phsyn_do, **ndot_phsyn_si;
	double **ndot_phssc_up, **ndot_phssc_do, **ndot_phssc_si;
	double **ndot_phecd_up, **ndot_phecd_do, **ndot_phecd_si;
	double **ndot_phblr_up, **ndot_phblr_do, **ndot_phblr_si;
	double **ndot_phdt_up, **ndot_phdt_do, **ndot_phdt_si;
	double **ndot_ph_dummy, **n_phup, **n_phdo; 
	double **ndot_phup_up, **ndot_phup_do, **ndot_phup_si;
	double **ndot_phdo_up, **ndot_phdo_do, **ndot_phdo_si;	
	double **n_phgg, **ndot_syn_uni, **ndot_phgg;
	double **ndot_phgg_up, **ndot_phgg_do, **ndot_phgg_si;
	double **tau_ssa_uni, **ndot_syn_unifeed, **tau_ssa_unifeed;
	double **n_phsynfeed, **n_phsynep, **n_phggup, **n_phggdo;
	double **ndot_phsynep_up, **ndot_phsynep_do, **ndot_phsynep_si;
	double **ndot_phsynfeed_up, **ndot_phsynfeed_do, **ndot_phsynfeed_si;
	double **ndot_phggup_up, **ndot_phggup_do, **ndot_phggup_si;
	double **ndot_phggdo_up, **ndot_phggdo_do, **ndot_phggdo_si;
	double **ndot_phfeed, **ndot_phfeed_up, **ndot_phfeed_do;
	double **n_phfeed, **n_phfeedup, **n_phfeeddo, **ndot_phfeed_si;
	double **ndot_phfeedup_up, **ndot_phfeedup_do, **ndot_phfeedup_si;
	double **ndot_phfeeddo_up, **ndot_phfeeddo_do, **ndot_phfeeddo_si;
	double **ndot_syn_tot, **tau_ssa_tot, **ndot_phsyntot_up;
	double **ndot_syn_unifeedtot, **tau_ssa_unifeedtot;
	

	// ARRAY[ZON][ASZ][QUANTIZE]
	//
	double ***ndot_syn_ord, n_phsynord[ZON][ASZ][QUANTIZE]; 
	double ***ndot_synord_up, ***ndot_synord_do, ***ndot_synord_si;
	double ***tau_ssa_ord; 

	// variables in the eve file are declared here
	//
	double q, Z, R_z, g0, Theta_obs, nu1, blf_0, zn_rs, alpha, Ldisk, MBH;
	double nu2, eta, L_w, t_w, t_v, Gammai, Gammao, endtime, zn_fs, R_o;
	double R_i, M_o, del_i, del_o, epsilelec, epsilbfield, zetae, gammacap;
	double L_BLR, tau_BLR, Rin_BLR, Rout_BLR, fcov, etacc, L_DT, FCOVDT;
	double Rin_DT, Rout_DT, Ldiskfrac;
	double thetaxy, thetaz, Bflag, b_ord;

	// the array is initialized with the address of the variables in the
	// required order
	//
	double* arr[VARNUM] = { &L_w, &t_w, &t_v, &Gammai, &Gammao, &blf_0, &g0,
				&M_o, &del_i, &del_o, &R_o, &R_i, &gammacap, &epsilelec,
				&epsilbfield, &zetae, &alpha, &nu1, &nu2, &q, &R_z, &zn_rs, &zn_fs,
				&eta, &Theta_obs, &Ldisk, &MBH, &etacc, &L_BLR, &Rin_BLR, &Rout_BLR,
				&tau_BLR, &fcov, &L_DT, &Rin_DT, &Rout_DT, &Ldiskfrac, &FCOVDT, &Z, &endtime,
	                        &thetaxy, &thetaz, &Bflag, &b_ord};


	// Other variables used in the program
	//
	double yminus1, yplus1, ipluso, gminus1, itimeso;
	double coeff4, coeff3, coeff2, coeff1, coeff0;
	double isqplusosq, gtimesy, gsq, bobiterm;
	double osqplusyisq, oplusyi;

	int i = 0, k = 0, j = 0, u = 0, time_count = 0, zn = 0;
	int l = 0, write_count = 0, FLUX = TRUE, Izindex = 0;
	int p = 0;

	int Nfs = 0, timecount_fs = 0, accf = TRUE;
	int timeobs_upround_fs = 0, timeobs_siround_fs = 0;
	int LDISK_fs = 1, LBLR_fs = 1, LDT_fs = 1;
	int THOMSONECDFS = FALSE, THOMSONECDFSZN[ZON];
	int THOMSONBLRFS = FALSE, THOMSONBLRFSZN[ZON];
	int THOMSONDTFS = FALSE, THOMSONDTFSZN[ZON];

	int Nrs = 0, timecount_rs = 0, accr = TRUE;
	int timeobs_siround_rs = 0, timeobs_doround_rs = 0, GAUSASZ;
	int LDISK_rs = 1, LBLR_rs = 1, LDT_rs = 1;
	int THOMSONECDRS = FALSE, THOMSONECDRSZN[ZON];
	int THOMSONBLRRS = FALSE, THOMSONBLRRSZN[ZON];
	int THOMSONDTRS = FALSE, THOMSONDTRSZN[ZON];

	start = time(NULL);

    /* Start parallel region */
    #pragma omp parallel private(nthreads, tid)
      {

      /* Obtain thread number */
      tid = omp_get_thread_num();

      /* Only master thread does this */
      if (tid == 0)
        {
        printf("Thread %d getting environment info...\n", tid);

        /* Get environment information */
        procs = omp_get_num_procs();
        nthreads = omp_get_num_threads();
        maxt = omp_get_max_threads();
        inpar = omp_in_parallel();
        dynamic = omp_get_dynamic();
        nested = omp_get_nested();

        /* Print environment information */
        printf("Number of processors = %d\n", procs);
        printf("Number of threads = %d\n", nthreads);
        printf("Max threads = %d\n", maxt);
        printf("In parallel? = %d\n", inpar);
        printf("Dynamic threads enabled? = %d\n", dynamic);
        printf("Nested parallelism supported? = %d\n", nested);

        }

      }  /* Done */

    //
    // USING MALLOC() FOR LARGE ARRAYS TO AVOID HEAP OVERFLOW
    //
    epsil = (double*) malloc(sizeof(double) * FREQ_GRID);
    nuFnu = (double*) malloc(sizeof(double) * FREQ_GRID);
    syn_nu = (double*) malloc(sizeof(double) * FREQ_GRID);
    epsilscatter = (double*) malloc(sizeof(double) * FREQ_GRID);
    nuFnusyn_fs = (double*) malloc(sizeof(double) * FREQ_GRID);
    nuFnussc_fs = (double*) malloc(sizeof(double) * FREQ_GRID);
    nuFnuecd_fs = (double*) malloc(sizeof(double) * FREQ_GRID);
    nuFnublr_fs = (double*) malloc(sizeof(double) * FREQ_GRID);
    nuFnudt_fs = (double*) malloc(sizeof(double) * FREQ_GRID);
    nuFnu_fs_up = (double*) malloc(sizeof(double) * FREQ_GRID);
    nuFnu_fs_si = (double*) malloc(sizeof(double) * FREQ_GRID);
    nuFnusyn_rs = (double*) malloc(sizeof(double) * FREQ_GRID);
    nuFnussc_rs = (double*) malloc(sizeof(double) * FREQ_GRID);
    nuFnuecd_rs = (double*) malloc(sizeof(double) * FREQ_GRID);
    nuFnublr_rs = (double*) malloc(sizeof(double) * FREQ_GRID);
    nuFnudt_rs = (double*) malloc(sizeof(double) * FREQ_GRID);
    nuFnu_rs_do = (double*) malloc(sizeof(double) * FREQ_GRID);
    nuFnu_rs_si = (double*) malloc(sizeof(double) * FREQ_GRID);

    ecdloss = (double*) malloc(sizeof(double) * EGRID);
    ecblrloss = (double*) malloc(sizeof(double) * EGRID);
    ecdthloss = (double*) malloc(sizeof(double) * EGRID);
    esynloss_fs = (double*) malloc(sizeof(double) * EGRID);
    Q_inj_rs = (double*) malloc(sizeof(double) * EGRID);
    ne_curnt_rs = (double*) malloc(sizeof(double) * EGRID);
    n_e0_rs = (double*) malloc(sizeof(double) * EGRID);
    n_e_p0_rs = (double*) malloc(sizeof(double) * EGRID);
    esynloss_rs = (double*) malloc(sizeof(double) * EGRID);
    ecblrthloss = (double*) malloc(sizeof(double) * EGRID);
    ecblrthapprox = (double*) malloc(sizeof(double) * EGRID);
    ec_dtthloss = (double*) malloc(sizeof(double) * EGRID);
    ec_dtloss = (double*) malloc(sizeof(double) * EGRID);
    gamma_e = (double*) malloc(sizeof(double) * EGRID);
    a_i = (double*) malloc(sizeof(double) * EGRID);
    Q_inj_fs = (double*) malloc(sizeof(double) * EGRID);
    ne_curnt_fs = (double*) malloc(sizeof(double) * EGRID);
    n_e0_fs = (double*) malloc(sizeof(double) * EGRID);
    n_e_p0_fs = (double*) malloc(sizeof(double) * EGRID);

    e_numden = (double*) malloc(sizeof(double) * ZON);
    gam_min = (double*) malloc(sizeof(double) * ZON);
    gam_max = (double*) malloc(sizeof(double) * ZON);
    nu_gmin = (double*) malloc(sizeof(double) * ZON);
    gammadot_synmin = (double*) malloc(sizeof(double) * ZON);
    gammadot_synmax = (double*) malloc(sizeof(double) * ZON);
    gammadot_sscmax = (double*) malloc(sizeof(double) * ZON);
    gammadot_blrmax = (double*) malloc(sizeof(double) * ZON);
    gammadot_blrmin = (double*) malloc(sizeof(double) * ZON);
    tdelsi_rs = (double*) malloc(sizeof(double) * ZON);
    gammadot_sscmin = (double*) malloc(sizeof(double) * ZON);
    totgammadotmax = (double*) malloc(sizeof(double) * ZON);
    gammadot_ecdmax = (double*) malloc(sizeof(double) * ZON);
    gammadot_ecdmin = (double*) malloc(sizeof(double) * ZON);
    tdelsi_fs = (double*) malloc(sizeof(double) * ZON);
    gammadot_dtmax = (double*) malloc(sizeof(double) * ZON);
    gammadot_dtmin = (double*) malloc(sizeof(double) * ZON);
    nu_gmin_rand = (double*) malloc(sizeof(double) * ZON);
    nu_gminfeed = (double*) malloc(sizeof(double) * ZON);

    // ARRAY[ZON][ASZ]/[EGRID]
    // ARRAY[ZON][ASZ][QUANTIZE]
    //
    power_syn = (double**) malloc(sizeof(double) * ZON);
    eloss = (double**) malloc(sizeof(double) * ZON);
    b_i = (double**) malloc(sizeof(double) * ZON);
    c_i = (double**) malloc(sizeof(double) * ZON);
    u_i = (double**) malloc(sizeof(double) * ZON);
    n_e_p = (double**) malloc(sizeof(double) * ZON);
    tau_ssa = (double**) malloc(sizeof(double) * ZON);
    n_e = (double**) malloc(sizeof(double) * ZON);
    r_i = (double**) malloc(sizeof(double) * ZON);
    n_ph = (double**) malloc(sizeof(double) * ZON);
    tau_gg = (double**) malloc(sizeof(double) * ZON);
    ndot_ph = (double**) malloc(sizeof(double) * ZON);
    ndot_ssc = (double**) malloc(sizeof(double) * ZON);
    ndot_gg = (double**) malloc(sizeof(double) * ZON);
    ndot_syn_ep = (double**) malloc(sizeof(double) * ZON);
    esscloss = (double**) malloc(sizeof(double) * ZON);
    nuLnu_syn = (double**) malloc(sizeof(double) * ZON);
    nuLnu_ssc = (double**) malloc(sizeof(double) * ZON);
    power_ssc = (double**) malloc(sizeof(double) * ZON);
    ndot_ph_up = (double**) malloc(sizeof(double) * ZON);
    ndot_ph_do = (double**) malloc(sizeof(double) * ZON);
    ndot_ph_si = (double**) malloc(sizeof(double) * ZON);
    ndot_ecd = (double**) malloc(sizeof(double) * ZON);
    nuLnu_ecd = (double**) malloc(sizeof(double) * ZON);
    power_ecd = (double**) malloc(sizeof(double) * ZON);
    nuLnu_blr = (double**) malloc(sizeof(double) * ZON);
    ndot_blr = (double**) malloc(sizeof(double) * ZON);
    power_blr = (double**) malloc(sizeof(double) * ZON);
    ndot_dt = (double**) malloc(sizeof(double) * ZON);
    nuLnu_dt = (double**) malloc(sizeof(double) * ZON);
    power_dt = (double**) malloc(sizeof(double) * ZON);
    n_phsyn = (double**) malloc(sizeof(double) * ZON); 
    n_phssc = (double**) malloc(sizeof(double) * ZON);
    n_phecd = (double**) malloc(sizeof(double) * ZON);
    n_phblr = (double**) malloc(sizeof(double) * ZON);
    n_phdt = (double**) malloc(sizeof(double) * ZON);
    ndot_phsyn_up = (double**) malloc(sizeof(double) * ZON);
    ndot_phsyn_do = (double**) malloc(sizeof(double) * ZON); 
    ndot_phsyn_si = (double**) malloc(sizeof(double) * ZON);
    ndot_phssc_up = (double**) malloc(sizeof(double) * ZON);
    ndot_phssc_do = (double**) malloc(sizeof(double) * ZON);
    ndot_phssc_si = (double**) malloc(sizeof(double) * ZON);
    ndot_phecd_up = (double**) malloc(sizeof(double) * ZON);
    ndot_phecd_do = (double**) malloc(sizeof(double) * ZON);
    ndot_phecd_si = (double**) malloc(sizeof(double) * ZON);
    ndot_phblr_up = (double**) malloc(sizeof(double) * ZON);
    ndot_phblr_do = (double**) malloc(sizeof(double) * ZON);
    ndot_phblr_si = (double**) malloc(sizeof(double) * ZON);
    ndot_phdt_up = (double**) malloc(sizeof(double) * ZON);
    ndot_phdt_do = (double**) malloc(sizeof(double) * ZON);
    ndot_phdt_si = (double**) malloc(sizeof(double) * ZON);
    ndot_ph_dummy = (double**) malloc(sizeof(double) * ZON);
    n_phup = (double**) malloc(sizeof(double) * ZON);
    n_phdo = (double**) malloc(sizeof(double) * ZON);
    ndot_phup_up = (double**) malloc(sizeof(double) * ZON);
    ndot_phup_do = (double**) malloc(sizeof(double) * ZON);
    ndot_phup_si = (double**) malloc(sizeof(double) * ZON);
    ndot_phdo_up = (double**) malloc(sizeof(double) * ZON);
    ndot_phdo_do = (double**) malloc(sizeof(double) * ZON);
    ndot_phdo_si = (double**) malloc(sizeof(double) * ZON);
    n_phgg = (double**) malloc(sizeof(double) * ZON);
    ndot_syn_uni = (double**) malloc(sizeof(double) * ZON);
    ndot_phgg = (double**) malloc(sizeof(double) * ZON);
    ndot_phgg_up = (double**) malloc(sizeof(double) * ZON);
    ndot_phgg_do = (double**) malloc(sizeof(double) * ZON);
    ndot_phgg_si = (double**) malloc(sizeof(double) * ZON);
    tau_ssa_uni = (double**) malloc(sizeof(double) * ZON);
    ndot_syn_ord = (double***) malloc(sizeof(double) * ZON);
    ndot_synord_up = (double***) malloc(sizeof(double) * ZON); 
    ndot_synord_do = (double***) malloc(sizeof(double) * ZON);
    ndot_synord_si = (double***) malloc(sizeof(double) * ZON);
    tau_ssa_ord = (double***) malloc(sizeof(double) * ZON);
    ndot_syn_unifeed = (double**) malloc(sizeof(double) * ZON);
    tau_ssa_unifeed = (double**) malloc(sizeof(double) * ZON);
    n_phsynep = (double**) malloc(sizeof(double) * ZON);
    n_phsynfeed = (double**) malloc(sizeof(double) * ZON);
    n_phggup = (double**) malloc(sizeof(double) * ZON);
    n_phggdo = (double**) malloc(sizeof(double) * ZON);
    ndot_phsynep_up = (double**) malloc(sizeof(double) * ZON);
    ndot_phsynep_do = (double**) malloc(sizeof(double) * ZON);
    ndot_phsynep_si = (double**) malloc(sizeof(double) * ZON);
    ndot_phsynfeed_up = (double**) malloc(sizeof(double) * ZON);
    ndot_phsynfeed_do = (double**) malloc(sizeof(double) * ZON);
    ndot_phsynfeed_si = (double**) malloc(sizeof(double) * ZON);
    ndot_phggup_up = (double**) malloc(sizeof(double) * ZON);
    ndot_phggup_do = (double**) malloc(sizeof(double) * ZON);
    ndot_phggup_si = (double**) malloc(sizeof(double) * ZON);
    ndot_phggdo_up = (double**) malloc(sizeof(double) * ZON);
    ndot_phggdo_do = (double**) malloc(sizeof(double) * ZON);
    ndot_phggdo_si = (double**) malloc(sizeof(double) * ZON);
    ndot_phfeed = (double**) malloc(sizeof(double) * ZON);
    ndot_phfeed_up = (double**) malloc(sizeof(double) * ZON);
    ndot_phfeed_do = (double**) malloc(sizeof(double) * ZON);
    n_phfeed = (double**) malloc(sizeof(double) * ZON);
    n_phfeedup = (double**) malloc(sizeof(double) * ZON);
    n_phfeeddo = (double**) malloc(sizeof(double) * ZON);
    ndot_phfeed_si = (double**) malloc(sizeof(double) * ZON);
    ndot_phfeedup_up = (double**) malloc(sizeof(double) * ZON); 
    ndot_phfeedup_do = (double**) malloc(sizeof(double) * ZON); 
    ndot_phfeedup_si = (double**) malloc(sizeof(double) * ZON);
    ndot_phfeeddo_up = (double**) malloc(sizeof(double) * ZON); 
    ndot_phfeeddo_do = (double**) malloc(sizeof(double) * ZON); 
    ndot_phfeeddo_si = (double**) malloc(sizeof(double) * ZON);
    ndot_syn_tot = (double**) malloc(sizeof(double) * ZON); 
    tau_ssa_tot = (double**) malloc(sizeof(double) * ZON);
    ndot_syn_unifeedtot = (double**) malloc(sizeof(double) * ZON);
    tau_ssa_unifeedtot = (double**) malloc(sizeof(double) * ZON);

    for (i = 0; i < ZON; i++) {

        power_syn[i] = (double*) malloc(sizeof(double) * ASZ);
        eloss[i] = (double*) malloc(sizeof(double) * EGRID);
        b_i[i] = (double*) malloc(sizeof(double) * EGRID);
        c_i[i] = (double*) malloc(sizeof(double) * EGRID);
        u_i[i] = (double*) malloc(sizeof(double) * EGRID);
        n_e_p[i] = (double*) malloc(sizeof(double) * EGRID);
        tau_ssa[i] = (double*) malloc(sizeof(double) * ASZ);
        n_e[i] = (double*) malloc(sizeof(double) * EGRID);
        r_i[i] = (double*) malloc(sizeof(double) * EGRID);
        n_ph[i] = (double*) malloc(sizeof(double) * ASZ);
        tau_gg[i] = (double*) malloc(sizeof(double) * ASZ);
        ndot_ph[i] = (double*) malloc(sizeof(double) * ASZ);
        ndot_ssc[i] = (double*) malloc(sizeof(double) * ASZ);
        ndot_gg[i] = (double*) malloc(sizeof(double) * EGRID);
        ndot_syn_ep[i] = (double*) malloc(sizeof(double) * ASZ);
        esscloss[i] = (double*) malloc(sizeof(double) * EGRID);
        nuLnu_syn[i] = (double*) malloc(sizeof(double) * ASZ);
        nuLnu_ssc[i] = (double*) malloc(sizeof(double) * ASZ);
        power_ssc[i] = (double*) malloc(sizeof(double) * ASZ);
        ndot_ph_up[i] = (double*) malloc(sizeof(double) * ASZ);
        ndot_ph_do[i] = (double*) malloc(sizeof(double) * ASZ);
        ndot_ph_si[i] = (double*) malloc(sizeof(double) * ASZ);
        ndot_ecd[i] = (double*) malloc(sizeof(double) * ASZ);
        nuLnu_ecd[i] = (double*) malloc(sizeof(double) * ASZ);
        power_ecd[i] = (double*) malloc(sizeof(double) * ASZ);
	nuLnu_blr[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_blr[i] = (double*) malloc(sizeof(double) * ASZ);
	power_blr[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_dt[i] = (double*) malloc(sizeof(double) * ASZ);
	nuLnu_dt[i] = (double*) malloc(sizeof(double) * ASZ);
	power_dt[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_dt[i] = (double*) malloc(sizeof(double) * ASZ);
	nuLnu_dt[i] = (double*) malloc(sizeof(double) * ASZ);
	power_dt[i] = (double*) malloc(sizeof(double) * ASZ);
	n_phsyn[i] = (double*) malloc(sizeof(double) * ASZ); 
	n_phssc[i] = (double*) malloc(sizeof(double) * ASZ);
	n_phecd[i] = (double*) malloc(sizeof(double) * ASZ);
	n_phblr[i] = (double*) malloc(sizeof(double) * ASZ);
	n_phdt[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phsyn_up[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phsyn_do[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phsyn_si[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phssc_up[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phssc_do[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phssc_si[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phecd_up[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phecd_do[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phecd_si[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phblr_up[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phblr_do[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phblr_si[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phdt_up[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phdt_do[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phdt_si[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_ph_dummy[i] = (double*) malloc(sizeof(double) * ASZ);
	n_phup[i] = (double*) malloc(sizeof(double) * ASZ);
	n_phdo[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phup_up[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phup_do[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phup_si[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phdo_up[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phdo_do[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phdo_si[i] = (double*) malloc(sizeof(double) * ASZ);
	n_phgg[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_syn_uni[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phgg[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phgg_up[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phgg_do[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phgg_si[i] = (double*) malloc(sizeof(double) * ASZ);
	tau_ssa_uni[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_syn_ord[i] = (double**) malloc(sizeof(double) * ASZ);
	ndot_synord_up[i] = (double**) malloc(sizeof(double) * ASZ);
	ndot_synord_do[i] = (double**) malloc(sizeof(double) * ASZ);
	ndot_synord_si[i] = (double**) malloc(sizeof(double) * ASZ);
	tau_ssa_ord[i] = (double**) malloc(sizeof(double) * ASZ);
	ndot_syn_unifeed[i] = (double*) malloc(sizeof(double) * ASZ);
	tau_ssa_unifeed[i] = (double*) malloc(sizeof(double) * ASZ);
	n_phsynep[i] = (double*) malloc(sizeof(double) * ASZ);
	n_phsynfeed[i] = (double*) malloc(sizeof(double) * ASZ);
	n_phggup[i] = (double*) malloc(sizeof(double) * ASZ);
	n_phggdo[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phsynep_up[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phsynep_do[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phsynep_si[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phsynfeed_up[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phsynfeed_do[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phsynfeed_si[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phggup_up[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phggup_do[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phggup_si[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phggdo_up[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phggdo_do[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phggdo_si[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phfeed[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phfeed_up[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phfeed_do[i] = (double*) malloc(sizeof(double) * ASZ);
	n_phfeed[i] = (double*) malloc(sizeof(double) * ASZ);
	n_phfeedup[i] = (double*) malloc(sizeof(double) * ASZ);
	n_phfeeddo[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phfeed_si[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phfeedup_up[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phfeedup_do[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phfeedup_si[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phfeeddo_up[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phfeeddo_do[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_phfeeddo_si[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_syn_tot[i] = (double*) malloc(sizeof(double) * ASZ);
	tau_ssa_tot[i] = (double*) malloc(sizeof(double) * ASZ);
	ndot_syn_unifeedtot[i] = (double*) malloc(sizeof(double) * ASZ);
	tau_ssa_unifeedtot[i] = (double*) malloc(sizeof(double) * ASZ);

	for (j = 0; j < ASZ; j++) {

	  ndot_syn_ord[i][j] = (double*) malloc(sizeof(double) * QUANTIZE);
	  ndot_synord_up[i][j] = (double*) malloc(sizeof(double) * QUANTIZE);
	  ndot_synord_do[i][j] = (double*) malloc(sizeof(double) * QUANTIZE);
	  ndot_synord_si[i][j] = (double*) malloc(sizeof(double) * QUANTIZE);
	  tau_ssa_ord[i][j] = (double*) malloc(sizeof(double) * QUANTIZE);

	}

    }


#ifdef BLRDT_UPH

    lineinten_main = (double*) malloc(sizeof(double) * 4000);
    continten_main = (double*) malloc(sizeof(double) * 4000);
    etaph_linemit = (double*) malloc(sizeof(double) * 4000); 
    etaphstar_linemit = (double*) malloc(sizeof(double) * 4000);

    uphblrster_eta = (double**) malloc(sizeof(double) * ED);
    uphblrster_conteta = (double**) malloc(sizeof(double) * ED);
    uphblrster_lineeta = (double**) malloc(sizeof(double) * ED);    
    for (i = 0; i < ED; i++) {

      uphblrster_eta[i] = (double*) malloc(sizeof(double) * 4000);
      uphblrster_conteta[i] = (double*) malloc(sizeof(double) * 4000);
      uphblrster_lineeta[i] = (double*) malloc(sizeof(double) * 4000);

    }

    contep = (double**) malloc(sizeof(double) * 4000);
    lineep = (double**) malloc(sizeof(double) * 4000);
    nph_cont = (double**) malloc(sizeof(double) * 4000); 
    nph_line = (double**) malloc(sizeof(double) * 4000);
    for (i = 0; i < 4000; i++) {

      contep[i] = (double*) malloc(sizeof(double) * CONTGRID);
      lineep[i] = (double*) malloc(sizeof(double) * LINEGRID);
      nph_cont[i] = (double*) malloc(sizeof(double) * CONTGRID);
      nph_line[i] = (double*) malloc(sizeof(double) * LINEGRID);

    }

#endif
    

	// Command Line Arguments used in order to keep track of the simulation
	// number
	//
	sprintf(evefile, "eve_cylfsrs%s.txt", sRunId);
	sprintf(logfile, "log_cyl%s.dat", sRunId);
	sprintf(freqfile, "freq_com_cyl%s.dat", sRunId);
	sprintf(nfnfile, "nuFnu_com_cyl%s_t", sRunId);
	sprintf(thetafile, "disk_theta%s.dat", sRunId);
	sprintf(ecdthlossfile, "disk_Thloss%s_", sRunId);
	sprintf(ecdlossfile, "disk_loss%s_", sRunId);
	sprintf(blrjfile, "blr_j%s.dat", sRunId);
	sprintf(blrIzfile, "blr_Iz%s_zht", sRunId);
	sprintf(blrIzfile2, "blr_Iz2%s_zht", sRunId);
	sprintf(blrIzfile3, "blr_Iz3%s_zht", sRunId);
	sprintf(blrIznphfile, "blr_Iznph%s_zht", sRunId);
	sprintf(blrnphfile_lin, "blr_nphlin%s_zht", sRunId);
	sprintf(blrnphfile_gausp, "blr_nphgausp%s_zht", sRunId);
	sprintf(blrnphfile_gausm, "blr_nphgausm%s_zht", sRunId);
	sprintf(blrgamdotfile, "blr_gamdot%s_zht", sRunId);
	sprintf(ndotblr_file, "ndotblr%s_zht", sRunId);
	sprintf(ndotecd_file, "ndotecd%s_zht", sRunId);
	sprintf(ndotecdth_file, "ndotecdth%s_zht", sRunId);
	sprintf(nfnecd_file, "nuFnucomp_fs%s_zht", sRunId);
	sprintf(nfnblr_file, "nuFnucomp_rs%s_zht", sRunId);
	sprintf(blrupnph_file, "upnphblr%s", sRunId);
	sprintf(dtepnph_file, "epnphdt%s_zht", sRunId);	
	sprintf(dtupnph_file, "upnphdt%s", sRunId);

	sprintf(fsinitfile, "fsinit_cyl%s.dat", sRunId);
	sprintf(fssafile, "fssa_cyl%s_t", sRunId);
	sprintf(fsynnfnfile, "fsyn_cyl%s_t", sRunId);
	sprintf(fselossfile, "fseloss_cyl%s_t", sRunId);
	sprintf(fsnewedenfile, "fsenew_cyl%s_t", sRunId);
	sprintf(fsinjfile, "fsinj_cyl%s.dat", sRunId);
	sprintf(fsphotfile, "fsphot_cyl%s_t", sRunId);
	sprintf(fsscnfnfile, "fssc_cyl%s_t", sRunId);
	sprintf(fsecdnfnfile, "fsecd_cyl%s_t", sRunId);
	sprintf(fsblrnfnfile, "fsblr_cyl%s_t", sRunId);
	sprintf(fsdtnfnfile, "fsdt_cyl%s_t", sRunId);
	sprintf(fstauggfile, "fstau_cyl%s_t", sRunId);
	sprintf(fscomphotfile, "fscomphot_cyl%s_t", sRunId);

	sprintf(rsinitfile, "rsinit_cyl%s.dat", sRunId);
	sprintf(rsynnfnfile, "rsyn_cyl%s_t", sRunId);
	sprintf(rselossfile, "rseloss_cyl%s_t", sRunId);
	sprintf(rsnewedenfile, "rsenew_cyl%s_t", sRunId);
	sprintf(rsinjfile, "rsinj_cyl%s.dat", sRunId);
	sprintf(rsphotfile, "rsphot_cyl%s_t", sRunId);
	sprintf(rsscnfnfile, "rssc_cyl%s_t", sRunId);
	sprintf(rsecdnfnfile, "rsecd_cyl%s_t", sRunId);
	sprintf(rsblrnfnfile, "rsblr_cyl%s_t", sRunId);
	sprintf(rsdtnfnfile, "rsdt_cyl%s_t", sRunId);
	sprintf(rscomphotfile, "rscomphot_cyl%s_t", sRunId);

	if (argc < 3 || argc > 6) {

	  printf ("\n\nIncorrect command line arguments.\n\nSyntax: ./shockrad <RUNID> <Output Directory> [ <restrict-wall-time (0/1)> <max-wall-time> <file-write-time> ]\n\nrestrict-wall-time: 1 to stop program at max-wall-time, 0 to disable\nmax-wall-time: number of seconds to run the simulation\nfile-write-time: seconds to set aside in the end to write data to files before stopping the simulation\n\n");

	  exit(-1);
	}

	if(argc == 6) {
	  
	  restrict_wall_time = atoi(argv[3]);

	  if(restrict_wall_time != 0 && restrict_wall_time != 1) {
	    
	    printf ("\n\nIncorrect command line arguments.\n\nSyntax: ./shockrad <RUNID> <Output Directory> [ <restrict-wall-time (0/1)> <max-wall-time> <file-write-time> ]\n\nrestrict-wall-time: 1 to stop program at max-wall-time, 0 to disable\nmax-wall-time: number of seconds to run the simulation\nfile-write-time: seconds to set aside in the end to write data to files before stopping the simulation\n\n");
	    exit(-1);
	  }

	  if(restrict_wall_time == 1) {
	    
	    max_wall_time = atof(argv[4]);
	    file_write_time = atof(argv[5]);
	  }
	}

        time(&timer);
        tm_info = localtime(&timer);

	strftime(startbuffer, 25, "%Y:%m:%d %H:%M:%S", tm_info);
        printf("\n\n\nSimulation started: %s\n\n\n", startbuffer);


	printf("\nRestrict wall time = %d", restrict_wall_time);
	printf("\nmax wall time = %e", max_wall_time);
	printf("\nfile write time = %e", file_write_time);

	printf("\n\nThe number of entries in the eve file should match the VARNUM number\n");

	// signal handler to catch ctrl-c and write the data into files before exiting the program.
	//
	signal(SIGINT, handleInterrupt);

	// Initializing all the arrays being used in the program to zero.
	//
	for (j = 0; j < EGRID; j++) {
		a_i[j] = ecdloss[j] = 0.0;
	}

	for (zn = 0; zn < ZON; zn++) {
	  for (i = 0; i < FREQ_GRID; i++) {
	    n_ph[zn][i] = ndot_ph[zn][i] = ndot_ssc[zn][i] = ndot_gg[zn][i] = 0.0;
	    ndot_ecd[zn][i] = ndot_blr[zn][i] = ndot_dt[zn][i] = 0.0;
	    tau_gg[zn][i] = ndot_syn_ep[zn][i] = tau_ssa[zn][i] = 0.0;
	    ndot_ph_up[zn][i] = ndot_ph_do[zn][i] = ndot_ph_si[zn][i] = 0.0;
	    n_phsyn[zn][i] = n_phssc[zn][i] = n_phecd[zn][i] = n_phblr[zn][i] = 0.0;
	    n_phdt[zn][i] = ndot_phsyn_up[zn][i] = ndot_phsyn_do[zn][i] = 0.0;
	    ndot_phsyn_si[zn][i] = ndot_phssc_up[zn][i] = ndot_phssc_do[zn][i] = 0.0;
	    ndot_phssc_si[zn][i] = ndot_phecd_up[zn][i] = ndot_phecd_do[zn][i] = 0.0;
	    ndot_phecd_si[zn][i] = ndot_phblr_up[zn][i] = ndot_phblr_do[zn][i] = 0.0;
	    ndot_phblr_si[zn][i] = ndot_phdt_up[zn][i] = ndot_phdt_do[zn][i] = 0.0;
	    ndot_phdt_si[zn][i] = ndot_ph_dummy[zn][i] = n_phup[zn][i] = 0.0;
	    n_phdo[zn][i] = ndot_phup_up[zn][i] = ndot_phup_do[zn][i] = 0.0;
	    ndot_phup_si[zn][i] = ndot_phdo_up[zn][i] = ndot_phdo_do[zn][i] = 0.0;
	    ndot_phdo_si[zn][i] = 0.0;
            n_phgg[zn][i] = ndot_syn_uni[zn][i] = 0.0;
	    ndot_phgg_up[zn][i] = ndot_phgg_do[zn][i] = ndot_phgg_si[zn][i] = 0.0;
	    ndot_phgg[zn][i] = tau_ssa_uni[zn][i] = 0.0;
	    ndot_syn_unifeed[zn][i] = tau_ssa_unifeed[zn][i] = n_phsynep[zn][i] = 0.0;
	    n_phsynfeed[zn][i] = n_phggup[zn][i] = n_phggdo[zn][i] = 0.0;
	    ndot_phsynep_up[zn][i] = ndot_phsynep_do[zn][i] = ndot_phsynep_si[zn][i] = 0.0;
	    ndot_phsynfeed_up[zn][i] = ndot_phsynfeed_do[zn][i] = ndot_phsynfeed_si[zn][i] = 0.0;
	    ndot_phggup_up[zn][i] = ndot_phggup_do[zn][i] = ndot_phggup_si[zn][i] = 0.0;
	    ndot_phggdo_up[zn][i] = ndot_phggdo_do[zn][i] = ndot_phggdo_si[zn][i] = 0.0;
	    ndot_phfeed[zn][i] = ndot_phfeed_up[zn][i] = ndot_phfeed_do[zn][i] = 0.0; 
	    n_phfeed[zn][i] = n_phfeedup[zn][i] = n_phfeeddo[zn][i] = 0.0;
	    ndot_phfeed_si[zn][i] = ndot_phfeedup_up[zn][i] = 0.0;
	    ndot_phfeedup_do[zn][i] = ndot_phfeedup_si[zn][i] = ndot_phfeeddo_up[zn][i] = 0.0;
	    ndot_phfeeddo_do[zn][i] = ndot_phfeeddo_si[zn][i] = 0.0;
	    ndot_syn_tot[zn][i] = tau_ssa_tot[zn][i] = 0.0;
	    ndot_syn_unifeedtot[zn][i] = tau_ssa_unifeedtot[zn][i] = 0.0;
	   
	    for (p = 0; p < QUANTIZE; p++) {
	      
	      ndot_syn_ord[zn][i][p] = n_phsynord[zn][i][p] = 0.0;
	      ndot_synord_up[zn][i][p] = ndot_synord_si[zn][i][p] = 0.0;
	      tau_ssa_ord[zn][i][p]  = 0.0;

	    }

	  }
	  
	  for (j = 0; j < EGRID; j++) {
	    n_e_p[zn][j] = esscloss[zn][j] = 0.0;
	  }

	}

	for (zn = 0; zn < ZON; zn++) {
		gam_min[zn] = gam_max[zn] = nu_gmin[zn] = gammadot_synmin[zn] = 0.0;
		gammadot_synmax[zn] = gammadot_sscmax[zn] = gammadot_sscmin[zn] = 0.0;
		THOMSONECDFSZN[zn] = THOMSONECDRSZN[zn] = THOMSONBLRFSZN[zn] = 0;
		THOMSONBLRRSZN[zn] = THOMSONDTFSZN[zn] = THOMSONDTRSZN[zn] = 0;
		nu_gmin_rand[zn] = nu_gminfeed[zn] = 0.0;
		//nu_gmin_disord[zn] = 0.0;
	}

	for (i = 0; i < GS; i++) {
		raddist[i] = thetar[i] = jlineemissiv[i] = jcontemissiv[i] = 0.0;
		blrrad[i] = 0.0;
	}

	for (i = 0; i < MAXPOINT; i++) {
		for (j = 0; j < FREQ_GRID; j++) {
			nph_contgausp[i][j] = 0.0;
			nph_linegausp[i][j] = nph_contgausm[i][j] = nph_linegausm[i][j] = 0.0;
		}
	}

	for (i = 0; i < ETASTEP; i++) {
		for (j = 0; j < FREQ_GRID; j++) {
			nph_linelin[i][j] = nph_contlin[i][j] = 0.0;
		}
	}

	for (i = 0; i < QUANTIZE; i++) {
	  sineord[i] = 1.0;
	  nu_gam1fs_ord[i] = nu_gam1rs_ord[i] = 0.0;
	}

	for (zn = 0; zn < ZON; zn++) {
	  for (i = 0; i < QUANTIZE; i++) {
	    nu_gminord[zn][i] = 0.0;
	  }
	}


	// the file that has the data
	//
	fp_eve = fopen(evefile, "rb");
	if (!fp_eve) {

		printf("\nCould not open the eve file: %s \n %s\n\n", evefile,
				strerror(errno));
		return -1;
	}

	fp_log = fopen(logfile, "a");
	if (!fp_log) {

		printf("\nCould not open the output file log. \n");
		return -1;
	}

	fp_freq = fopen(freqfile, "w");
	if (!fp_freq) {

		printf("\nCould not open the output file freq. \n");
		return -1;
	}

	fp_fsinj = fopen(fsinjfile, "w");
	if (!fp_fsinj) {

		printf("\nCould not open the output file fsinj. \n");
		return -1;
	}

	fp_rsinj = fopen(rsinjfile, "w");
	if (!fp_rsinj) {

		printf("\nCould not open the output file rsinj. \n");
		return -1;
	}

	// go into a loop and read the values from the file and assign them to
	// the corresponding variables
	//
	for (u = 0; u < VARNUM && !feof(fp_eve); u++) {

		while ( !feof(fp_eve)  ) {

			c1 = fgetc(fp_eve);

			if(c1 == ':') {

                break;
			}
		}

		while ( !feof(fp_eve)  ) {

			c1 = fgetc(fp_eve);

			if(c1 != ' ') {

                break;
			}
		}

		fseek(fp_eve, -1, SEEK_CUR);
		fscanf(fp_eve, "%lf\n", arr[u]);
	}

	fclose(fp_eve);

	printf("\nL_w = %e, t_eject = %e, t_int = %e", L_w, t_w, t_v);
	printf("\ngam_index = %e, M_o = %e", gammacap, M_o);
	printf("\nInner, Outer shell BLF = %e %e", Gammai, Gammao);
	printf("\nInner, Outer shell width = %e %e", del_i, del_o);
	printf("\nInner, Outer shell position = %e %e", R_i, R_o);
	printf("\nepsil_e = %e, epsil_b = %e", epsilelec, epsilbfield);
	printf("\nzeta_e = %e, g0 = %e, alpha = %e", zetae, g0, alpha);
	printf("\nZone radius = %e, blf_0 = %e", R_z, blf_0);
	printf("\nZone number FS = %e, RS = %e", zn_fs, zn_rs);
	printf("\nenergy index = %e, Observing angle = %e", q, Theta_obs);
	printf("\nInitial, Final syn. frequency = %e, %e", nu1, nu2);
	printf("\nLdisk = %e, MBH = %e, eta_acc = %e", Ldisk, MBH, etacc);
	printf("\nLblr = %e, tau_blr = %e, fcov = %e", L_BLR, tau_BLR, fcov);
	printf("\nRin_BLR = %e, Rout_BLR = %e", Rin_BLR, Rout_BLR);
	printf("\nL_DT = %e, DTcov = %e, DTTEMP = %e\n", L_DT, FCOVDT, DTTEMP);
	printf("\nredshift = %e, escape parameter = %e", Z, eta);
	printf("\ntheta_xy = %e, theta_z = %e\n", thetaxy, thetaz); 
	printf("\nBflag = %e, Ordered fraction = %e\n", Bflag, b_ord);
	printf("\nSimulation stop time = %e\n", endtime);

	fprintf(fp_log, "\nL_w = %e, t_eject = %e, t_int = %e", L_w, t_w, t_v);
	fprintf(fp_log, "\ngam_index = %e, M_o = %e", gammacap, M_o);
	fprintf(fp_log, "\nInner, Outer shell BLF = %e %e", Gammai, Gammao);
	fprintf(fp_log, "\nInner, Outer shell width = %e %e", del_i, del_o);
	fprintf(fp_log, "\nInner, Outer shell position = %e %e", R_i, R_o);
	fprintf(fp_log, "\nepsil_e = %e, epsil_b = %e", epsilelec, epsilbfield);
	fprintf(fp_log, "\nzeta_e = %e, g0 = %e, alpha = %e", zetae, g0, alpha);
	fprintf(fp_log, "\nZone radius = %e, blf_0 = %e", R_z, blf_0);
	fprintf(fp_log, "\nZone number FS = %e, RS = %e", zn_fs, zn_rs);
	fprintf(fp_log, "\nenergy index = %e, Observing angle = %e", q, Theta_obs);
	fprintf(fp_log, "\nInitial, Final syn. frequency = %e, %e", nu1, nu2);
	fprintf(fp_log, "\nLdisk = %e, MBH = %e, eta_acc = %e", Ldisk, MBH, etacc);
	fprintf(fp_log, "\nLblr = %e, tau_blr = %e, fcov = %e", L_BLR, tau_BLR, fcov);
	fprintf(fp_log, "\nRin_BLR = %e, Rout_BLR = %e, L_DT = %e, FCOVDT = %e", Rin_BLR,
		Rout_BLR, L_DT, FCOVDT);
	fprintf(fp_log, "\nRin_DT = %e, Rout_DT = %e, Ldiskfrac = %e", Rin_DT,
		Rout_DT, Ldiskfrac);
	fprintf(fp_log, "\nredshift = %e, escape parameter = %e", Z, eta);
	fprintf(fp_log, "\nSimulation stop time = %e", endtime);
	fprintf(fp_log, "\ntheta_xy = %e, theta_z = %e in comoving frame\n", 
		thetaxy, thetaz);
	fprintf(fp_log, "\nField flag = %e, Fraction of ordered component = %e\n",
		Bflag, b_ord);

        if ((Bflag == 0.0) && (b_ord != 0.0)) {

	  printf("Error, when Bflag = 0, b_ord can only be 0. Change its value\n");
	  fclose(fp_log);
	  exitSimulation (start, logfile, -1, sOutputDirectory, nfnfile);
        }

	dl = dL_calc(Z, H0);
	mu = cos(pi * Theta_obs / 180.0);
	te_esc = eta * (R_z / c);

	// Velocity of outer & inner shell calculated here.
	//
	beta_o = beta(Gammao);
	beta_i = beta(Gammai);

	// Mass, volume and density of the inner shell declared here.
	//
	M_i = ((L_w * t_w) - (M_o * Gammao * 9.0e20)) / (Gammai * 9.0e20);
	// Inner shell rest mass
	V_i = pi * SQR(R_z) * del_i; // In the lab frame
	robar_i = M_i / (V_i * Gammai); // In comoving frame of inner shell.
	fprintf(fp_log, "M_i = %e, V_i = %e, robar_i = %e\n", M_i, V_i, robar_i);

	// Volume and density of the outer shell declared here.
	//
	V_o = pi * SQR(R_z) * del_o; // In the lab frame.
	robar_o = M_o / (V_o * Gammao); // In comoving frame of outer shell.
	fprintf(fp_log, "M_o = %e, V_o = %e, robar_o = %e\n", M_o, V_o, robar_o);

	// BLF & internal energy of merged shell and conversion efficiency,
	// in Lab frame
	//
	gmt1 = (M_i * Gammai) + (M_o * Gammao);
	gmt2 = (M_i / Gammai) + (M_o / Gammao);
	Gamma_m = sqrt((gmt1 / gmt2));

	eint1 = Gammai - Gamma_m;
	eint2 = Gammao - Gamma_m;
	Eint_m = ((M_i * eint1) + (M_o * eint2)) * 9.0e20;

	etat1 = ((M_i + M_o) * Gamma_m) / ((M_i * Gammai) + (M_o * Gammao));
	eta_IS = 1.0 - etat1;
	fprintf(fp_log, "Gamma_m = %e, Eint_m = %e, eta_IS = %e\n", Gamma_m, Eint_m,
			eta_IS);

	// Minimum RLF & BLF for the shocked fluid calculated here.
	//
	Y = robar_i / robar_o;
	fprintf(fp_log, "Y = %e\n", Y);

	// Value of accuracy within which BLF of the shocked fluid should be
	// calculated.
	xacc = 1.0e-3;
	x0 = blf_0;

	yminus1 = Y - 1.0;
	isqplusosq = SQR(Gammai) + SQR(Gammao);
	gminus1 = gammacap - 1.0;
	yplus1 = Y + 1.0;
	gtimesy = gammacap * Y;
	gsq = SQR(gammacap);
	itimeso = Gammai * Gammao;
	ipluso = Gammai + Gammao;
	bobiterm = 1.0 - (beta_o * beta_i);
	osqplusyisq = SQR(Gammao) + SQR((Y * Gammai));
	oplusyi = Gammao + (Y * Gammai);

	coeff4 =
	  gsq
	  * (SQR(yminus1)
	     + (4.0 * Y
		* (isqplusosq
		   - (2.0 * SQR(itimeso) * bobiterm))));
	
	coeff3 = 2.0 * gammacap * (1.0 - gammacap)
	  * ((yplus1 * oplusyi) - (2.0 * Y * itimeso * ipluso * bobiterm));
	
	coeff2 = (2.0 * gtimesy * (2.0 - (3.0 * gammacap)) * isqplusosq)
	  + (2.0 * gammacap * (gammacap - 2.0) * osqplusyisq)
	  - (2.0 * Y * itimeso * bobiterm
	     * (SQR(gminus1) - (4.0 * gsq * itimeso)))
	  + (4.0 * gtimesy * gminus1) + ((1.0 - gsq) * (1.0 + SQR(Y)));
	
	coeff1 =
	  2.0 * (1.0 - gammacap)
	  * ((Y * ipluso
	      * (1.0
		 + (gammacap
		    * (itimeso
		       * (1.0
			  - (2.0 * beta_o
			     * beta_i))
		       - 1.0))))
	     - ((1.0 + gammacap) * (Gammao + (SQR(Y) * Gammai)))
	     + (gammacap
		* (pow(Gammao, 3.0)
		   + (SQR(Y) * pow(Gammai, 3.0)))));
	
	coeff0 = (gsq * (pow(Gammao, 4.0) + (SQR(Y) * pow(Gammai, 4.0))))
	  + (2.0 * gtimesy
	     * ((gminus1 * isqplusosq)
		- (gammacap * (SQR(itimeso) + 1.0)) + 2.0))
	  + ((1.0 - gsq) * osqplusyisq)
	  - (2.0 * SQR(gminus1) * Y * itimeso * beta_o * beta_i) - (2.0 * Y);

	//printf("\n coeff4 = %e, coeff3 = %e, coeff2 = %e\n", coeff4, coeff3, coeff2);
	//printf("coeff1 = %e, coeff0 = %e\n\n", coeff1, coeff0);

	Gamma_sh = myrtnewt(funcal, x0, xacc, coeff4, coeff3, coeff2, coeff1,
			coeff0);

	if ((Gammao > Gamma_sh) || (Gamma_sh > Gammai)) {

		printf("Error, change the value of blf_0 as Gamma_sh = %e\n", Gamma_sh);
		fclose(fp_log);
		exitSimulation (start, logfile, -1, sOutputDirectory, nfnfile);
	}

	Beta_sh = beta(Gamma_sh);
	mu_prime = (mu - Beta_sh) / (1.0 - (Beta_sh * mu));
	fprintf(fp_log, "Gamma_sh = %e, Beta_sh = %e, mu = %e, mu_prime = %e\n",
			Gamma_sh, Beta_sh, mu, mu_prime);
	D = 1.0 / (Gamma_sh * (1.0 - (Beta_sh * mu)));
	D_Z = D / (1.0 + Z);
	fprintf(fp_log, "D = %e, D_Z = %e, dl = %e\n", D, D_Z, dl);

	// BLF, velocity and energy of FS/RS calculated here.
	//
	Gamma_fs = (Gammao * Gamma_sh) * (1.0 - (Beta_sh * beta_o));
	beta_fs = beta(Gamma_fs);
	coeff_fs = 1837.0 * (epsilelec / zetae) * (Gamma_fs - 1.0);
	coeffq_fs = ((q - 2.0) / (q - 1.0)) * coeff_fs;

	Gamma_rs = (Gammai * Gamma_sh) * (1.0 - (Beta_sh * beta_i));
	beta_rs = beta(Gamma_rs);
	coeff_rs = 1837.0 * (epsilelec / zetae) * (Gamma_rs - 1.0);
	coeffq_rs = ((q - 2.0) / (q - 1.0)) * coeff_rs;

	fprintf(fp_log, "beta_fs = %e, beta_o = %e\n", beta_fs, beta_o);
	fprintf(fp_log, "beta_rs = %e, beta_i = %e\n", beta_rs, beta_i);
	fprintf(fp_log, "Gamma_fs = %e, Gamma_rs = %e\n", Gamma_fs, Gamma_rs);

	energies(robar_o, robar_i, gammacap, Gamma_fs, Gamma_rs, &U_fs, &U_rs);
	fprintf(fp_log, "U_fs = %e, U_rs = %e\n", U_fs, U_rs); //In comoving frame

	del_fs = widthfs(gammacap, Gamma_fs, del_o, Gammao); //In comoving frame.
	del_rs = widthrs(gammacap, Gamma_rs, del_i, Gammai);
	fprintf(fp_log, "shell width_fs/rs = %e, %e\n", del_fs, del_rs);

	Nrs = zn_rs;
	Nfs = zn_fs;

	znht_fs = del_fs / zn_fs; //In the comoving frame.
	znht_rs = del_rs / zn_rs;

	// Volume of a FS zone in terms of dl
	Voldl_fs = (SQR(R_z) * znht_fs * 0.25) / SQR(dl);
	Volcyl_fs = pi * SQR(R_z) * znht_fs;
	fprintf(fp_log, "znht_fs = %e, Vol_cyl = %e, Vol_dl = %e\n", znht_fs,
			Volcyl_fs, Voldl_fs);

	tph_esc_fs = esctimecalc(R_z, znht_fs);
	lmean_fs = tph_esc_fs * c;
	P_up_fs = (0.5 * R_z) / (R_z + znht_fs);
	P_side_fs = znht_fs / (R_z + znht_fs);
	fprintf(fp_log, "tph_esc_fs = %e, P_up_fs = %e, P_side_fs = %e, lmean_fs = %e\n",
		tph_esc_fs, P_up_fs, P_side_fs, lmean_fs);

	// Volume of a RS zone in terms of dl
	Voldl_rs = (SQR(R_z) * znht_rs * 0.25) / SQR(dl);
	Volcyl_rs = pi * SQR(R_z) * znht_rs;
	fprintf(fp_log, "zn_ht_rs = %e, Vol_cyl_rs = %e, Vol_dl_rs = %e\n", znht_rs,
			Volcyl_rs, Voldl_rs);

	tph_esc_rs = esctimecalc(R_z, znht_rs);
	lmean_rs = tph_esc_rs * c;
	P_up_rs = (0.5 * R_z) / (R_z + znht_rs);
	P_side_rs = znht_rs / (R_z + znht_rs);
	fprintf(fp_log, "tph_esc_rs = %e, P_up_rs = %e, P_side_rs = %e, lmean_rs = %e\n",
		tph_esc_rs, P_up_rs, P_side_rs, lmean_rs);

	Ntot = Nfs + Nrs;

	if (Ntot > ZON) {

	  printf("Either reduce the number of zones, Ntot, or increase the array size of all the arrays, i.e increase ZON\n");
	  fclose(fp_log);
	  exitSimulation (start, logfile, -1, sOutputDirectory, nfnfile);
	  
	} else {

	  fprintf(fp_log, "Total number of zones Ntot = %d\n", Ntot);
	}

	// Position of the collision of the two shells, in the lab frame.
	//
	z_c = colrad(Gammai, Gammao, t_v, R_o, R_i, del_o);
	fprintf(fp_log, "In lab frame, collision radius = %e\n", z_c);

	// In comoving frame of the shock
	tcr_fs_zn = znht_fs / (c * beta_fs);
	tcr_rs_zn = znht_rs / (c * beta_rs);
	fprintf(fp_log, "tcr_fs_zn = %e, tcr_rs_zn = %e\n", tcr_fs_zn, tcr_rs_zn);

	if ((tcr_fs_zn < 1.0e-1) || (tcr_rs_zn < 1.0e-1)) {

		printf("\nError, time value shouldn't be so less\n");
		fclose(fp_log);
		exitSimulation (start, logfile, -1, sOutputDirectory, nfnfile);
	}

	// In comoving frame of the shock
	tcr_fs_com = (znht_fs * zn_fs) / (c * beta_fs);
	tcr_rs_com = (znht_rs * zn_rs) / (c * beta_rs);
	fprintf(fp_log, "FS & RS crossing time = %e %e\n", tcr_fs_com, tcr_rs_com);

	B_rs = sqrt((8.0 * pi * epsilbfield * U_rs));
	B_fs = sqrt((8.0 * pi * epsilbfield * U_fs));
	fprintf(fp_log, "B_fs = %e, B_rs = %e\n", B_fs, B_rs);	

	if ((B_fs < 1.0e-3) || (B_rs < 1.0e-3)) {

	  printf("Error, magnetic field value can't be that low for inner jet\n");
	  fclose(fp_log);
	  exitSimulation (start, logfile, -1, sOutputDirectory, nfnfile);
	}

	/*
	Bdisord_fs = b_disord * B_fs;
	Bdisord_rs = b_disord * B_rs;
	Bord_fs = (1.0 - b_disord) * B_fs;
	Bord_rs = (1.0 - b_disord) * B_rs;
	*/

	gamma_max_rs = 4.654954e+07 * sqrt((alpha / B_rs));
	gamma_max_fs = 4.654954e+07 * sqrt((alpha / B_fs));
	fprintf(fp_log, "gamma_max_fs = %e, gamma_max_rs = %e\n", gamma_max_fs,
			gamma_max_rs);

	gmax_upper_fs = 4.654954e+07 / sqrt(B_fs); // Value from sqrt(3e/SIGMAT).
	gmax_upper_rs = 4.654954e+07 / sqrt(B_rs);
	fprintf(fp_log, "gamma_upper_fs = %e, gamma_upper_rs = %e\n", gmax_upper_fs,
			gmax_upper_rs);

	if ((gamma_max_fs > gemax) || (gamma_max_rs > gemax)) {

	  printf("Error, gamma_max values for RS/FS region should be lower than the highest value of the electron energy grid\n");
	  fclose(fp_log);
	  exitSimulation (start, logfile, -1, sOutputDirectory, nfnfile);
	}

	if ((gamma_max_fs > gmax_upper_fs) || (gamma_max_rs > gmax_upper_rs)) {

	  printf("Error, acceleration timescale should be smaller than the syn. cooling timescale\n");
	  fclose(fp_log);
	  exitSimulation (start, logfile, -1, sOutputDirectory, nfnfile);
	}

	amax_fs = pow(gamma_max_fs, (1.0 - q)) * (gamma_max_fs - coeffq_fs);
	bmax_fs = log(gamma_max_fs) + (coeff_fs / gamma_max_fs);
	dmax_fs = gamma_max_fs - (coeff_fs * log(gamma_max_fs));
	gamma_b_fs = myrtnewtq(funcalq, g0, xacc, q, coeff_fs, coeffq_fs, dmax_fs,
			bmax_fs, amax_fs);

	amax_rs = pow(gamma_max_rs, (1.0 - q)) * (gamma_max_rs - coeffq_rs);
	bmax_rs = log(gamma_max_rs) + (coeff_rs / gamma_max_rs);
	dmax_rs = gamma_max_rs - (coeff_rs * log(gamma_max_rs));
	gamma_b_rs = myrtnewtq(funcalq, g0, xacc, q, coeff_rs, coeffq_rs, dmax_rs,
			       bmax_rs, amax_rs);
	
	fprintf(fp_log, "gamma_b_fs = %e, gamma_b_rs = %e\n", gamma_b_fs,
		gamma_b_rs);
	
	if ((gamma_b_fs < ge) || (gamma_b_rs < ge) || (gamma_b_fs > gamma_max_fs)
	    || (gamma_b_rs > gamma_max_rs)) {
	  printf("Error, gamma1 value is either too low or too high\n");
	  fclose(fp_log);
	  exitSimulation (start, logfile, -1, sOutputDirectory, nfnfile);
	}

	rL_fs = sqrt((SQR(gamma_max_fs) - 1)) * Energy / (ECHARGE * B_fs);
	rL_rs = sqrt((SQR(gamma_max_rs) - 1)) * Energy / (ECHARGE * B_rs);
	fprintf(fp_log, "Larmor radii, fs = %e, rs = %e\n", rL_fs, rL_rs);

	if ((rL_fs > znht_fs) || (rL_rs > znht_rs)) {

	  printf("Error, larmor radius of the fastest particle should be smaller than the zone width in both regions\n");
	  fclose(fp_log);
	  exitSimulation (start, logfile, -1, sOutputDirectory, nfnfile);
	}


	// Calculate the angle between the magnetic field and the aberrated 
	// line of sight in the plasma frame in order to modulate the 
	// synchrotron and subsequently the SSC emission accordingly. 
	//
	thetaobs_rad = pi * Theta_obs / 180.0;
	thetaxy_rad = pi * thetaxy / 180.0;  // plasma frame value
	thetaz_rad = pi * thetaz / 180.0;    // plasma frame value
	printf("\nthetaobs_rad = %e, thetaxy_rad = %e, thetaz_rad = %e\n", 
	       thetaobs_rad, thetaxy_rad, thetaz_rad);

	// plasma frame value
	//
	if (Bflag == 0.0) {

	  printf("\n Randomly oriented field requested\n");

	} else if (Bflag == 1.0) {
	  //sinezi = parallelgeom(D, thetaobs_rad);
	  //etassczival = sqrt(1.0 - SQR(sinezi));
	  parallelgeom(D, thetaobs_rad, Gamma_sh, Beta_sh, &sinezi, &etassczival);
	  sinezifeed = 0.0;
	   
	  printf("\nParallel geometry sinezi = %e, coszi = %e, sinezifeed = %e\n", 
		 sinezi, etassczival, sinezifeed);

	} else if (Bflag == 2.0) {
	  //sinezi = transgeom(D, thetaobs_rad, thetaxy_rad);
	  //etassczival = sqrt(1.0 - SQR(sinezi));
	  transgeom(D, thetaobs_rad, thetaxy_rad, &sinezi, &etassczival);
	  sinezifeed = 1.0;

	  printf("\nTransverse geometry sinezi = %e, coszi = %e, sinezifeed = %e\n", 
		 sinezi, etassczival, sinezifeed);

	} else if (Bflag == 3.0) {
	  //sinezi = obliquegeom(D, thetaobs_rad, thetaxy_rad, thetaz_rad, Gamma_sh, Beta_sh);
	  //etassczival = sqrt(1.0 - SQR(sinezi));
	  obliquegeom(D, thetaobs_rad, thetaxy_rad, thetaz_rad, Gamma_sh, Beta_sh, &sinezi, &etassczival);
	  sinezifeed = sqrt(1.0 - SQR(cos(thetaz_rad)));

	  printf("\nOblique geometry sinezi = %e, coszi = %e, sinezifeed = %e\n", 
		 sinezi, etassczival, sinezifeed);

	} else if (Bflag == 4.0) {
	  toroidgeom(D, thetaobs_rad, sineord, etassczipval);
	  sinezifeed = 1.0;

	  printf("\nToroidal geometry with QUANTIZE = %d, sinezifeed = %e\n", 
		 QUANTIZE, sinezifeed);

	} else if (Bflag == 5.0) {	  
	  helixgeom(D, Gamma_sh, Beta_sh, thetaobs_rad, thetaz_rad, sineord, etassczipval);
	  sinezifeed = sqrt(1.0 - SQR(cos(thetaz_rad)));

	  printf("\nHelical geometry with QUANTIZE = %d, sinezifeed = %e\n", 
		 QUANTIZE, sinezifeed);

	} else {

	  printf("\nAny other orientation of the B-field not supported in this model\n");

	}
	
	if (Bflag == 4.0) {
	  
	  for (j = 0; j < QUANTIZE; j++) {
	    
	    //etassczipval[j] = sqrt(1.0 - SQR(toroidsine[j]));
	    printf("Toroidal geometry sinezi[%d] = %e, coszi[%d] = %e\n", 
		   j, sineord[j], j, etassczipval[j]);	    
	    
	  }

	} else if (Bflag == 5.0) {

	  for (j = 0; j < QUANTIZE; j++) {

	    //etassczipval[j] = sqrt(1.0 - SQR(helixsine[j]));
	    printf("Helical geometry sinezi[%d] = %e, coszi[%d] = %e\n", 
		   j, sineord[j], j, etassczipval[j]);

	  }

	}


	// Calculate the characteristic frequency of gamma_min to normalize the
	// synchrotron power obtained from the approximate expression at frequencies
	// below nu_gamma1.
	//
	// Corresponding to a completely randomly oriented field:
	nu_gamma1_fsrand = NU_C * B_fs * SQR(gamma_b_fs);
	nu_gamma1_rsrand = NU_C * B_rs * SQR(gamma_b_rs);
	fprintf(fp_log, "nu_gamma1_fsrand = %e, nu_gamma1_rsrand = %e\n", 
		nu_gamma1_fsrand, nu_gamma1_rsrand);

	// Corresponding to the ordered component of the field:
	if ((Bflag != 0.0) && (Bflag != 4.0) && (Bflag != 5.0)) {

	  nu_gamma1_fs = NU_C * B_fs * sinezi * SQR(gamma_b_fs);
	  nu_gamma1_rs = NU_C * B_rs * sinezi * SQR(gamma_b_rs);

	  nu_gamma1_fsfeed = NU_C * B_fs * sinezifeed * SQR(gamma_b_fs);
	  nu_gamma1_rsfeed = NU_C * B_rs * sinezifeed * SQR(gamma_b_rs);

	  fprintf(fp_log, "sinezi = %e, nu_gamma1_fs = %e, nu_gamma1_rs = %e\n", 
		  sinezi, nu_gamma1_fs, nu_gamma1_rs);

	  fprintf(fp_log, "Bflag = %e, sinezifeed = %e, nu_gamma1_fsfeed = %e, nu_gamma1_rsfeed = %e\n",
                  Bflag, sinezifeed, nu_gamma1_fsfeed, nu_gamma1_rsfeed);

	} else if ((Bflag == 4.0) || (Bflag == 5.0)) {

	  for (j = 0; j < QUANTIZE; j++) {

	    nu_gam1fs_ord[j] = NU_C * B_fs * sineord[j] * SQR(gamma_b_fs);
	    nu_gam1rs_ord[j] = NU_C * B_rs * sineord[j] * SQR(gamma_b_rs);

	    fprintf(fp_log, "ordersine[%d] = %e, nu_gam1fs_order[%d] = %e, nu_gam1rs_order[%d] = %e\n",
		    j, sineord[j], j, nu_gam1fs_ord[j], j, nu_gam1rs_ord[j]);      
	  }

          nu_gamma1_fsfeed = NU_C * B_fs * sinezifeed * SQR(gamma_b_fs);
          nu_gamma1_rsfeed = NU_C * B_rs * sinezifeed * SQR(gamma_b_rs);	  

	  fprintf(fp_log, "Bflag = %e, sinezifeed = %e, nu_gamma1_fsfeed = %e, nu_gamma1_rsfeed = %e\n",
                  Bflag, sinezifeed, nu_gamma1_fsfeed, nu_gamma1_rsfeed);
	  
	} else {

	  nu_gamma1_fs = nu_gamma1_fsrand;
	  nu_gamma1_rs = nu_gamma1_rsrand;

	}


	// Calculate synchrotron break energy for electrons in FS and RS regions
	//
	nu_B_fs = 1.29095e-9 * SQR(B_fs);
	gamma_break_fs = 1.0 / (te_esc * nu_B_fs);

	nu_B_rs = 1.29095e-9 * SQR(B_rs);
	gamma_break_rs = 1.0 / (te_esc * nu_B_rs);
	fprintf(fp_log, "syn. gam_break_fs = %e, gam_break_rs = %e\n",
			gamma_break_fs, gamma_break_rs);

	nui = nu1;
	fprintf(fp_log, "nui = %e\n", nui);

	steps = pow((gemax / ge), 1.0 / ((double) (EGRID - 1)));
	stepsnu = pow((nu2 / nu1), 1.0 / ((double) (FREQ_GRID - 1)));
	fprintf(fp_log, "steps = %e, stepsnu = %e\n", steps, stepsnu);

	// Defining particle energy grid in the following loop
	//
	for (i = 0; i < EGRID; i++) {

		gamma_e[i] = ge;
		ge *= steps;

	}

	// Defining frequency and epsilon (h*nu/(m_e*c^2)) array
	//
	for (i = 0; i < FREQ_GRID; i++) {

		syn_nu[i] = nui;
		epsil[i] = 8.082107e-21 * syn_nu[i]; // Value from h/(m_e*c^2)
		epsilscatter[i] = epsil[i];
		nui *= stepsnu;
		fprintf(fp_freq, "%e %e %e\n", syn_nu[i], epsil[i], epsilscatter[i]);
	}

	fclose(fp_freq);

	// Determine injection function, as given in ApJ, 2002, 581, 127, for both
	// forward and reverse shocks in the comoving frame.
	//
	num_fs = epsilelec * U_fs / (zn_fs * tcr_fs_zn);
	num_rs = epsilelec * U_rs / (zn_rs * tcr_rs_zn);

	if (q == 2.0) {

		Q0_fs = num_fs / (Energy * log(gamma_max_fs / gamma_b_fs));
		Q0_rs = num_rs / (Energy * log(gamma_max_rs / gamma_b_rs));

	} else {

	  Q0_fs = (num_fs * (2.0 - q))
	    / (Energy
	       * ((pow(gamma_max_fs, (2.0 - q)))
		  - (pow(gamma_b_fs, (2.0 - q)))));

	  Q0_rs = (num_rs * (2.0 - q))
	    / (Energy
	       * ((pow(gamma_max_rs, (2.0 - q)))
		  - (pow(gamma_b_rs, (2.0 - q)))));

	}

	fprintf(fp_log, "num_fs = %e, Q0_fs = %e, num_rs = %e, Q0_rs = %e\n",
			num_fs, Q0_fs, num_rs, Q0_rs);

	// Defining particle injection spectrum.
	//
	for (j = 0; j < EGRID; j++) {

		if (gamma_e[j] < gamma_b_fs) {

		  Q_inj_fs[j] = SQR((gamma_e[j]/gamma_b_fs)) / pow(gamma_b_fs, q);

		} else if (gamma_e[j] < gamma_max_fs) {

		  Q_inj_fs[j] = 1.0 / pow(gamma_e[j], q);

		} else {

		  Q_inj_fs[j] = 1.0
		    / (pow(gamma_max_fs, q)
		       * SQR((gamma_e[j] / gamma_max_fs))
		       * SQR((gamma_e[j] / gamma_max_fs)));
		}

		Q_inj_fs[j] *= Q0_fs * sqrt((SQR(gamma_e[j]) - 1.0)) / gamma_e[j];
		fprintf(fp_fsinj, "%e %e\n", gamma_e[j], Q_inj_fs[j]);

		if (gamma_e[j] < gamma_b_rs) {
		  
		  Q_inj_rs[j] = SQR((gamma_e[j]/gamma_b_rs)) / pow(gamma_b_rs, q);
		  
		} else if (gamma_e[j] < gamma_max_rs) {

		  Q_inj_rs[j] = 1.0 / pow(gamma_e[j], q);
		  
		} else {
		  
		  Q_inj_rs[j] = 1.0
		    / (pow(gamma_max_rs, q)
		       * SQR((gamma_e[j] / gamma_max_rs))
		       * SQR((gamma_e[j] / gamma_max_rs)));
		}

		Q_inj_rs[j] *= Q0_rs * sqrt((SQR(gamma_e[j]) - 1.0)) / gamma_e[j];
		fprintf(fp_rsinj, "%e %e\n", gamma_e[j], Q_inj_rs[j]);

	}

	fclose(fp_fsinj);
	fclose(fp_rsinj);


	// Calculate the appropriate time delay for each zone.
	//
	if (mu_prime > 0.0) {

		for (zn = Nrs; zn < Ntot; zn++) {

			tdelsi_fs[zn] = (((double) (Ntot - zn - 1)) * znht_fs * mu_prime)
							/ c;
			/* In the comoving frame */
			fprintf(fp_log, "\ntdelsi_fs[%d] = %e, no. = %d\n", zn,
					tdelsi_fs[zn], (Ntot - zn - 1));

		}

		for (zn = (Nrs - 1); zn >= 0; zn--) {

			tdelsi_rs[zn] = (((((double) (Nrs - zn)) * znht_rs)
					  + (znht_fs * zn_fs)) * mu_prime) / c;
			/* In the comoving frame */
			fprintf(fp_log, "\ntdelsi_rs[%d] = %e, no. = %d\n", zn,
					tdelsi_rs[zn], (Nrs - zn - 1));

		}

	} else {

		mu_prime_abs = fabs(mu_prime);
		printf("Out of superluminal cone, mu, %e < beta_sh, %e\n", mu, Beta_sh);
		printf("mu_prime_abs = %e\n", mu_prime_abs);

		for (zn = (Nrs - 1); zn >= 0; zn--) {

			tdelsi_rs[zn] = (((double) (zn)) * znht_rs * mu_prime_abs) / c;
			/* In the comoving frame */
			fprintf(fp_log, "\ntdelsi_rs[%d] = %e, no. = %d\n", zn,
					tdelsi_rs[zn], zn);

		}

		for (zn = Nrs; zn < Ntot; zn++) {

			tdelsi_fs[zn] = (((((double) (zn - Nfs + 1)) * znht_fs)
					+ (znht_rs * zn_rs)) * mu_prime_abs) / c;
			/* In the comoving frame */
			fprintf(fp_log, "\ntdelsi_fs[%d] = %e, no. = %d\n", zn,
					tdelsi_fs[zn], (zn - Nfs));

		}

	}


	 //Calculate distance of emission region from 0.03 - 3 pc in 100
	 //logarithmic grid points.
	 //
	 logstepz = pow((3.0 / 0.03), 1.0/((double) (ED - 1)));
	 emitdist_min = 9.25704e16;

	 for (i = 0; i < ED; i++) {

	 emit_dist[i] = emitdist_min;
	 Izero[i] = L_BLR / (JCOEFF * SQR(emit_dist[i]));

	 emitdist_min *= logstepz;

	 printf("zht[%d] = %e, Izero[%d] = %e\n",
	 i, emit_dist[i], i, Izero[i]);

	 }


#ifdef BLRDT_UPH

	 // Calculate etaph array from -1 to 1 for 4000 points to be used only
	 // for obtaining BLR intensity profiles
	 //
	 stepetaph = 2.0 / ((double) (4000));
	 etaph_linemitmin = -1.0;

	 for (i = 0; i < 4000; i++) {
	   
	   etaph_linemit[i] = etaph_linemitmin;
	   etaphstar_linemit[i] = (etaph_linemit[i] + Beta_sh)
	     / (1.0 + (Beta_sh * etaph_linemit[i]));

	   etaph_linemitmin += stepetaph;

	 }

#endif

	 //Calculate eta array to be used in the program for BLR calculations
	 //
	 fp_2401 = fopen("abs_weight_24-10.dat", "w");
	 if (!fp_2401) {
	   
	   printf("\nCould not open the output file abs_weight_2401. \n");
	   return -1;
	 }

	// Calculate the abscissas and weights of the Guass-Legendre n-point quadrature
	// formula, for given lower and upper limits and n
	//
	gauleg(-1.0, 0.0, absarr-1, weights-1, 24);

	for (i = 0; i < MAXPOINT; i++) {
		fprintf(fp_2401, "%e %e\n", absarr[i], weights[i]);

	}

	fclose(fp_2401);

	integ_sample = qgaus(integ_exam, 0.0, 1.0, absarr-1, weights-1);
	printf("integ_sample = %e\n", integ_sample);

	//Calculate eta points from 1 to 0 on a linear grid
	//
	fp_eta = fopen("eta_points_11.dat", "w");
	if (!fp_eta) {

		printf("\nCould not open the output file eta_points_11. \n");
		return -1;
	}

	lingrid(0.0, 1.0, etaph_lin);

	for (i = 0; i < ETASTEP; i++) {
		
		etaphstar_lin[i] = (etaph_lin[i] + Beta_sh)
						/ (1.0 + (Beta_sh * etaph_lin[i]));
		printf("etaph_lin[%d] = %e, etaphstar_lin[%d] = %e\n", i, etaph_lin[i],
				i, etaphstar_lin[i]);

		fprintf(fp_eta, "%e %e\n", etaph_lin[i], etaphstar_lin[i]);

	}

	fclose(fp_eta);

	gausgrid(-1.0, 0.0, 24, absarr-1, etaph_gausp, etaph_gausm, &GAUSASZ);

	printf("GAUSASZ = %d, 2MAXPOINT = %d\n", GAUSASZ, (2 * MAXPOINT));

	if ((MAXPOINT) < GAUSASZ) {

		printf("GAUSASZ = %d should be less than 2 * MAXGRID = %d\n", GAUSASZ,
				MAXPOINT);

		exitSimulation (start, logfile, -1, sOutputDirectory, nfnfile);
	}

	for (i = 0; i < GAUSASZ; i++) {

		etaphstar_gausp[i] = (etaph_gausp[i] + Beta_sh)
						/ (1.0 + (Beta_sh * etaph_gausp[i]));

		etaphstar_gausm[i] = (etaph_gausm[i] + Beta_sh)
						/ (1.0 + (Beta_sh * etaph_gausm[i]));
		printf("etaph_gausp[%d] = %e, etaphstar_gausp[%d] = %e\n", i,
				etaph_gausp[i], i, etaphstar_gausp[i]);

		printf("etaph_gausm[%d] = %e, etaphstar_gausm[%d] = %e\n", i,
				etaph_gausm[i], i, etaphstar_gausm[i]);

	}

	//Declaring phi grid to be used in disk, BLR, & dusty torus calculations.
	//
	phimin = 0.0;
	phistep = (pi - phimin) / ((double) (PHIGRID));

	for (i = 0; i < PHIGRID; i++) {

	    phi_main[i] = phimin;

	    printf("phi_main[%d] = %e\n", i, phi_main[i]);

	    phimin += phistep;

	}

	printf("\n");


	// Calculate disk temperature for various disk radii, both in lab frame,
	// for use in ECD calculations
	//
	if (Ldisk != 0.0) {

		stepepdisk = pow((EPDMAX / EPDMIN), 1.0 / ((double) (DISKGRID - 1)));
		diskepmin = EPDMIN;

		for (i = 0; i < DISKGRID; i++) {

			disklabep[i] = diskepmin;
			printf("disklabep[%d] = %e\n", i, disklabep[i]);

			diskepmin *= stepepdisk;
		}

		fp_theta = fopen(thetafile, "w");
		if (!fp_theta) {

			printf("\nCould not open the output file disk_theta. \n");
			return -1;
		}

		for (i = 0; i < RGRID; i++) {
			disktemp(MBH, Ldisk, etacc, raddist, thetar);
			fprintf(fp_theta, "%e %e\n", raddist[i], thetar[i]);

		}

		fclose(fp_theta);

	}


	// Calculate BLR photon density for various angles, both in comoving
	// frame, for use in ECBLR calculations
	//
	if (L_BLR != 0.0) {

		fp_blrj = fopen(blrjfile, "w");
		if (!fp_blrj) {

			printf("\nCould not open the output file blr_j. \n");
			return -1;
		}

		fp_blrin = fopen("blr_lines.txt", "r");
		if (!fp_blrin) {

			printf("\nCould not open the input BLR line file. \n");
			return -1;
		}

		for (j = 0; ((j < LINEGRID) && (!feof(fp_blrin))); j++) {

			fscanf(fp_blrin, "%lf %lf\n", &linelablam[j], &lineNV[j]);

			linelabfreq[j] = c / (linelablam[j] * 1e-8);
			linelabep[LINEGRID - j - 1] = (PLANCKCONST * linelabfreq[j])
							/ Energy;

			printf("linelablam[%d] = %e, linelabep[%d] = %e\n",
			       j, linelablam[j], (LINEGRID-j-1), linelabep[LINEGRID-j-1]);

		}

		fclose(fp_blrin);

		for (j = 0; j < (LINEGRID / 2); j++) {

			temp1 = lineNV[j];
			temp2 = lineNV[LINEGRID - j - 1];
			lineNV[LINEGRID - j - 1] = temp1;
			lineNV[j] = temp2;

		}

		printf("\n");

		for (j = 0; j < LINEGRID; j++) {

			nphlinelab[j] = LINECOEFF * lineNV[j] / linelabep[j];

			printf("linelabep[%d] = %e, nphlinelab[%d] = %e, lineNV[%d] = %e\n",
					j, linelabep[j], j, nphlinelab[j], j, lineNV[j]);
		}

		contsteps = pow((CONTMAX / CONTMIN), 1.0 / ((double) (CONTGRID - 1)));
		contepmin = CONTMIN;

		for (i = 0; i < CONTGRID; i++) {

			contlabep[i] = contepmin;
			contepmin *= contsteps;

			printf("contlabep[%d] = %e\n", i, contlabep[i]);
		}

		eplabinteg_main = epsummain(Gamma_sh, 0.0, contlabep);

		printf("eplabinteg_main = %e\n", eplabinteg_main);

		for (i = 0; i < CONTGRID; i++) {

			expterm_main = contlabep[i] / AVTEMP;

			if (expterm_main >= 150.0) {
			  nphcontlab[i] = 0.0;

			} else if ((expterm_main >= 24.0) && (expterm_main < 150.0)) {
			  contd1_main = exp(expterm_main);
			  nphcontlab[i] = CONTCOEFF * SQR(contlabep[i]) /
			    (contd1_main * eplabinteg_main);

			} else if ((expterm_main < 24.0) && (expterm_main > 1.0e-4)) {
			  contd1_main = exp(expterm_main) - 1.0;
			  nphcontlab[i] = CONTCOEFF * SQR(contlabep[i]) /
			    (contd1_main * eplabinteg_main);

			} else {
			  contd1_main = expterm_main;
			  nphcontlab[i] = CONTCOEFF * SQR(contlabep[i]) /
			    (contd1_main * eplabinteg_main);

			}

			//printf("contlabep[%d] = %e, nphcontlab[%d] = %e\n", i, contlabep[i], i, nphcontlab[i]);

		}


#ifdef BLRDT_UPH

		//For calculation of uph BLR purposes only
		//
		for (i = 0; i < 4000; i++) {

		  bphterm_main = 1.0 + (Beta_sh * etaph_linemit[i]);

		  for (k = 0; k < LINEGRID; k++) {

		    lineep[i][k] = linelabep[k] / (Gamma_sh * bphterm_main);

		  }

		  for (k = 0; k < CONTGRID; k++) {

		    contep[i][k] = contlabep[k] / (Gamma_sh * bphterm_main);
		  }

		}

#endif

		for (i = 0; i < GAUSASZ; i++) {

			bphtermp = 1.0 + (Beta_sh * etaph_gausp[i]);

			for (k = 0; k < LINEGRID; k++) {
				lineepp[i][k] = linelabep[k] / (Gamma_sh * bphtermp);

				//printf("lineepp[%d][%d] = %e, linelabep[%d] = %e\n", i, k, lineepp[i][k], k, linelabep[k]);

			}

			for (k = 0; k < CONTGRID; k++) {
				contepp[i][k] = contlabep[k] / (Gamma_sh * bphtermp);

				//printf("contepp[%d][%d] = %e, contlabep[%d] = %e\n", i, k, contepp[i][k], k, contlabep[k]);

			}

			bphtermm = 1.0 + (Beta_sh * etaph_gausm[i]);

			for (k = 0; k < LINEGRID; k++) {
				lineepm[i][k] = linelabep[k] / (Gamma_sh * bphtermm);

				//printf("lineepm[%d][%d] = %e, linelabep[%d] = %e\n", i, k, lineepm[i][k], k, linelabep[k]);

			}

			for (k = 0; k < CONTGRID; k++) {
				contepm[i][k] = contlabep[k] / (Gamma_sh * bphtermm);

				//printf("contepm[%d][%d] = %e, contlabep[%d] = %e\n", i, k, contepm[i][k], k, contlabep[k]);

			}

		}

		for (i = 0; i < ETASTEP; i++) {

			bphterm_lin = 1.0 + (Beta_sh * etaph_lin[i]);
			for (k = 0; k < LINEGRID; k++) {
				lineeplin[i][k] = linelabep[k] / (Gamma_sh * bphterm_lin);

				//printf("lineeplin[%d][%d] = %e, linelabep[%d] = %e\n", i, k, lineeplin[i][k], k, linelabep[k]);

			}

			for (k = 0; k < CONTGRID; k++) {
				conteplin[i][k] = contlabep[k] / (Gamma_sh * bphterm_lin);

				//printf("conteplin[%d][%d] = %e, contlabep[%d] = %e\n", i, k, conteplin[i][k], k, contlabep[k]);

			}

		}


		for (i = 0; i < PHIGRID; i++) {

			for (j = 0; j < GAUSASZ; j++) {

			  coskip_main[i][j] = (etaph_gausp[j] * mu_prime)
			    + (sqrt(
				    ((1.0 - SQR(etaph_gausp[j]))
				     * (1.0 - SQR(mu_prime))))
			       * cos(phi_main[i]));
			  
			  if (coskip_main[i][j] < -0.99999999)
			    coskip_main[i][j] = -0.99999999;
			  if (coskip_main[i][j] > 0.99999999)
			    coskip_main[i][j] = 0.99999999;
			  
			  coskim_main[i][j] = (etaph_gausm[j] * mu_prime)
			    + (sqrt(
				    ((1.0 - SQR(etaph_gausm[j]))
				     * (1.0 - SQR(mu_prime))))
			       * cos(phi_main[i]));

			  if (coskim_main[i][j] < -0.99999999)
			    coskim_main[i][j] = -0.99999999;
			  if (coskim_main[i][j] > 0.99999999)
			    coskim_main[i][j] = 0.99999999;

			}

			for (j = 0; j < ETASTEP; j++) {
			  
			  coskilin_main[i][j] = (etaph_lin[j] * mu_prime)
			    + (sqrt(
				    ((1.0 - SQR(etaph_lin[j]))
				     * (1.0 - SQR(mu_prime))))
			       * cos(phi_main[i]));

			  if (coskilin_main[i][j] < -0.99999999)
			    coskilin_main[i][j] = -0.99999999;
			  if (coskilin_main[i][j] > 0.99999999)
			    coskilin_main[i][j] = 0.99999999;
			  
			}
			
		}


		// Calculate BLR emissivity for line and continuum emission.
		//
		blremissiv(L_BLR, tau_BLR, fcov, Rin_BLR, Rout_BLR, blrrad,
				jlineemissiv, jcontemissiv, &lineinteg, &continteg);
		for (i = 0; i < GS; i++) {

			fprintf(fp_blrj, "%e %e %e\n", blrrad[i], jlineemissiv[i],
					jcontemissiv[i]);

		}

		fclose(fp_blrj);
		//printf("lineinteg = %e, continteg = %e\n", lineinteg, continteg);


#ifdef BLRDT_UPH

		// Calculate BLR intensity profiles using 4000 etaph points for 
		// illustration purposes only.
		//
                sprintf(filename_blrupnph, "%s%s", blrupnph_file, ".dat");
                fp_blrupnph = fopen(filename_blrupnph, "w");
                if (!fp_blrupnph) {
                  printf("\nCould not open the output file BLR upnph. \n");
                  return -1;
                }

		for (i = 0; i < ED; i++) {
		  
		  sprintf(filename_blrIz, "%s%d%s", blrIzfile, i, ".dat");
		  fp_blrIz = fopen(filename_blrIz, "w");
		  if(!fp_blrIz) {
		    
		    printf("\nCould not open the output file blr_Iz. \n");
		    return -1;
		  }

		  for (j = 0; j < 4000; j++) {
		    
		    blrintens(L_BLR, tau_BLR, fcov, emit_dist[i], Rin_BLR, Rout_BLR, etaph_linemit[j], Gamma_sh, lineinteg, continteg, &lineI, &contI);
		    fprintf(fp_blrIz, "%e %e %e %e %e %e\n", etaph_linemit[j], etaphstar_linemit[j], 
			    max(1.0e-40, (lineI / Izero[i])), max(1.0e-40, (contI / Izero[i])), max(1.0e-40, lineI), max(1.0e-40, contI));
		    
		  }

		  fclose(fp_blrIz);

		  /*
                  nphblrtotuph(lineNV, Gamma_sh, lineinteg, continteg, L_BLR, tau_BLR, fcov, emit_dist[i], Rin_BLR, Rout_BLR, &lineI, &contI, lineep, contep, etaph_linemit, nph_cont, nph_line, eplabinteg_main, blrIznphfile, lineinten_main, continten_main, i);
		  */

		  nphblrtot(lineNV, Gamma_sh, lineinteg, continteg, L_BLR, tau_BLR, fcov, emit_dist[i], Rin_BLR, Rout_BLR, -1.0, 0.0, 0.0, 1.0, &lineI, &contI, lineepp, contepp, lineepm, contepm, lineeplin, conteplin, etaph_lin, etaph_gausp, etaph_gausm, absarr-1, weights-1, nph_contlin, nph_linelin, nph_contgausp, nph_linegausp, nph_contgausm, nph_linegausm, eplabinteg_main, blrIznphfile, lineinten_linmain, continten_linmain, lineinten_gauspmain, continten_gauspmain, lineinten_gausmmain, continten_gausmmain, i, GAUSASZ);

		  // Calculate total photon energy density as a function of distance 
		  // for the BLR.
		  //
		  for (j = 0; j < GAUSASZ; j++) {

		    uphblrster_contgausp = ne_num(contepp[j], nph_contgausp[j], CONTGRID);
		    uphblrster_linegausp = ne_num(lineepp[j], nph_linegausp[j], LINEGRID);
		    
		    uphblrster_contgausm = ne_num(contepm[j], nph_contgausm[j], CONTGRID);
		    uphblrster_linegausm = ne_num(lineepm[j], nph_linegausm[j], LINEGRID);
		    
		    uphblrcom_gauscont[i] += (weights[j] * (uphblrster_contgausp + uphblrster_contgausm));
		    uphblrcom_gausline[i] += (weights[j] * (uphblrster_linegausp + uphblrster_linegausm));

		  }

		  uph_blrcomgaus[i] = uphblrcom_gauscont[i] + uphblrcom_gausline[i];
		  uph_blrcomgaus[i] *= 0.5;

		  for (l = 0; l < ETASTEP; l++) {

		    //for (l = 0; l < 4000; l++) { 
		    
		    uphblrster_conteta[i][l] = ne_num(conteplin[l], nph_contlin[l], CONTGRID);
		    uphblrster_lineeta[i][l] = ne_num(lineeplin[l], nph_linelin[l], LINEGRID);

		    /*
                    uphblrster_conteta[i][l] = ne_num(contep[l], nph_cont[l], CONTGRID);
                    uphblrster_lineeta[i][l] = ne_num(lineep[l], nph_line[l], LINEGRID);
		    */
		    
		    uphblrster_eta[i][l] = uphblrster_conteta[i][l] + uphblrster_lineeta[i][l];

		  }

		  uph_blrcometa[i] = ne_numtot(etaph_lin, uphblrster_eta[i], ETASTEP);
		  //uph_blrcometa[i] = ne_numtot(etaph_linemit, uphblrster_eta[i], 4000);

		  uph_blrcomtot[i] = Energy * (uph_blrcometa[i] + uph_blrcomgaus[i]);
		  //uph_blrcomtot[i] = Energy * uph_blrcometa[i];

                 fprintf(fp_blrupnph, "%e %e %e %e\n", emit_dist[i], max(1.0e-40, (Energy * uph_blrcomgaus[i])),
                         max(1.0e-40, (Energy * uph_blrcometa[i])), max(1.0e-40, uph_blrcomtot[i]));

		 /*
                 fprintf(fp_blrupnph, "%e %e %e\n", emit_dist[i], max(1.0e-40, (Energy * uph_blrcometa[i])), 
			 max(1.0e-40, uph_blrcomtot[i]));
		 */

		}

		fclose(fp_blrupnph);		

#endif

	}


	//Calculate DT photon energy in lab frame and use in calculating
	//comoving DT photon energies for all incoming photon angles.
	//
	if (L_DT != 0.0) {

	  dtstepsize = pow((EPDTMAX / EPDTMIN), 1.0/((double) (DTGRID - 1)));
	  dt_epmin = EPDTMIN;

	  printf("\n");

	  for (i = 0; i < DTGRID; i++) {

	    dtlabpe[i] = dt_epmin;

	    printf("dtlabpe[%d] = %e\n", i, dtlabpe[i]);

	    dt_epmin *= dtstepsize;

	  }

	  //Calculate the DT area for a given disk fraction and DT covering 
	  //fraction
	  //
	  Ldt = Ldiskfrac * Ldisk * 1.0e46;
	  dt_area = Ldt / (SIGMASB * SQR(DTTEMP * DTTEMP) * FCOVDT);

	  printf("Ldt = %e, SIGMASB = %e, DTTEMP = %e, FCOVDT = %e\n", 
		 Ldt, SIGMASB, DTTEMP, FCOVDT);

	  printf("\nLdt = %e, dt_area = %e\n", Ldt, dt_area);

	
#ifdef BLRDT_UPH

	  // Calculate total photon energy density as a function of distance 
	  // for the dusty torus.
	  //
	  sprintf(filename_dtupnph, "%s%s", dtupnph_file, ".dat");
	  fp_dtupnph = fopen(filename_dtupnph, "w");
	  if (!fp_dtupnph) {
	    printf("\nCould not open the output file DT upnph. \n");
	    return -1;
	  }

	  for (k = 0; k < ED; k++) {

	    // Calculate min. and max. cosine of angles, as a function of z, for dusty
	    // torus photons between which they will enter the jet, in comoving frame.
	    //
	    etadet(Gamma_sh, emit_dist[k], &mincos, &maxcos, etadtarr);
	    
	    // Calculate the corresponding comoving energies of incoming DT photons
	    // for all angles between min. & max. cosine angles.
	    //
	    nphdt(Gamma_sh, etadtarr, dtlabpe, dtep, nph_dt, FCOVDT);

	    sprintf(filename_dtepnph, "%s%d%s", dtepnph_file, k, ".dat");
	    fp_dtepnph = fopen(filename_dtepnph, "w");
	    if (!fp_dtepnph) {
	      printf("\nCould not open the output file DT epnph. \n");
	      return -1;
	    }
	    
	    for (i = 0; i < ETADT; i++) {

	      inten_dtcom[i] = c * Energy * ne_num(dtep[i], nph_dt[i], DTGRID);
	      
	      fprintf(fp_dtepnph, "%e %e %e\n", etadtarr[i], max(1.0e-40, inten_dtcom[i]), max(1.0e-40, (inten_dtcom[i]/(c*Energy))));
	      
	    }
	    
	    fclose(fp_dtepnph);

	    for (i = 0; i < ETADT; i++) {
	      
	      uph_dtcom[k] = ne_num(dtep[i], nph_dt[i], DTGRID);
	      uph_dtcomster[k][i] = uph_dtcom[k];

	    }
	      
	    uph_dtcomtot[k] = Energy * ne_numtot(etadtarr, uph_dtcomster[k], ETADT);
	       	
	    fprintf(fp_dtupnph, "%e %e %e\n", emit_dist[k], max(1.0e-40, (Energy * uph_dtcom[k])), 
		    max(1.0e-40, uph_dtcomtot[k])); 

	  }

	  fclose(fp_dtupnph);

#endif

	}

	//exitSimulation (start, logfile, -1, sOutputDirectory, nfnfile);


	// Main loop starts here
	//
	printf("Simulation starts\n");

	printf("Start time = %e seconds\n", start);

	fprintf(fp_log, "\ntime count = %d, timecount_fs = %d, timecount_rs = %d\n",
			time_count, timecount_fs, timecount_rs);

	fprintf(fp_log, "Start time = %e seconds\n",
			(double) (start) / (double) (CLOCKS_PER_SEC));
	fclose(fp_log);

	// Calculate maximum electron energy loss rate corresponding to a 
	// randomly oriented field due to synchrotron mechanism here:
	//
	gammadotsynmax_fs = fabs(gammasyndot(gamma_max_fs, B_fs));
	gammadotsynmax_rs = fabs(gammasyndot(gamma_max_rs, B_rs));

	// Calculate electron energy loss rate for the entire electron
	// grid due to synchrotron mechanism for a randomly oriented field:
	//
	for (k = 0; k < EGRID; k++) {

	  esynloss_fs[k] = fabs(gammasyndot(gamma_e[k], B_fs));
	  esynloss_rs[k] = fabs(gammasyndot(gamma_e[k], B_rs));

	}

	volfs_param = 5.355265e-47 * Volcyl_fs;
	volrs_param = 5.355265e-47 * Volcyl_rs;

	lmeanfs_param = 6.234375e-26 * lmean_fs; // Value = 3*SIGMAT/32
	lmeanrs_param = 6.234375e-26 * lmean_rs;

	pupfs_param = P_up_fs / tph_esc_fs;
	psidefs_param = P_side_fs / tph_esc_fs;

	puprs_param = P_up_rs / tph_esc_rs;
	psiders_param = P_side_rs / tph_esc_rs;

	voldlfs_param = 5.355265e-47 * Voldl_fs * SQR((D * D)) * 1.0e23 / (1.0 + Z);
	voldlrs_param = 5.355265e-47 * Voldl_rs * SQR((D * D)) * 1.0e23 / (1.0 + Z);	


	while (simtime <= endtime && FLUX == TRUE) {       

		if (time_count == 0) {
			f = Nrs + 1; // Initializing the zone.
			accf = TRUE;

		} else if ((f < Ntot) && (simtime > time_fs)) {
			f = f + 1; // Moving to the next zone.
			timecount_fs = 0;

		} else if ((f == Ntot) && (simtime > time_fs)) {
			accf = FALSE; // When in the last zone & shock has crossed it.

		}

		if (time_count == 0) {
			r = Nrs - 1; // Initializing the zone.
			accr = TRUE;

		} else if ((r > 0) && (simtime > time_rs)) {
			r = r - 1; // Moving to the next zone.
			timecount_rs = 0;

		} else if ((r == 0) && (simtime > time_rs)) {
			accr = FALSE; // When in the last zone & shock has crossed it.

		}


		if (fmod(write_count, 10) == 0) {

			fp_log = fopen(logfile, "a");
			if (!fp_log) {

				printf("\nCould not open the output file log. \n");
				return -1;
			}

			fprintf(fp_log,
					"\nf = %d, Nfs = %d, time_fs initially = %e, timecount_fs = %d, accf = %d\n",
					f, Nfs, time_fs, timecount_fs, accf);
			fprintf(fp_log,
					"r = %d, Nrs = %d, time_rs initially = %e, timecount_rs = %d, accr = %d\n",
					r, Nrs, time_rs, timecount_rs, accr);

		}

		if (files) {

			sprintf(filename0_fs, "%s%d%s", fssafile, time_count, ".dat");
			fp_fssa = fopen(filename0_fs, "w");
			if (!fp_fssa) {

				printf("\nCould not open the output file ssa. \n");
				return -1;
			}

			sprintf(filename1_fs, "%s%d%s", fsynnfnfile, time_count, ".dat");
			fp_fsynnfn = fopen(filename1_fs, "a");
			if (!fp_fsynnfn) {

				printf("\nCould not open the output file fsynnfn. \n");
				return -1;
			}

			sprintf(filename1_rs, "%s%d%s", rsynnfnfile, time_count, ".dat");
			fp_rsynnfn = fopen(filename1_rs, "a");
			if (!fp_rsynnfn) {

				printf("\nCould not open the output file rsynnfn. \n");
				return -1;
			}

		}


		// Calculate particle injection timescale only when acceleration is on
		//
		if (accf == TRUE) {
			dtinj_fs = te_esc;
			dt_j_fs = 1.0e20;

			if (timecount_fs != 0) {
				for (j = 0; j < EGRID; j++) {

					if ((gamma_e[j] >= gamma_b_fs)
							&& (gamma_e[j] <= gamma_max_fs)) {
						dt_j_fs = ne_curnt_fs[j] / Q_inj_fs[j];

					}

					if (dt_j_fs < dtinj_fs)
						dtinj_fs = dt_j_fs;

					// In order to get the smallest possible injection timescale from
					// the entire particle density and the injection function array.

				}

			}

		}

		if (accr == TRUE) {
			dtinj_rs = te_esc;
			dt_j_rs = 1.0e20;

			if (timecount_rs != 0) {
				for (j = 0; j < EGRID; j++) {

					if ((gamma_e[j] >= gamma_b_rs)
							&& (gamma_e[j] <= gamma_max_rs)) {
						dt_j_rs = ne_curnt_rs[j] / Q_inj_rs[j];

					}

					if (dt_j_rs < dtinj_rs)
						dtinj_rs = dt_j_rs;

					// In order to get the smallest possible injection timescale from
					// the entire particle density and the injection function array.

				}

			}

		}

		
		// Calculate disk energy and the corresponding etaph and cosineki angles
		// in PF.
		//
		if ((Ldisk != 0.0) && ((LDISK_fs != 0) || (LDISK_rs != 0))) {

		  for (j = 0; j < RGRID; j++) {
		    
		    rzterm_main = 0.5 * SQR((raddist[j] / z_c));
		    if (rzterm_main < 1.0e-6) {
		      hypt[j] = z_c * (1.0 + rzterm_main);
		      etaphd[j] = (1.0 - (Beta_sh * rzterm_main) - Beta_sh) /
			(1.0 + rzterm_main - Beta_sh);
		      
		    } else {
		      xterm_main = SQR(raddist[j]) + SQR(z_c);
		      hypt[j] = sqrt(xterm_main);
		      etaphd[j] = (z_c - (Beta_sh * hypt[j])) / 
			(hypt[j] - (Beta_sh * z_c));
		      
		    }
		    
		    if (etaphd[j] < -0.99999999) etaphd[j] = -0.99999999;                                                                      
		    if (etaphd[j] > 0.99999999) etaphd[j] = 0.99999999;                                                                        
		    
		    for (i = 0; i < DISKGRID; i++) {
		      
		      p_engy[j][i] = disklabep[i] / (1.0 + (Beta_sh * etaphd[j]));
		      
		    }

		  }
		  
		}


		// Calculate photon density from BLR entering the emission region in PF.
		// for both linear and gaussian grid.
		//
		if ((L_BLR != 0.0) && ((LBLR_fs != 0) || (LBLR_rs != 0))) {

		  nphblrtot(lineNV, Gamma_sh, lineinteg, continteg, L_BLR, tau_BLR, fcov,
			    z_c, Rin_BLR, Rout_BLR, -1.0, 0.0, 0.0, 1.0, &lineI, &contI,
			    lineepp, contepp, lineepm, contepm, lineeplin, conteplin,
			    etaph_lin, etaph_gausp, etaph_gausm, absarr-1, weights-1,
			    nph_contlin, nph_linelin, nph_contgausp, nph_linegausp,
			    nph_contgausm, nph_linegausm, eplabinteg_main, blrIznphfile,
			    lineinten_linmain, continten_linmain, lineinten_gauspmain,
			    continten_gauspmain, lineinten_gausmmain, continten_gausmmain,
			    Izindex, GAUSASZ);

		}


		// Calculate min. and max. cosine of angles, as a function of z, for dusty
                // torus photons between which they will enter the jet, in comoving frame.
                //
		if ((L_DT != 0.0) && ((LDT_fs != 0) || (LDT_rs != 0))) {

		  etadet(Gamma_sh, z_c, &mincos, &maxcos, etadtarr);
		  printf("\nz_c = %e, mincos = %e, maxcos = %e\n", z_c, mincos, maxcos);

		  // Calculate the corresponding comoving energies of incoming DT photons
		  // for all angles between min. & max. cosine angles.
		  //
		  nphdt(Gamma_sh, etadtarr, dtlabpe, dtep, nph_dt, FCOVDT);

		}


		// Calculate time step here
		//
		dtcool_fs = 1.0e+10;
		for (zn = Nrs; zn < f; zn++) {

		  if ((accf == FALSE) || ((accf == TRUE) && (zn != (f - 1)))) {
		    
		    if (gam_max[zn] == 0.0) {
		      
		      gam_max[zn] = gamma_max_fs;
		      
		    }
		    
		    gammadot_synmax[zn] = fabs(gammasyndot(gam_max[zn], B_fs));
		    
		    gammadot_sscmax[zn] = fabs(gammasscdot(epsil, n_phgg[zn], gam_max[zn], FREQ_GRID));
		      
		    if ((Ldisk != 0.0) && (LDISK_fs != 0)) {
		      
		      gammadot_ecdmax[zn] = fabs(gammaecdiskdot(gam_max[zn], Gamma_sh, z_c, thetar, raddist));
		      
		    } else {

		      gammadot_ecdmax[zn] = 0.0;
		    }
		    
		    if ((L_BLR != 0.0) && (LBLR_fs != 0)) {
		      
		      gammadot_blrmax[zn] = fabs(
						 gamdotblr(gam_max[zn], -1.0, 0.0, lineepp, lineepm,
							   contepp, contepm, lineeplin, conteplin,
							   nph_linelin, nph_contlin, etaph_lin,
							   nph_linegausp, nph_linegausm, weights-1,
							   etaph_gausp, etaph_gausm, nph_contgausp,
							   nph_contgausm, GAUSASZ));
		      
		    } else {

		      gammadot_blrmax[zn] = 0.0;
		    }
		    
		    if ((L_DT != 0.0) && (LDT_fs != 0)) {
		      
		      gammadot_dtmax[zn] = fabs(gamdotdt(gam_max[zn], etadtarr, dtep, nph_dt));
		      
		    } else {

		      gammadot_dtmax[zn] = 0.0;
		    }
		    
		  } else {

		    gammadot_synmax[zn] = gammadotsynmax_fs;
		    
		    if ((Ldisk != 0.0) && (LDISK_fs != 0)) {
		      
		      gammadot_ecdmax[zn] = fabs(gammaecdiskdot(gamma_max_fs, Gamma_sh, z_c, thetar, raddist));
		      
		    } else {

		      gammadot_ecdmax[zn] = 0.0;
		    }
		    
		    if ((L_BLR != 0.0) && (LBLR_fs != 0)) {
		      
		      gammadot_blrmax[zn] = fabs(
						 gamdotblr(gamma_max_fs, -1.0, 0.0, lineepp, lineepm,
							   contepp, contepm, lineeplin, conteplin,
							   nph_linelin, nph_contlin, etaph_lin,
							   nph_linegausp, nph_linegausm, weights-1,
							   etaph_gausp, etaph_gausm, nph_contgausp,
							   nph_contgausm, GAUSASZ));
		      
		    } else {

		      gammadot_blrmax[zn] = 0.0;
		    }

		    if ((L_DT != 0.0) && (LDT_fs != 0)) {
		      
		      gammadot_dtmax[zn] = fabs(gamdotdt(gamma_max_fs, etadtarr, dtep, nph_dt));
		      
		    } else {

		      gammadot_dtmax[zn] = 0.0;
		    }

		    if (timecount_fs == 0) {

		      gammadot_sscmax[zn] = 0.0;
		      
		    } else {
		      
		      gammadot_sscmax[zn] = fabs(gammasscdot(epsil, n_phgg[zn], gamma_max_fs, FREQ_GRID));
			
		    }
		  }
		  
		  totgammadotmax[zn] = gammadot_synmax[zn] + gammadot_sscmax[zn] + gammadot_ecdmax[zn] 
		    + gammadot_blrmax[zn] + gammadot_dtmax[zn];
		  
		  if ((accf == FALSE) || ((accf == TRUE) && (zn != (f - 1)))) {
		    dtcool_zn_fs = gam_max[zn] / totgammadotmax[zn];
		    
		    if (fmod(write_count, 10) == 0) {
		      fprintf(fp_log,
			      "totgammadotmax[%d] = %e, gam_max[%d] = %e, syndot = %e, sscdot = %e, ecddot = %e, blrdot = %e, dtdot = %e\n",
			      zn, totgammadotmax[zn], zn, gam_max[zn],
			      gammadot_synmax[zn], gammadot_sscmax[zn],
			      gammadot_ecdmax[zn], gammadot_blrmax[zn], gammadot_dtmax[zn]);
		    }
		    
		  } else {
		    dtcool_zn_fs = gamma_max_fs / totgammadotmax[zn];
		    
		    if (fmod(write_count, 10) == 0) {
		      fprintf(fp_log,
			      "totgammadotmax[%d] = %e, gamma_max_fs = %e, syndot = %e, sscdot = %e, ecddot = %e, blrdot = %e, dtdot = %e\n",
			      zn, totgammadotmax[zn], gamma_max_fs,
			      gammadot_synmax[zn], gammadot_sscmax[zn],
			      gammadot_ecdmax[zn], gammadot_blrmax[zn], gammadot_dtmax[zn]);
		    }
		  }
		  
		  dtcool_fs = min(dtcool_fs, dtcool_zn_fs);
		  
		}
		
		dtcool_rs = 1.0e+10;
		for (zn = (Nrs - 1); zn >= r; zn--) {

		  if ((accr == FALSE) || ((accr == TRUE) && (zn != r))) {
		    
		    if (gam_max[zn] == 0.0) {

		      gam_max[zn] = gamma_max_rs;
		    }
		    
		    gammadot_synmax[zn] = fabs(gammasyndot(gam_max[zn], B_rs));
		    
		    gammadot_sscmax[zn] = fabs(gammasscdot(epsil, n_phgg[zn], gam_max[zn], FREQ_GRID));
		      
		    if ((Ldisk != 0.0) && (LDISK_rs != 0)) {
		      
		      gammadot_ecdmax[zn] = fabs(gammaecdiskdot(gam_max[zn], Gamma_sh, z_c, thetar, raddist));
		      
		    } else {

		      gammadot_ecdmax[zn] = 0.0;
		    }
		    
		    if ((L_BLR != 0.0) && (LBLR_rs != 0)) {
		      
		      gammadot_blrmax[zn] = fabs(
						 gamdotblr(gam_max[zn], -1.0, 0.0, lineepp, lineepm,
							   contepp, contepm, lineeplin, conteplin,
							   nph_linelin, nph_contlin, etaph_lin,
							   nph_linegausp, nph_linegausm, weights-1,
							   etaph_gausp, etaph_gausm, nph_contgausp,
							   nph_contgausm, GAUSASZ));
		      
		    } else {

		      gammadot_blrmax[zn] = 0.0;
		    }
		    
		    if ((L_DT != 0.0) && (LDT_rs != 0)) {
		      
		      gammadot_dtmax[zn] = fabs(gamdotdt(gam_max[zn], etadtarr, dtep, nph_dt));
		      
		    } else {

		      gammadot_dtmax[zn] = 0.0;
		    }
		    
		  } else {
		    
		    gammadot_synmax[zn] = gammadotsynmax_rs;
		    
		    if ((Ldisk != 0.0) && (LDISK_rs != 0)) {
		      
		      gammadot_ecdmax[zn] = fabs(gammaecdiskdot(gamma_max_rs, Gamma_sh, z_c, thetar, raddist));
		      
		    } else {

		      gammadot_ecdmax[zn] = 0.0;
		    }
		    
		    if ((L_BLR != 0.0) && (LBLR_rs != 0)) {
		      
		      gammadot_blrmax[zn] = fabs(
						 gamdotblr(gamma_max_rs, -1.0, 0.0, lineepp, lineepm,
							   contepp, contepm, lineeplin, conteplin,
							   nph_linelin, nph_contlin, etaph_lin,
							   nph_linegausp, nph_linegausm, weights-1,
							   etaph_gausp, etaph_gausm, nph_contgausp,
							   nph_contgausm, GAUSASZ));
		      
		    } else {

		      gammadot_blrmax[zn] = 0.0;
		    }
		    
		    if ((L_DT != 0.0) && (LDT_rs != 0)) {
		      
		      gammadot_dtmax[zn] = fabs(gamdotdt(gamma_max_rs, etadtarr, dtep, nph_dt));
		      
		    } else {

		      gammadot_dtmax[zn] = 0.0;
		    }
		    
		    if (timecount_rs == 0) {

		      gammadot_sscmax[zn] = 0.0;

		    } else {

		      gammadot_sscmax[zn] = fabs(gammasscdot(epsil, n_phgg[zn], gamma_max_rs, FREQ_GRID));
			
		    }
		    
		  }
		  
		  totgammadotmax[zn] = gammadot_synmax[zn] + gammadot_sscmax[zn] + gammadot_ecdmax[zn] 
		    + gammadot_blrmax[zn] + gammadot_dtmax[zn];
		  
		  if ((accr == FALSE) || ((accr == TRUE) && (zn != r))) {
		    dtcool_zn_rs = gam_max[zn] / totgammadotmax[zn];
		    
		    if (fmod(write_count, 10) == 0) {
		      fprintf(fp_log,
			      "totgammadotmax[%d] = %e, gam_max[%d] = %e, syndot = %e, sscdot = %e, ecddot = %e, blrdot = %e, dtdot = %e\n",
			      zn, totgammadotmax[zn], zn, gam_max[zn],
			      gammadot_synmax[zn], gammadot_sscmax[zn],
			      gammadot_ecdmax[zn], gammadot_blrmax[zn], gammadot_dtmax[zn]);
		    }

		  } else {
		    dtcool_zn_rs = gamma_max_rs / totgammadotmax[zn];
		    
		    if (fmod(write_count, 10) == 0) {
		      fprintf(fp_log,
			      "totgammadotmax[%d] = %e, gamma_max_rs = %e, syndot = %e, sscdot = %e, ecddot = %e, blrdot = %e, dtdot = %e\n",
			      zn, totgammadotmax[zn], gamma_max_rs,
			      gammadot_synmax[zn], gammadot_sscmax[zn],
			      gammadot_ecdmax[zn], gammadot_blrmax[zn], gammadot_dtmax[zn]);
		    }
		  }
		  
		  dtcool_rs = min(dtcool_rs, dtcool_zn_rs);
		  
		}

		if (accf == TRUE) {
		  dt_fs = min(tph_esc_shfs, min(tcr_fs_zn, min(tph_esc_fs, min(dtcool_fs, min(dtinj_fs, te_esc)))));

		  if (fmod(write_count, 10) == 0) {

				fprintf(fp_log,
						"dtcool_fs = %e, tph_esc_fs = %e, dtinj_fs = %e\n",
						dtcool_fs, tph_esc_fs, dtinj_fs);
				fprintf(fp_log, "te_esc = %e, tcr_fs_zn = %e\n", te_esc,
						tcr_fs_zn);
				fprintf(fp_log, "dt_fs = %e\n", dt_fs);

			}

		} else {
			dt_fs = min(tph_esc_fs, min(dtcool_fs, te_esc));

			if (fmod(write_count, 10) == 0) {

				fprintf(fp_log,
						"dtcool_fs = %e, tph_esc_fs = %e, te_esc = %e, dt_fs = %e\n",
						dtcool_fs, tph_esc_fs, te_esc, dt_fs);

			}

		}

		if (accr == TRUE) {
		  dt_rs = min(tph_esc_shrs, min(tcr_rs_zn, min(tph_esc_rs, min(dtcool_rs, min(dtinj_rs, te_esc)))));

		  if (fmod(write_count, 10) == 0) {

				fprintf(fp_log,
						"dtcool_rs = %e, tph_esc_rs = %e, dtinj_rs = %e\n",
						dtcool_rs, tph_esc_rs, dtinj_rs);
				fprintf(fp_log, "te_esc = %e, tcr_rs_zn = %e\n", te_esc,
						tcr_rs_zn);
				fprintf(fp_log, "dt_rs = %e\n", dt_rs);

			}

		} else {
			dt_rs = min(tph_esc_rs, min(dtcool_rs, te_esc));

			if (fmod(write_count, 10) == 0) {

				fprintf(fp_log,
						"dtcool_rs = %e, tph_esc_rs = %e, te_esc = %e, dt_rs = %e\n",
						dtcool_rs, tph_esc_rs, te_esc, dt_rs);

			}

		}

		dt = min(dt_fs, dt_rs);
		if ((dt == tph_esc_fs) || (dt == tph_esc_rs))
			dt = 0.25 * dt;
		else
			dt = 0.5 * dt;

                if (dt < 10.0) {

		  printf("\nTime step, dt = %e, too small at time counter = %d\n", dt, time_count);
                  printf("\ndtcool_fs = %e, tph_esc_fs = %e, dtinj_fs = %e\n",
                         dtcool_fs, tph_esc_fs, dtinj_fs);
		  printf("te_esc = %e, tcr_fs_zn = %e, tcr_rs_zn = %e\n",
                         te_esc, tcr_fs_zn, tcr_rs_zn);

                  printf("\ndtcool_rs = %e, tph_esc_rs = %e, dtinj_rs = %e\n",
                         dtcool_rs, tph_esc_rs, dtinj_rs);
                  printf("dt_fs = %e, dt_rs = %e\n", dt_fs, dt_rs);

                  exitSimulation (start, logfile, -1, sOutputDirectory, nfnfile);

                }

		simtime += dt;

		if (f != Ntot) {
			time_fs = ((double) (f - Nrs)) * tcr_fs_zn;

		} else
			time_fs = ((double) (Nfs)) * tcr_fs_zn;

		if (r != 0) {
			time_rs = ((double) (Nrs - r)) * tcr_rs_zn;

		} else
			time_rs = ((double) (Nrs)) * tcr_rs_zn;

		if (fmod(write_count, 10) == 0) {

			fprintf(fp_log, "dt = %e, time = %e, time_fs = %e, time_rs = %e\n",
					dt, simtime, time_fs, time_rs);

			fclose(fp_log);

		}


		// Calculate emitting volume related parameters as long as the
		// shocks are in the system.
		//
		if (accf == TRUE) {

		  dsh_fs = (beta_fs * c * simtime) - (((double) (f - Nrs - 1)) * znht_fs);
		  Volcylsh_fs = pi * SQR(R_z) * dsh_fs;
		  volshfs_param = 5.355265e-47 * Volcylsh_fs;
		  Voldlsh_fs = (SQR(R_z) * dsh_fs * 0.25) / SQR(dl);
		  voldlshfs_param = 5.355265e-47 * Voldlsh_fs * SQR((D * D)) * 1.0e23 / (1.0 + Z);

		  dshfs = beta_fs * c * simtime;
		  Voldlshfs = (SQR(R_z) * dshfs * 0.25) / SQR(dl);
		  voldlshfsparam = 5.355265e-47 * Voldlshfs * SQR((D * D)) * 1.0e23 / (1.0 + Z);

		  if (fmod(write_count, 10) == 0) {

		    fp_log = fopen(logfile, "a");
                    if (!fp_log) {

		      printf("\nCould not open the output file log. \n");
                      return -1;
                    }

		    fprintf(fp_log, "f = %d, dsh_fs = %e, Volcylsh_fs = %e, Voldlsh_fs = %e\n",
			    f, dsh_fs, Volcylsh_fs, Voldlsh_fs);
		    fprintf(fp_log, "dshfs = %e, Voldlshfs = %e, voldlshfsparam = %e\n",
			    dshfs, Voldlshfs, voldlshfsparam);

		  }

		  tph_esc_shfs = esctimecalc(R_z, dsh_fs);
		  lmean_shfs = tph_esc_shfs * c;
		  lmeanshfs_param = 6.234375e-26 * lmean_shfs; // Value = 3*SIGMAT/32

		  P_up_shfs = (0.5 * R_z) / (R_z + dsh_fs);
		  pupshfs_param = P_up_shfs / tph_esc_shfs;
		  P_side_shfs = dsh_fs / (R_z + dsh_fs);
		  psideshfs_param = P_side_shfs / tph_esc_shfs;

		  if (fmod(write_count, 10) == 0) {

		    fprintf(fp_log, "f = %d, tph_esc_shfs = %e, P_up_shfs = %e, P_side_shfs = %e, lmean_shfs = %e\n",
			    f, tph_esc_shfs, P_up_shfs, P_side_shfs, lmean_shfs);

		    fclose(fp_log);

		  }

		}


		if (accr == TRUE) {

		  dsh_rs = (beta_rs * c * simtime) - (((double) (Nrs - r - 1)) * znht_rs);
		  Volcylsh_rs = pi * SQR(R_z) * dsh_rs;
		  volshrs_param = 5.355265e-47 * Volcylsh_rs;
		  Voldlsh_rs = (SQR(R_z) * dsh_rs * 0.25) / SQR(dl);
		  voldlshrs_param = 5.355265e-47 * Voldlsh_rs * SQR((D * D)) * 1.0e23 / (1.0 + Z);

		  if (mu_prime < 0.0) {

		    dshrs = beta_rs * c * simtime;
		    Voldlshrs = (SQR(R_z) * dshrs * 0.25) / SQR(dl);
		    voldlshrsparam = 5.355265e-47 * Voldlshrs * SQR((D * D)) * 1.0e23 / (1.0 + Z);

		  }

		  if (fmod(write_count, 10) == 0) {

                    fp_log = fopen(logfile, "a");
                    if (!fp_log) {

                      printf("\nCould not open the output file log. \n");
                      return -1;
                    }

		    fprintf(fp_log, "r = %d, dsh_rs = %e, Volcylsh_rs = %e, Voldlsh_rs = %e\n",
			    r, dsh_rs, Volcylsh_rs, Voldlsh_rs);

		    if (mu_prime < 0.0) {

		      fprintf(fp_log, "dshrs = %e, Voldlshrs = %e, voldlshrsparam = %e\n",
			      dshrs, Voldlshrs, voldlshrsparam);

		    }

		  }

		  tph_esc_shrs = esctimecalc(R_z, dsh_rs);
		  lmean_shrs = tph_esc_shrs * c;
		  lmeanshrs_param = 6.234375e-26 * lmean_shrs;

		  P_up_shrs = (0.5 * R_z) / (R_z + dsh_rs);
		  pupshrs_param = P_up_shrs / tph_esc_shrs;
		  P_side_shrs = dsh_rs / (R_z + dsh_rs);
		  psideshrs_param = P_side_shrs / tph_esc_shrs;

		  if (fmod(write_count, 10) == 0) {

		    fprintf(fp_log, "r = %d, tph_esc_shrs = %e, P_up_shrs = %e, P_side_shrs = %e, lmean_shrs = %e\n",
			    r, tph_esc_shrs, P_up_shrs, P_side_shrs, lmean_shrs);

		    fclose(fp_log);

		  }

		}


		// Calculate the initial electron density in both regions.
		//
		if (timecount_fs == 0) {

		  if (time_count == 0) {

		    fp_fsinit = fopen(fsinitfile, "w");
		    if (!fp_fsinit) {

		      printf("\nCould not open the initial density fs file.\n");
		      return -1;
		    }

		  }

		  for (j = 0; j < EGRID; j++) {

		    n_e0_fs[j] = Q_inj_fs[j] * dt;
		    n_e_p0_fs[j] = 2.0 * n_e0_fs[j];

		    if (time_count == 0) {
		      fprintf(fp_fsinit, "%e %e\n", n_e0_fs[j], n_e_p0_fs[j]);
		    }

		  }

		  if (time_count == 0) {
		    fclose(fp_fsinit);
		  }

		}

		if (timecount_rs == 0) {

		  if (time_count == 0) {

		    fp_rsinit = fopen(rsinitfile, "w");
		    if (!fp_rsinit) {

		      printf("\nCould not open the initial density rs file.\n");
		      return -1;
		    }

		  }

		  for (j = 0; j < EGRID; j++) {

		    n_e0_rs[j] = Q_inj_rs[j] * dt;
		    n_e_p0_rs[j] = 2.0 * n_e0_rs[j];

		    if (time_count == 0) {
		      fprintf(fp_rsinit, "%e %e\n", n_e0_rs[j], n_e_p0_rs[j]);
		    }

		  }

		  if (time_count == 0) {
		    fclose(fp_rsinit);
		  }

		}

		// Calculate the total electron number density
		//
		if (time_count == 0) {
		  e_numden_fs = ne_numtot(gamma_e, n_e0_fs, EGRID);
		  e_numden_rs = ne_numtot(gamma_e, n_e0_rs, EGRID);
		  //printf("e_numden_fs = %e, e_numden_rs = %e, time count = %d\n", e_numden_fs, e_numden_rs, time_count);
		  
		}

		// Calculate the characteristic frequency of gamma_min to normalize the
		// synchrotron power obtained from the approximate expression at frequencies
		// below nu_gamma1, when gamma_b is changing with time.
		//
		for (zn = Nrs; zn < f; zn++) {
		  if ((accf == FALSE) || ((accf == TRUE) && (zn != (f - 1)))) {
		    if (gam_min[zn] == 0.0) {
		      gam_min[zn] = gamma_b_fs;
		    }
		    
		    // Corresponding to a completely randomly oriented field:
		    //
		    nu_gmin_rand[zn] = NU_C * B_fs * SQR(gam_min[zn]);

		    // Corresponding to the ordered component of the field:
		    //
		    if ((Bflag != 0.0) && (Bflag != 4.0) && (Bflag != 5.0)) {
		      
		      nu_gmin[zn] = NU_C * B_fs * sinezi * SQR(gam_min[zn]);
		      nu_gminfeed[zn] = NU_C * B_fs * sinezifeed * SQR(gam_min[zn]);
		      
		    } else if ((Bflag == 4.0) || (Bflag == 5.0)) {
		      
		      for (j = 0; j < QUANTIZE; j++) {
			
			nu_gminord[zn][j] = NU_C * B_fs * sineord[j] * SQR(gam_min[zn]);
		      }

		      nu_gminfeed[zn] = NU_C * B_fs * sinezifeed * SQR(gam_min[zn]);
		      
		    } else {
		      
		      nu_gmin[zn] = nu_gmin_rand[zn];
		      
		    }
		    
		  } else {
		    
		    nu_gmin_rand[zn] = nu_gamma1_fsrand;
		    
		    if ((Bflag != 0.0) && (Bflag != 4.0) && (Bflag != 5.0)) {
		      
		      nu_gmin[zn] = nu_gamma1_fs;
		      nu_gminfeed[zn] = nu_gamma1_fsfeed;
		      
		    } else if ((Bflag == 4.0) || (Bflag == 5.0)) {
		      
		      for (j = 0; j < QUANTIZE; j++) {
			
			nu_gminord[zn][j] = nu_gam1fs_ord[j];
		      }

		      nu_gminfeed[zn] = nu_gamma1_fsfeed;
		      
		    } else {
		      
		      nu_gmin[zn] = nu_gmin_rand[zn];
		      
		    }
		    
		  }
		  
		  if (fmod(write_count, 10) == 0) {
		    fp_log = fopen(logfile, "a");
		    if (!fp_log) {
		      printf("\nCould not open the output file log. \n");
		      return -1;
		    }
		    
		    if ((Bflag != 4.0) && (Bflag != 5.0)) {
		      
		      fprintf(fp_log, "nu_gmin_rand[%d] = %e, nu_gmin[%d] = %e, gam_min[%d] = %e, accf = %d\n", 
			      zn, nu_gmin_rand[zn], zn, nu_gmin[zn], zn, gam_min[zn], accf);
		      
		    } else {
		      
		      for (j = 0; j < QUANTIZE; j++) {
			
			fprintf(fp_log, "nu_gmin_rand[%d] = %e, nu_gminorder[%d][%d] = %e, gam_min[%d] = %e, accf = %d\n", 
				zn, nu_gmin_rand[zn], zn, j, nu_gminord[zn][j], zn, gam_min[zn], accf);
			
		      }
		      
		    } 

		    fprintf(fp_log, "nu_gminfeed[%d] = %e\n", zn, nu_gminfeed[zn]);
		    
		    fclose(fp_log);
		  }
		}
		
		for (zn = (Nrs - 1); zn >= r; zn--) {
		  if ((accr == FALSE) || ((accr == TRUE) && (zn != r))) {
		    if (gam_min[zn] == 0.0) {
		      gam_min[zn] = gamma_b_rs;
		    }
		    
		    // Corresponding to a completely randomly oriented field:
		    //
		    nu_gmin_rand[zn] = NU_C * B_rs * SQR(gam_min[zn]);

		    // Corresponding to the ordered component of the field:
		    //		    
		    if ((Bflag != 0.0) && (Bflag != 4.0) && (Bflag != 5.0)) {

		      nu_gmin[zn] = NU_C * B_rs * sinezi * SQR(gam_min[zn]);
		      nu_gminfeed[zn] = NU_C * B_rs * sinezifeed * SQR(gam_min[zn]);
		      
		    } else if ((Bflag == 4.0) || (Bflag == 5.0)) {

		      for (j = 0; j < QUANTIZE; j++) {
			
			nu_gminord[zn][j] = NU_C * B_rs * sineord[j] * SQR(gam_min[zn]);
		      }

		      nu_gminfeed[zn] = NU_C * B_rs * sinezifeed * SQR(gam_min[zn]);
		      
		    } else {
			
			nu_gmin[zn] = nu_gmin_rand[zn];

		    }
		    
		  } else {
		    
		    nu_gmin_rand[zn] = nu_gamma1_rsrand;

		    if ((Bflag != 0.0) && (Bflag != 4.0) && (Bflag != 5.0)) {

		      nu_gmin[zn] = nu_gamma1_rs;
		      nu_gminfeed[zn] = nu_gamma1_rsfeed;
		      
		    } else if ((Bflag == 4.0) || (Bflag == 5.0)) {

		      for (j = 0; j < QUANTIZE; j++) {
			
			nu_gminord[zn][j] = nu_gam1rs_ord[j];
		      }

		      nu_gminfeed[zn] = nu_gamma1_rsfeed;
		      
		    } else {

		      nu_gmin[zn] = nu_gmin_rand[zn];

		    }

		  }
		  
		  if (fmod(write_count, 10) == 0) {
		    fp_log = fopen(logfile, "a");
		    if (!fp_log) {
		      printf("\nCould not open the output file log. \n");
		      return -1;
		    }
		    
		    if ((Bflag != 4.0) && (Bflag != 5.0)) {
		      
		      fprintf(fp_log, "nu_gmin_rand[%d] = %e, nu_gmin[%d] = %e, gam_min[%d] = %e, accr = %d\n", 
			      zn, nu_gmin_rand[zn], zn, nu_gmin[zn], zn, gam_min[zn], accr);
		      
		    } else {
		      
		      for (j = 0; j < QUANTIZE; j++) {
			
			fprintf(fp_log, "nu_gmin_rand[%d] = %e, nu_gminord[%d][%d] = %e, gam_min[%d] = %e, accr = %d\n", 
				zn, nu_gmin_rand[zn], zn, j, nu_gminord[zn][j], zn, gam_min[zn], accr);
			
		      }		      
		      
		    } 

		    fprintf(fp_log, "nu_gminfeed[%d] = %e\n", zn, nu_gminfeed[zn]);

		    fclose(fp_log);
		  }
		}

		// Calculate synchrotron spectrum including the synchrotron self
		// absorption process here
		//
		//#pragma omp parallel for private(k)
		for (zn = Nrs; zn < f; zn++) {

#pragma omp parallel for private(p)
		  for (k = 0; k < FREQ_GRID; k++) {
		    
		    if ((timecount_fs == 0) && (zn == (f - 1))) {

		      // Calculate synchrotron radiation for a randomly oriented field:
		      //
		      ndot_syn_ep[zn][k] = ndotsyn(syn_nu[k], nu_gmin_rand[zn], B_fs, 
                                                   gamma_e, n_e_p0_fs, EGRID);
		      
		      tau_ssa[zn][k] = lmean_shfs * alpha_ssa(syn_nu[k], nu_gmin_rand[zn], 
		                                              B_fs, gamma_e, n_e_p0_fs);
		      
		      // Calculate synchrotron radiation for the ordered component:
		      //
		      if ((Bflag != 0.0) && (Bflag != 4.0) && (Bflag != 5.0)) {
			
                        ndot_syn_uni[zn][k] = ndotsynuni(syn_nu[k], nu_gmin[zn], B_fs,
							 gamma_e, n_e_p0_fs, EGRID, sinezi);
			
                        tau_ssa_uni[zn][k] = lmean_shfs * alpha_ssauni(syn_nu[k], nu_gmin[zn], 
								       B_fs, gamma_e, n_e_p0_fs, sinezi);
			
                        ndot_syn_unifeed[zn][k] = ndotsynuni(syn_nu[k], nu_gminfeed[zn], B_fs,
							     gamma_e, n_e_p0_fs, EGRID, sinezifeed);
			
                        tau_ssa_unifeed[zn][k] = lmean_shfs * alpha_ssauni(syn_nu[k], nu_gminfeed[zn],
									   B_fs, gamma_e, n_e_p0_fs, sinezifeed);
			
		      } else if ((Bflag == 4.0) || (Bflag == 5.0)) {
			
			ndot_syn_uni[zn][k] = 0.0;
			tau_ssa_uni[zn][k] = 0.0;
			
			for (p = 0; p < QUANTIZE; p++) {
			  
			  ndot_syn_uni[zn][k] += ndotsynuni(syn_nu[k], nu_gminord[zn][p], B_fs,
							    gamma_e, n_e_p0_fs, EGRID, sineord[p]);
			  
			  ndot_syn_ord[zn][k][p] = ndotsynuni(syn_nu[k], nu_gminord[zn][p], B_fs,
							      gamma_e, n_e_p0_fs, EGRID, sineord[p]);
			  
			  tau_ssa_uni[zn][k] += (lmean_shfs * alpha_ssauni(syn_nu[k], nu_gminord[zn][p], 
									   B_fs, gamma_e, n_e_p0_fs, sineord[p]));
			  
			  tau_ssa_ord[zn][k][p] = lmean_shfs * alpha_ssauni(syn_nu[k], nu_gminord[zn][p],
									    B_fs, gamma_e, n_e_p0_fs, sineord[p]);
			  /*
			    printf("In FS, time_count = %d, ndot_syn_ord[%d][%d][%d] = %e, tau_ssa_ord[%d][%d][%d] = %e\n",
			    time_count, zn, k, p, ndot_syn_ord[zn][k][p], zn, k, p, tau_ssa_ord[zn][k][p]);
			  */
			}
			
                        ndot_syn_unifeed[zn][k] = ndotsynuni(syn_nu[k], nu_gminfeed[zn], B_fs,
							     gamma_e, n_e_p0_fs, EGRID, sinezifeed);
			
                        tau_ssa_unifeed[zn][k] = lmean_shfs * alpha_ssauni(syn_nu[k], nu_gminfeed[zn],
									   B_fs, gamma_e, n_e_p0_fs, sinezifeed);
			
			/*
			  printf("In FS, time_count = %d, ndot_syn_uni[%d][%d] = %e, tau_ssa_uni[%d][%d] = %e\n", 
			  time_count, zn, k, ndot_syn_uni[zn][k], zn, k, tau_ssa_uni[zn][k]);
			*/
			
		      } else {
			
			ndot_syn_uni[zn][k] = 0.0;
			tau_ssa_uni[zn][k] = 0.0;
			ndot_syn_unifeed[zn][k] = 0.0;
			tau_ssa_unifeed[zn][k] = 0.0;
			
		      }
		      
                    } else {
		      
		      // For randomly oriented field:
		      //
		      ndot_syn_ep[zn][k] = ndotsyn(syn_nu[k], nu_gmin_rand[zn], B_fs, 
						   gamma_e, n_e_p[zn], EGRID);
		      
		      if ((accf == TRUE) && (zn == (f - 1))) {
			
			tau_ssa[zn][k] = lmean_shfs * alpha_ssa(syn_nu[k], nu_gmin_rand[zn], 
								B_fs, gamma_e, n_e_p[zn]);
		      } else {
			
			tau_ssa[zn][k] = lmean_fs * alpha_ssa(syn_nu[k], nu_gmin_rand[zn], 
							      B_fs, gamma_e, n_e_p[zn]);
		      }
		      
		      // For the ordered component of the field:
		      //
		      if ((Bflag != 0.0) && (Bflag != 4.0) && (Bflag != 5.0)) {
			
                        ndot_syn_uni[zn][k] = ndotsynuni(syn_nu[k], nu_gmin[zn], B_fs,
							 gamma_e, n_e_p[zn], EGRID, sinezi);
			
			ndot_syn_unifeed[zn][k] = ndotsynuni(syn_nu[k], nu_gminfeed[zn], B_fs,
							     gamma_e, n_e_p[zn], EGRID, sinezifeed);
			
			if ((accf == TRUE) && (zn == (f - 1))) {
			  
			  tau_ssa_uni[zn][k] = lmean_shfs * alpha_ssauni(syn_nu[k], nu_gmin[zn], 
									 B_fs, gamma_e, n_e_p[zn], sinezi);
			  
                          tau_ssa_unifeed[zn][k] = lmean_shfs * alpha_ssauni(syn_nu[k], nu_gminfeed[zn],
									     B_fs, gamma_e, n_e_p[zn], sinezifeed);
			} else {
			  
			  tau_ssa_uni[zn][k] = lmean_fs * alpha_ssauni(syn_nu[k], nu_gmin[zn], 
								       B_fs, gamma_e, n_e_p[zn], sinezi);
			  
                          tau_ssa_unifeed[zn][k] = lmean_fs * alpha_ssauni(syn_nu[k], nu_gminfeed[zn],
									   B_fs, gamma_e, n_e_p[zn], sinezifeed);
			}
			
		      } else if ((Bflag == 4.0) || (Bflag == 5.0)) {
			
			ndot_syn_uni[zn][k] = 0.0;
			tau_ssa_uni[zn][k] = 0.0;
			
			for (p = 0; p < QUANTIZE; p++) {
			  
			  ndot_syn_uni[zn][k] += ndotsynuni(syn_nu[k], nu_gminord[zn][p], B_fs,
							    gamma_e, n_e_p[zn], EGRID, sineord[p]);
			  
			  ndot_syn_ord[zn][k][p] = ndotsynuni(syn_nu[k], nu_gminord[zn][p], B_fs,
							      gamma_e, n_e_p[zn], EGRID, sineord[p]);
			  
			  if ((accf == TRUE) && (zn == (f - 1))) {
			    
			    tau_ssa_uni[zn][k] += (lmean_shfs * alpha_ssauni(syn_nu[k], nu_gminord[zn][p], 
									     B_fs, gamma_e, n_e_p[zn], sineord[p]));
			    
			    tau_ssa_ord[zn][k][p] = lmean_shfs * alpha_ssauni(syn_nu[k], nu_gminord[zn][p],
									      B_fs, gamma_e, n_e_p[zn], sineord[p]);
			  } else {
			    
			    tau_ssa_uni[zn][k] += (lmean_fs * alpha_ssauni(syn_nu[k], nu_gminord[zn][p], B_fs, gamma_e,
									   n_e_p[zn], sineord[p]));
			    
			    tau_ssa_ord[zn][k][p] = lmean_fs * alpha_ssauni(syn_nu[k], nu_gminord[zn][p], B_fs, gamma_e,
									    n_e_p[zn], sineord[p]);
			  }
			  
			}
			
                        ndot_syn_unifeed[zn][k] = ndotsynuni(syn_nu[k], nu_gminfeed[zn], B_fs,
							     gamma_e, n_e_p[zn], EGRID, sinezifeed);
			
                        if ((accf == TRUE) && (zn == (f - 1))) {
			  
                          tau_ssa_unifeed[zn][k] = lmean_shfs * alpha_ssauni(syn_nu[k], nu_gminfeed[zn],
									     B_fs, gamma_e, n_e_p[zn], sinezifeed);
			} else {
			  
                          tau_ssa_unifeed[zn][k] = lmean_fs * alpha_ssauni(syn_nu[k], nu_gminfeed[zn],
									   B_fs, gamma_e, n_e_p[zn], sinezifeed);
                        }
			
			
			/*
			  printf("In FS, time_count = %d, ndot_syn_uni[%d][%d] = %e, tau_ssa_uni[%d][%d] = %e\n", 
			  time_count, zn, k, ndot_syn_uni[zn][k], zn, k, tau_ssa_uni[zn][k]);
			*/
			
		      } else {
			
			ndot_syn_uni[zn][k] = 0.0;
			tau_ssa_uni[zn][k] = 0.0;
			ndot_syn_unifeed[zn][k] = 0.0;
			tau_ssa_unifeed[zn][k] = 0.0;
			
		      }
		      
		    }

		    ndot_syn_tot[zn][k] = (b_ord * ndot_syn_uni[zn][k]) + 
		      ((1.0 - b_ord) * ndot_syn_ep[zn][k]);
		    
		    tau_ssa_tot[zn][k] = (b_ord * tau_ssa_uni[zn][k]) + 
		      ((1.0 - b_ord) * tau_ssa[zn][k]); 

		    if (Bflag != 0.0) {
		      
		      ndot_syn_unifeedtot[zn][k] = (b_ord * ndot_syn_unifeed[zn][k]) +
			((1.0 - b_ord) * ndot_syn_ep[zn][k]);
		      
		      tau_ssa_unifeedtot[zn][k] = (b_ord * tau_ssa_unifeed[zn][k]) +
			((1.0 - b_ord) * tau_ssa[zn][k]);

		    }
 
		    power_syn[zn][k] = 5.355265e-47 * syn_nu[k] * ndot_syn_ep[zn][k]; // Value = h^2/(m_e*c^2)
		    
		    if (tau_ssa[zn][k] > TAULIMIT) {
		      ndot_syn_ep[zn][k] *= (1.0 - exp(-tau_ssa[zn][k])) / tau_ssa[zn][k];
		      
		    }
		    
		    if (tau_ssa_tot[zn][k] > TAULIMIT) {
		      ndot_syn_tot[zn][k] *= (1.0 - exp(-tau_ssa_tot[zn][k])) / tau_ssa_tot[zn][k];
		      
		    }
		    
		    if (tau_ssa_uni[zn][k] > TAULIMIT) {
		      ndot_syn_uni[zn][k] *= (1.0 - exp(-tau_ssa_uni[zn][k])) / tau_ssa_uni[zn][k];
		      
		    }				
		    
		    if (Bflag != 0.0) {
		      
		      if (tau_ssa_unifeedtot[zn][k] > TAULIMIT) {
			ndot_syn_unifeedtot[zn][k] *= (1.0 - exp(-tau_ssa_unifeedtot[zn][k])) / tau_ssa_unifeedtot[zn][k];
			
		      }
		      
		    }
		    
		    if ((Bflag == 4.0) || (Bflag == 5.0)) {
		      
		      for (p = 0; p < QUANTIZE; p++) {
			
			if (tau_ssa_ord[zn][k][p] > TAULIMIT) {
			  ndot_syn_ord[zn][k][p] *= (1.0 - exp(-tau_ssa_ord[zn][k][p])) / tau_ssa_ord[zn][k][p];
			  
			}
			
		      }
		      
		    }
				
		  } // END PARALLEL FOR SECTION

		} // END PARALLEL FOR SECTION on zn

                if (files) {

		  for (zn = Nrs; zn < f; zn++) {
		    for (k = 0; k < FREQ_GRID; k++) {
		      
		      if (zn == (f - 1)) {
                        fprintf(fp_fssa, "%e %e %e %e %e %e %e %e %e\n", syn_nu[k], epsil[k],
                                ndot_syn_ep[zn][k], max(1.0e-40, power_syn[zn][k]),
                                tau_ssa[zn][k], ndot_syn_uni[zn][k], tau_ssa_uni[zn][k], 
				ndot_syn_tot[zn][k], tau_ssa_tot[zn][k]);
			
		      }
		      
		      fprintf(fp_fsynnfn, "%e %e %e %e %d\n", syn_nu[k], ndot_syn_ep[zn][k], 
			      ndot_syn_uni[zn][k], ndot_syn_tot[zn][k], zn);
		      
		    }
		  }
		}

	//#pragma omp parallel for private(k)
	        for (zn = (Nrs - 1); zn >= r; zn--) {
		  
#pragma omp parallel for private(p)
		  for (k = 0; k < FREQ_GRID; k++) {
		    
                    if ((timecount_rs == 0) && (zn == r)) {
		      
                      // Calculate synchrotron radiation for a randomly oriented field:                        
		      // 
		      ndot_syn_ep[zn][k] = ndotsyn(syn_nu[k], nu_gmin_rand[zn], B_rs,
						   gamma_e, n_e_p0_rs, EGRID);

		      tau_ssa[zn][k] = lmean_shrs * alpha_ssa(syn_nu[k], nu_gmin_rand[zn], 
							      B_rs, gamma_e, n_e_p0_rs);

                      // Calculate synchrotron radiation for the ordered component:                    
		      //         
		      if ((Bflag != 0.0) && (Bflag != 4.0) && (Bflag != 5.0)) {

                        ndot_syn_uni[zn][k] = ndotsynuni(syn_nu[k], nu_gmin[zn], B_rs,
							 gamma_e, n_e_p0_rs, EGRID, sinezi);

                        tau_ssa_uni[zn][k] = lmean_shrs * alpha_ssauni(syn_nu[k], nu_gmin[zn], 
								       B_rs, gamma_e, n_e_p0_rs, sinezi);

                        ndot_syn_unifeed[zn][k] = ndotsynuni(syn_nu[k], nu_gminfeed[zn], B_rs,
							     gamma_e, n_e_p0_rs, EGRID, sinezifeed);

                        tau_ssa_unifeed[zn][k] = lmean_shrs * alpha_ssauni(syn_nu[k], nu_gminfeed[zn],
									   B_rs, gamma_e, n_e_p0_rs, sinezifeed);

		      } else if ((Bflag == 4.0) || (Bflag == 5.0)) {

			ndot_syn_uni[zn][k] = 0.0;
			tau_ssa_uni[zn][k] = 0.0;
			
			for (p = 0; p < QUANTIZE; p++) {
			  
			  ndot_syn_uni[zn][k] += ndotsynuni(syn_nu[k], nu_gminord[zn][p], B_rs,
							    gamma_e, n_e_p0_rs, EGRID, sineord[p]);

			  ndot_syn_ord[zn][k][p] = ndotsynuni(syn_nu[k], nu_gminord[zn][p], B_rs,
							      gamma_e, n_e_p0_rs, EGRID, sineord[p]);
			  
			  tau_ssa_uni[zn][k] += (lmean_shrs * alpha_ssauni(syn_nu[k], nu_gminord[zn][p], 
									   B_rs, gamma_e, n_e_p0_rs, sineord[p]));

			  tau_ssa_ord[zn][k][p] = lmean_shrs * alpha_ssauni(syn_nu[k], nu_gminord[zn][p],
									    B_rs, gamma_e, n_e_p0_rs, sineord[p]);
			  /*
			  printf("In RS, ndot_syn_ord[%d][%d][%d] = %e, tau_ssa_ord[%d][%d][%d] = %e\n",
                                 zn, k, p, ndot_syn_ord[zn][k][p], zn, k, p, tau_ssa_ord[zn][k][p]);
			  */
			}

                        ndot_syn_unifeed[zn][k] = ndotsynuni(syn_nu[k], nu_gminfeed[zn], B_rs,
							     gamma_e, n_e_p0_rs, EGRID, sinezifeed);

                        tau_ssa_unifeed[zn][k] = lmean_shrs * alpha_ssauni(syn_nu[k], nu_gminfeed[zn],
									   B_rs, gamma_e, n_e_p0_rs, sinezifeed);

			/*
			printf("In RS, time_count = %d, ndot_syn_uni[%d][%d] = %e, tau_ssa_uni[%d][%d] = %e\n", 
			       time_count, zn, k, ndot_syn_uni[zn][k], zn, k, tau_ssa_uni[zn][k]);
			*/
			
		      } else {

			ndot_syn_uni[zn][k] = 0.0;
			tau_ssa_uni[zn][k] = 0.0;
			ndot_syn_unifeed[zn][k] = 0.0;
			tau_ssa_unifeed[zn][k] = 0.0;
			
		      }

                    } else {

		      // For randomly oriented field:
		      ndot_syn_ep[zn][k] = ndotsyn(syn_nu[k], nu_gmin_rand[zn], B_rs,
						   gamma_e, n_e_p[zn], EGRID);

		      if ((accr == TRUE) && (zn == r)) {
		      
			tau_ssa[zn][k] = lmean_shrs * alpha_ssa(syn_nu[k], nu_gmin_rand[zn],
								B_rs, gamma_e, n_e_p[zn]);
		      } else {

			tau_ssa[zn][k] = lmean_rs * alpha_ssa(syn_nu[k], nu_gmin_rand[zn],
							      B_rs, gamma_e, n_e_p[zn]);
		      }

                      // For the ordered component of the field:                                                         
		      //                    
		      if ((Bflag != 0.0) && (Bflag != 4.0) && (Bflag != 5.0)) {
			
                        ndot_syn_uni[zn][k] = ndotsynuni(syn_nu[k], nu_gmin[zn], B_rs,
							 gamma_e, n_e_p[zn], EGRID, sinezi);

			ndot_syn_unifeed[zn][k] = ndotsynuni(syn_nu[k], nu_gminfeed[zn], B_rs,
							     gamma_e, n_e_p[zn], EGRID, sinezifeed);
			
			if ((accr == TRUE) && (zn == r)) {
			  
			  tau_ssa_uni[zn][k] = lmean_shrs * alpha_ssauni(syn_nu[k], nu_gmin[zn], 
									 B_rs, gamma_e, n_e_p[zn], sinezi);
			  
			  tau_ssa_unifeed[zn][k] = lmean_shrs * alpha_ssauni(syn_nu[k], nu_gminfeed[zn],
									     B_rs, gamma_e, n_e_p[zn], sinezifeed);
			} else {
			  
			  tau_ssa_uni[zn][k] = lmean_rs * alpha_ssauni(syn_nu[k], nu_gmin[zn], 
								       B_rs, gamma_e, n_e_p[zn], sinezi);

			  tau_ssa_unifeed[zn][k] = lmean_rs * alpha_ssauni(syn_nu[k], nu_gminfeed[zn],
									   B_rs, gamma_e, n_e_p[zn], sinezifeed);
			}

		      } else if ((Bflag == 4.0) || (Bflag == 5.0)) {

			ndot_syn_uni[zn][k] = 0.0;
			tau_ssa_uni[zn][k] = 0.0;

			for (p = 0; p < QUANTIZE; p++) {

			  ndot_syn_uni[zn][k] += ndotsynuni(syn_nu[k], nu_gminord[zn][p], B_rs,
							    gamma_e, n_e_p[zn], EGRID, sineord[p]);

			  ndot_syn_ord[zn][k][p] = ndotsynuni(syn_nu[k], nu_gminord[zn][p], B_rs,
							      gamma_e, n_e_p[zn], EGRID, sineord[p]);

			  if ((accr == TRUE) && (zn == r)) {
			    
			    tau_ssa_uni[zn][k] += (lmean_shrs * alpha_ssauni(syn_nu[k], nu_gminord[zn][p], 
									     B_rs, gamma_e, n_e_p[zn], sineord[p]));

			    tau_ssa_ord[zn][k][p] = lmean_shrs * alpha_ssauni(syn_nu[k], nu_gminord[zn][p],
									      B_rs, gamma_e, n_e_p[zn], sineord[p]);
			  } else {

			    tau_ssa_uni[zn][k] += (lmean_rs * alpha_ssauni(syn_nu[k], nu_gminord[zn][p], 
									   B_rs, gamma_e, n_e_p[zn], sineord[p]));

			    tau_ssa_ord[zn][k][p] = lmean_rs * alpha_ssauni(syn_nu[k], nu_gminord[zn][p],
									    B_rs, gamma_e, n_e_p[zn], sineord[p]);
			  }
			  
			}

                        ndot_syn_unifeed[zn][k] = ndotsynuni(syn_nu[k], nu_gminfeed[zn], B_rs,
							     gamma_e, n_e_p[zn], EGRID, sinezifeed);

                        if ((accr == TRUE) && (zn == r)) {

                          tau_ssa_unifeed[zn][k] = lmean_shrs * alpha_ssauni(syn_nu[k], nu_gminfeed[zn],
									     B_rs, gamma_e, n_e_p[zn], sinezifeed);
                        } else {
			  
                          tau_ssa_unifeed[zn][k] = lmean_rs * alpha_ssauni(syn_nu[k], nu_gminfeed[zn],
									   B_rs, gamma_e, n_e_p[zn], sinezifeed);
                        }

			/*
			printf("In RS, time_count = %d, ndot_syn_uni[%d][%d] = %e, tau_ssa_uni[%d][%d] = %e\n", 
			       time_count, zn, k, ndot_syn_uni[zn][k], zn, k, tau_ssa_uni[zn][k]);
			*/
			
		      } else {

			ndot_syn_uni[zn][k] = 0.0;
			tau_ssa_uni[zn][k] = 0.0;
			ndot_syn_unifeed[zn][k] = 0.0;
			tau_ssa_unifeed[zn][k] = 0.0;
			
		      }

                    }

		    ndot_syn_tot[zn][k] = (b_ord * ndot_syn_uni[zn][k]) + 
		      ((1.0 - b_ord) * ndot_syn_ep[zn][k]);

		    tau_ssa_tot[zn][k] = (b_ord * tau_ssa_uni[zn][k]) + 
		      ((1.0 - b_ord) * tau_ssa[zn][k]);

		    if (Bflag != 0.0) {
		      
		      ndot_syn_unifeedtot[zn][k] = (b_ord * ndot_syn_unifeed[zn][k]) +
			((1.0 - b_ord) * ndot_syn_ep[zn][k]);

                    tau_ssa_unifeedtot[zn][k] = (b_ord * tau_ssa_unifeed[zn][k]) +
                      ((1.0 - b_ord) * tau_ssa[zn][k]);

		    }
		    
                    power_syn[zn][k] = 5.355265e-47 * syn_nu[k] * ndot_syn_ep[zn][k];

                    if (tau_ssa[zn][k] > TAULIMIT) {
		      ndot_syn_ep[zn][k] *= (1.0 - exp(-tau_ssa[zn][k])) / tau_ssa[zn][k];

                    }

                    if (tau_ssa_tot[zn][k] > TAULIMIT) {
                      ndot_syn_tot[zn][k] *= (1.0 - exp(-tau_ssa_tot[zn][k])) / tau_ssa_tot[zn][k];

                    }

		    if (tau_ssa_uni[zn][k] > TAULIMIT) {
		      ndot_syn_uni[zn][k] *= (1.0 - exp(-tau_ssa_uni[zn][k])) / tau_ssa_uni[zn][k];

		    }

		    if (Bflag != 0.0) {

		      if (tau_ssa_unifeedtot[zn][k] > TAULIMIT) {
			ndot_syn_unifeedtot[zn][k] *= (1.0 - exp(-tau_ssa_unifeedtot[zn][k])) / tau_ssa_unifeedtot[zn][k];
			
		      }

		    }

                    if ((Bflag == 4.0) || (Bflag == 5.0)) {

                      for (p = 0; p < QUANTIZE; p++) {

			if (tau_ssa_ord[zn][k][p] > TAULIMIT) {
                          ndot_syn_ord[zn][k][p] *= (1.0 - exp(-tau_ssa_ord[zn][k][p])) / tau_ssa_ord[zn][k][p];

			}

                      }

                    }

		  } // END PARALLEL FOR SECTION

		} // END PARALLEL FOR SECTION on zn

		if (files) {
		  for (zn = (Nrs - 1); zn >= r; zn--) {
		    for (k = 0; k < FREQ_GRID; k++) {

		      fprintf(fp_rsynnfn, "%e %e %e %e %e %e %e %d\n", syn_nu[k], 
			      ndot_syn_ep[zn][k], ndot_syn_uni[zn][k], ndot_syn_tot[zn][k], 
			      tau_ssa[zn][k], tau_ssa_uni[zn][k], tau_ssa_tot[zn][k], zn);

		    }
		  }
		}

		if (files) {
		  fclose(fp_fssa);
		  fclose(fp_fsynnfn);
		  fclose(fp_rsynnfn);
		}


		/*
		 sprintf(filename2_fs, "%s%d%s", fselossfile, time_count, ".dat");
		 fp_fseloss = fopen(filename2_fs, "a");
		 if(!fp_fseloss) {

		   printf("\nCould not open the output file fseloss. \n");
		   return -1;
		 }

		 sprintf(filename2_rs, "%s%d%s", rselossfile, time_count, ".dat");
		 fp_rseloss = fopen(filename2_rs, "a");
		 if(!fp_rseloss) {

		   printf("\nCould not open the output file rseloss. \n");
		   return -1;
		 }
		*/


		// Calculate electron energy loss due to radiative processes.
		// Losses due to Accrestion Disk calculated here for both full & Thomson
		// cross-sections depending on the value of gamma_max
		//
		#pragma omp parallel for
		 for (k = 0; k < EGRID; k++) {

		   if ((Ldisk != 0.0) && ((LDISK_fs != 0) || (LDISK_rs != 0))) {

		     /*
                    ecdthloss[k] = fabs(
                            gammaecdthdot(gamma_e[k], Gamma_sh, z_c, thetar,
                                    raddist));
		     */

                    ecdloss[k] = fabs(
                            gammaecdiskdot(gamma_e[k], Gamma_sh, z_c, thetar,
                                    raddist));

		   } else {

		     //ecdthloss[k] = 0.0;
		     ecdloss[k] = 0.0;

		   }

		 }
		 // END PARALLEL FOR SECTION

		// Losses due to BLR photons calculated here for both full & Thomson cross-section.
		//
		#pragma omp parallel for
		 for (k = 0; k < EGRID; k++) {

		   if ((L_BLR != 0.0) && ((LBLR_fs != 0) || (LBLR_rs != 0))) {

		     /*
		     ecblrthloss[k] = fabs(
					   gamdotblrth(
						       gamma_e[k], -1.0, 0.0, Gamma_sh,
						       eplabinteg_main, lineinten_linmain,
						       continten_linmain, etaph_lin,
						       lineinten_gauspmain, lineinten_gausmmain,
						       weights-1, etaph_gausp, etaph_gausm,
						       continten_gauspmain, continten_gausmmain, GAUSASZ
						       )
					   );
		     */
		     
		     ecblrloss[k] = fabs(
					 gamdotblr(
						   gamma_e[k], -1.0, 0.0, lineepp, lineepm,
						   contepp, contepm, lineeplin, conteplin,
						   nph_linelin, nph_contlin, etaph_lin,
						   nph_linegausp, nph_linegausm, weights-1,
						   etaph_gausp, etaph_gausm, nph_contgausp,
						   nph_contgausm, GAUSASZ
						   )
					 );
		     
		   } else {

		     //ecblrthloss[k] = 0.0;
		     ecblrloss[k] = 0.0;
		     
		   }
		 }
		 // END PARALLEL FOR SECTION

		 // Losses due to DT photons calculated here for both full & Thomson
		 // cross-sections.
		 //
	         #pragma omp parallel for
		 for (k = 0; k < EGRID; k++) {

		   if ((L_DT != 0.0) && ((LDT_fs != 0) || (LDT_rs != 0))) {

		     //ec_dtthloss[k] = fabs(gamdotdtth(gamma_e[k], Gamma_sh, etadtarr, FCOVDT));

		     ec_dtloss[k] = fabs(gamdotdt(gamma_e[k], etadtarr, dtep, nph_dt));

		   } else {

		     //ec_dtthloss[k] = 0.0;
		     ec_dtloss[k] = 0.0;
		     
		   }
		   
		 }
		 // END PARALLEL FOR SECTION for Losses due to DT photons


	// Total losses calculated here.
	//
	//#pragma omp parallel for private(k)
		 for (zn = Nrs; zn < f; zn++) {
		   
#pragma omp parallel for
		   for (k = 0; k < EGRID; k++) {
		     
		     if ((zn == (f - 1)) && (timecount_fs == 0)) {
		       
		       eloss[zn][k] = esynloss_fs[k] + ecdloss[k]
			 + ecblrloss[k] + ec_dtloss[k];
		       
		     } else {

		       esscloss[zn][k] = fabs(gammasscdot(epsil, n_phgg[zn], gamma_e[k], FREQ_GRID));			 
		       
		       /*
			 printf("\n In FS, time_count = %d, Bflag = %e, esscloss[%d][%d] = %e\n", 
			 time_count, Bflag, zn, k, esscloss[zn][k]);
		       */
		       
		       eloss[zn][k] = esynloss_fs[k] + esscloss[zn][k]
			 + ecdloss[k] + ecblrloss[k] + ec_dtloss[k];
		       
		     }
		     
		   }
		   // END PARALLEL FOR SECTION
		 }
		 // END PARALLEL FOR SECTION on zn
		 
		 /*
		   for (zn = Nrs; zn < f; zn++) {
		   
		   for (k = 0; k < EGRID; k++) {
		   
		   fprintf(fp_fseloss, "%e %e %e %e %e %e %e %e %e %e %d\n", gamma_e[k], eloss[zn][k], esynloss_fs[k],
		   max(1.0e-40, esscloss[zn][k]), max(1.0e-40, ecdloss[k]), max(1.0e-40, ecdthloss[k]),
		   max(1.0e-40, ecblrloss[k]), max(1.0e-40, ecblrthloss[k]),
		   max(1.0e-40, ec_dtloss[k]), max(1.0e-40, ec_dtthloss[k]), zn);
		   
		   }
		   
		   fprintf(fp_fseloss, "\n");
		   }
		   
		   fclose(fp_fseloss);
		 */
		 
		 //#pragma omp parallel for private(k)
		 for (zn = (Nrs - 1); zn >= r; zn--) {
		   
#pragma omp parallel for
		   for (k = 0; k < FREQ_GRID; k++) {
		     
		     if ((zn == r) && (timecount_rs == 0)) {
		       
		       eloss[zn][k] = esynloss_rs[k] + ecdloss[k]
			 + ecblrloss[k] + ec_dtloss[k];
		       
		     } else {

		       esscloss[zn][k] = fabs(gammasscdot(epsil, n_phgg[zn], gamma_e[k], FREQ_GRID));

		       /*
			 printf("\nIn RS, time_count = %d, Bflag = %e, esscloss[%d][%d] = %e\n", 
			 time_count, Bflag, zn, k, esscloss[zn][k]);
		       */
		       
		       eloss[zn][k] = esynloss_rs[k] + esscloss[zn][k]
			 + ecdloss[k] + ecblrloss[k] + ec_dtloss[k];
		       
		     }
		     
		   }
		   // END PARALLEL FOR SECTION
		 }
		 // END PARALLEL FOR SECTION on zn
		 
		 /*
		   for (zn = (Nrs - 1); zn >= r; zn--) {
		   
		   for (k = 0; k < EGRID; k++) {
		   
		   fprintf(fp_rseloss, "%e %e %e %e %e %e %e %e %e %e %d\n", gamma_e[k], eloss[zn][k], esynloss_rs[k],
		   max(1.0e-40, esscloss[zn][k]), max(1.0e-40, ecdloss[k]), max(1.0e-40, ecdthloss[k]),
		   max(1.0e-40, ecblrloss[k]), max(1.0e-40, ecblrthloss[k]),
		   max(1.0e-40, ec_dtloss[k]), max(1.0e-40, ec_dtthloss[k]), zn);
		   
		   }
		   
		   fprintf(fp_rseloss, "\n");
		   }
		   
		   fclose(fp_rseloss);
		 */
		 
		 if ((files)
		     && ((f != (Nrs + 1))
			 || ((f == (Nrs + 1)) && (timecount_fs != 0)))) {
		   
		   sprintf(filename4_fs, "%s%d%s", fstauggfile, time_count, ".dat");
		   fp_fstgg = fopen(filename4_fs, "w");
		   if (!fp_fstgg) {
		     
		     printf("\nCould not open the output file fstaugg. \n");
		     return -1;
		   }
		   
		   sprintf(filename5_fs, "%s%d%s", fsscnfnfile, time_count, ".dat");
		   fp_fsscnfn = fopen(filename5_fs, "a");
		   if (!fp_fsscnfn) {
		     
		     printf("\nCould not open the output file fsscnfn. \n");
		     return -1;
		   }
		   
		 }
		 
		 if ((files)
		     && ((r != (Nrs - 1))
			 || ((r == (Nrs - 1)) && (timecount_rs != 0)))) {
		   
		   sprintf(filename5_rs, "%s%d%s", rsscnfnfile, time_count, ".dat");
		   fp_rsscnfn = fopen(filename5_rs, "a");
		   if (!fp_rsscnfn) {
		     
		     printf("\nCould not open the output file rsscnfn. \n");
		     return -1;
		   }
		   
		 }
		 

		 // Calculate and update the SSC spectrum at every time step
		 //
		 //#pragma omp parallel for private(k)
		 for (zn = Nrs; zn < f; zn++) {
		   
		   if ((zn == (f - 1)) && (timecount_fs == 0))
		     continue;
		   
#pragma omp parallel for
		   for (k = 0; k < FREQ_GRID; k++) {
		     
		     ndot_ssc[zn][k] = ndotssc(epsilscatter[k], epsil, gamma_e, n_e_p[zn], 
					       n_phgg[zn], FREQ_GRID);
		       
		     //printf("\nIn FS, time_count = %d, Bflag = %e, ndot_ssc[%d][%d] = %e\n", time_count, Bflag, zn, k, ndot_ssc[zn][k]);
			 
		     power_ssc[zn][k] = 5.355265e-47 * syn_nu[k] * ndot_ssc[zn][k];
		     
		     if (epsilscatter[k] > 1.0e-9) {
		       // epsilscatter[k] greater enough that it gives meaningful value of tau_gg[zn][k].
		       
		       if ((accf == TRUE) && (zn == (f - 1))) {			

			 tau_gg[zn][k] = lmeanshfs_param * taugg(epsilscatter[k], epsil, n_phgg[zn], FREQ_GRID);			

		       } else {			

			 tau_gg[zn][k] = lmeanfs_param * taugg(epsilscatter[k], epsil, n_phgg[zn], FREQ_GRID);	

		       }
		       
		     } else {
		       tau_gg[zn][k] = 0.0;
		     }
		     
		     if (tau_gg[zn][k] > TAULIMIT) {
		       ndot_ssc[zn][k] *= ((1.0 - exp(-tau_gg[zn][k])) / tau_gg[zn][k]);
		       
		     }
		   }
		   // END PARALLEL FOR SECTION
		 }
		 // END PARALLEL FOR SECTION on zn
		 
		 
		 if (files) {
		   for (zn = Nrs; zn < f; zn++) {
		     
		     if ((zn == (f - 1)) && (timecount_fs == 0)) {
		       
		       continue;
		     }
		     
		     for (k = 0; k < FREQ_GRID; k++) {
		       
		       if (zn == (f - 1)) {
			 fprintf(fp_fstgg, "%e %e %e %e\n", epsilscatter[k],
				 ndot_ssc[zn][k], power_ssc[zn][k], tau_gg[zn][k]);
		       }
		       
		       fprintf(fp_fsscnfn, "%e %e %d\n", syn_nu[k], ndot_ssc[zn][k], zn);
		     }
		   }
		 }
		 

		 //#pragma omp parallel for private(k)
		 for (zn = (Nrs - 1); zn >= r; zn--) {
		   
		   if ((zn == r) && (timecount_rs == 0)) {
		     continue;
		   }
		   
#pragma omp parallel for
		   for (k = 0; k < FREQ_GRID; k++) {
		     
		     ndot_ssc[zn][k] = ndotssc(epsilscatter[k], epsil, gamma_e, n_e_p[zn],
					       n_phgg[zn], FREQ_GRID);
		     
		     //printf("\nIn RS, time_count = %d, Bflag = %e, ndot_ssc[%d][%d] = %e\n", time_count, Bflag, zn, k, ndot_ssc[zn][k]);
		     
		     power_ssc[zn][k] = 5.355265e-47 * syn_nu[k] * ndot_ssc[zn][k];
		     
		     if (epsilscatter[k] > 1.0e-9) {
		       // epsilscatter[k] greater enough that it gives meaningful value of tau_gg[zn][k].
		       
		       if ((accr == TRUE) && (zn == r)) {

			 tau_gg[zn][k] = lmeanshrs_param * taugg(epsilscatter[k], epsil, n_phgg[zn], FREQ_GRID);

		       } else {

			 tau_gg[zn][k] = lmeanrs_param * taugg(epsilscatter[k], epsil, n_phgg[zn], FREQ_GRID);

		       }
		       
		     } else {
		       tau_gg[zn][k] = 0.0;
		     }
		     
		     if (tau_gg[zn][k] > TAULIMIT) {
		       ndot_ssc[zn][k] *= ((1.0 - exp(-tau_gg[zn][k])) / tau_gg[zn][k]);
		       
		     }
		   }
		   // END PARALLEL FOR SECTION
		 }
		 // END PARALLEL FOR SECTION ON ZN
		 
		 if (files) {
		   for (zn = (Nrs - 1); zn >= r; zn--) {
		     
		     if ((zn == r) && (timecount_rs == 0)) {
		       continue;
		     }
		     
		     for (k = 0; k < FREQ_GRID; k++) {
		       
		       fprintf(fp_rsscnfn, "%e %e %d\n", syn_nu[k], ndot_ssc[zn][k], zn);
		     }
		   }
		   
		   if (((f != (Nrs + 1)) || ((f == (Nrs + 1)) && (timecount_fs != 0)))) {
		     fclose(fp_fsscnfn);
		     fclose(fp_fstgg);
		   }
		   
		   if (((r != (Nrs - 1)) || ((r == (Nrs - 1)) && (timecount_rs != 0)))) {
		     fclose(fp_rsscnfn);
		   }
		 }
		 
		 
		// Calculate and update the ECD spectrum at every time step if Ldisk, MBH,
		// and acc. efficiency are not zero.
		//
		if (Ldisk != 0.0) {

			if (files) {
			  sprintf(filename9_fs, "%s%d%s", fsecdnfnfile, time_count, ".dat");
				fp_fsecdnfn = fopen(filename9_fs, "a");
				if (!fp_fsecdnfn) {

					printf("\nCould not open the output file fsecd. \n");
					return -1;
				}

				sprintf(filename9_rs, "%s%d%s", rsecdnfnfile, time_count, ".dat");
				fp_rsecdnfn = fopen(filename9_rs, "a");
				if (!fp_rsecdnfn) {

					printf("\nCould not open the output file rsecd. \n");
					return -1;
				}

			}


	  if (LDISK_fs != 0) {
	    
			//#pragma omp parallel for private(k)
			for (zn = Nrs; zn < f; zn++) {

			  #pragma omp parallel for
                    for (k = 0; k < FREQ_GRID; k++) {

		if (((accf == FALSE) && (gam_max[zn] < GAMMAXLIM)) || ((accf == TRUE) && (gamma_max_fs < GAMMAXLIM)) || ((accf == TRUE) && (zn != (f - 1)) && (gam_max[zn] < GAMMAXLIM))) {

		  if ((timecount_fs == 0) && (zn == (f - 1))) {

		    ndot_ecd[zn][k] = ndotth(epsilscatter[k], Gamma_sh, z_c, mu_prime, thetar, raddist, gamma_e, n_e_p0_fs);

		  } else {

		    ndot_ecd[zn][k] = ndotth(epsilscatter[k], Gamma_sh, z_c, mu_prime, thetar, raddist, gamma_e, n_e_p[zn]);

		  }

                        } else {

                            if ((timecount_fs == 0) && (zn == (f - 1))) {

		    ndot_ecd[zn][k] = ndotecd(epsilscatter[k], Gamma_sh, z_c, mu_prime, thetar, raddist, gamma_e, n_e_p0_fs, hypt, etaphd, p_engy, phi_main, phistep);

                            } else {

		    ndot_ecd[zn][k] = ndotecd(epsilscatter[k], Gamma_sh, z_c, mu_prime, thetar, raddist, gamma_e, n_e_p[zn], hypt, etaphd, p_engy, phi_main, phistep);

                            }

			}

			if (ndot_ecd[zn][k] < 0.0) {

			  printf("In FS, ndotecd should not be negative for ndot_ecd[%d][%d] = %e, time_count = %d\n",
				 zn, k, ndot_ecd[zn][k], time_count);
			  exitSimulation (start, logfile, -1, sOutputDirectory, nfnfile);
			}


		power_ecd[zn][k] = 5.355265e-47 * syn_nu[k] * ndot_ecd[zn][k];

                        if (!((timecount_fs == 0) && (zn == (f - 1)))) {
			  if (tau_gg[zn][k] > TAULIMIT) {
		    ndot_ecd[zn][k] *= ((1.0 - exp(-tau_gg[zn][k])) / tau_gg[zn][k]);
			  }
                        }
                    }
                // END PARALLEL FOR SECTION
			}
		// END PARALLEL FOR SECTION ON ZN

	  } else {

	    for (zn = Nrs; zn < f; zn++) {
	      for (k = 0; k < FREQ_GRID; k++) {
		ndot_ecd[zn][k] = 0.0;
	      }
	    }

	  }

	  if (files) {

	    for (zn = Nrs; zn < f; zn++) {
	      for (k = 0; k < FREQ_GRID; k++) {
		
		fprintf(fp_fsecdnfn, "%e %e %e %d\n", syn_nu[k],
			ndot_ecd[zn][k], tau_gg[zn][k], zn);
		}
	      }
	    
	    fclose(fp_fsecdnfn);
	  }
	  

	  if (LDISK_rs != 0) {
	    
			//#pragma omp parallel for private(k)
			for (zn = (Nrs - 1); zn >= r; zn--) {

			  #pragma omp parallel for
	      for (k = 0; k < FREQ_GRID; k++) {

		if (((accr == FALSE) && (gam_max[zn] < GAMMAXLIM)) || ((accr == TRUE) && (gamma_max_rs < GAMMAXLIM)) || ((accr == TRUE) && (zn != r) && (gam_max[zn] < GAMMAXLIM))) {

		  if ((timecount_rs == 0) && (zn == r)) {

		    ndot_ecd[zn][k] = ndotth(epsilscatter[k], Gamma_sh, z_c, mu_prime, thetar, raddist, gamma_e, n_e_p0_rs);

		  } else {

		    ndot_ecd[zn][k] = ndotth(epsilscatter[k], Gamma_sh, z_c, mu_prime, thetar, raddist, gamma_e, n_e_p[zn]);
			      
		  }

		} else {
		  
		  if ((timecount_rs == 0) && (zn == r)) {
		    
		    ndot_ecd[zn][k] = ndotecd(epsilscatter[k], Gamma_sh, z_c, mu_prime, thetar, raddist, gamma_e, n_e_p0_rs, hypt, etaphd, p_engy, phi_main, phistep);
		    
		  } else {
		    
		    ndot_ecd[zn][k] = ndotecd(epsilscatter[k], Gamma_sh, z_c, mu_prime, thetar, raddist, gamma_e, n_e_p[zn], hypt, etaphd, p_engy, phi_main, phistep);
		    
		  }
		  
		}
		
		if (ndot_ecd[zn][k] < 0.0) {
		  
		  printf("In RS, ndotecd should not be negative for ndot_ecd[%d][%d] = %e, time_count = %d\n",
			 zn, k, ndot_ecd[zn][k], time_count);
			      exitSimulation (start, logfile, -1, sOutputDirectory, nfnfile);
		  
		}
		
		power_ecd[zn][k] = 5.355265e-47 * syn_nu[k] * ndot_ecd[zn][k];
			    
		if (!((timecount_rs == 0) && (zn == r))) {
		  if (tau_gg[zn][k] > TAULIMIT) {
		    ndot_ecd[zn][k] *= ((1.0 - exp(-tau_gg[zn][k])) / tau_gg[zn][k]);
		  }
		}
	      }
	      // END PARALLEL FOR SECTION
			}
	    // END PARALLEL FOR SECTION ON ZN
			
	  } else {

            for(zn = (Nrs - 1); zn >= r; zn--) {
	      for (k =0; k < FREQ_GRID; k++) {
		ndot_ecd[zn][k] = 0.0;
	      }
            }

          }


	  if (files) {
	    for (zn = (Nrs - 1); zn >= r; zn--) {
	      for (k = 0; k < FREQ_GRID; k++) {
		
		fprintf(fp_rsecdnfn, "%e %e %e %d\n", syn_nu[k],
			ndot_ecd[zn][k], tau_gg[zn][k], zn);
		  
		}
		
	      }
	    
	    fclose(fp_rsecdnfn);
	  }
	  
		} // END IF (Ldisk != 0.0)


		// Calculate and update the ECBLR spectrum at every time step if L_BLR
		// is not zero.
		//
		if (L_BLR != 0.0) {

			if (files) {
			  sprintf(filename12_fs, "%s%d%s", fsblrnfnfile, time_count, ".dat");
				fp_fsblrnfn = fopen(filename12_fs, "a");
				if (!fp_fsblrnfn) {

					printf("\nCould not open the output file fsblr. \n");
					return -1;
				}

				sprintf(filename12_rs, "%s%d%s", rsblrnfnfile, time_count, ".dat");
				fp_rsblrnfn = fopen(filename12_rs, "a");
				if (!fp_rsblrnfn) {

					printf("\nCould not open the output file rsblr. \n");
					return -1;
				}

			}


			// ECBLR FOR THE FORWARD SHOCK
			//
	  if (LBLR_fs != 0) {

			//#pragma omp parallel for private(k)
			for (zn = Nrs; zn < f; zn++) {

			  #pragma omp parallel for
                    for (k = 0; k < FREQ_GRID; k++) {

		if (((accf == FALSE) && (gam_max[zn] < GAMMAXLIM)) || ((accf == TRUE) && (gamma_max_fs < GAMMAXLIM)) || ((accf == TRUE) && (zn != (f - 1)) && (gam_max[zn] < GAMMAXLIM))) {

		  if ((timecount_fs == 0) && (zn == (f - 1))) {

		    ndot_blr[zn][k] = ndotblrth(epsilscatter[k], Gamma_sh,
						eplabinteg_main, -1.0, 0.0,
						weights-1, etaph_gausp, etaph_gausm,
						etaph_lin, lineNV, gamma_e, n_e_p0_fs,
						lineepp, lineepm, lineeplin,
						continten_gauspmain, continten_gausmmain,
						continten_linmain, lineinten_gauspmain,
						lineinten_gausmmain, lineinten_linmain,
						GAUSASZ, coskip_main, coskim_main,
						coskilin_main, phistep);

		  } else {

		    ndot_blr[zn][k] = ndotblrth(epsilscatter[k], Gamma_sh,
						eplabinteg_main, -1.0, 0.0,
						weights-1, etaph_gausp, etaph_gausm,
						etaph_lin, lineNV, gamma_e, n_e_p[zn],
						lineepp, lineepm, lineeplin,
						continten_gauspmain, continten_gausmmain,
						continten_linmain, lineinten_gauspmain,
						lineinten_gausmmain, lineinten_linmain,
						GAUSASZ, coskip_main, coskim_main,
						coskilin_main, phistep);
		  }
		  
		} else {

		  if ((timecount_fs == 0) && (zn == (f - 1))) {
		    
		    ndot_blr[zn][k] = ndotblr(epsilscatter[k],
					      mu_prime, gamma_e, n_e_p0_fs, -1.0, 0.0,
					      lineepp, lineepm, contepp, contepm, lineeplin,
					      conteplin, nph_linelin, nph_contlin, etaph_lin,
					      nph_linegausp, nph_linegausm, weights-1,
					      etaph_gausp, etaph_gausm, nph_contgausp,
					      nph_contgausm, GAUSASZ, coskip_main,
					      coskim_main, coskilin_main, phistep);
		    
		  } else {
		    
		    ndot_blr[zn][k] = ndotblr(epsilscatter[k],
					      mu_prime, gamma_e, n_e_p[zn], -1.0, 0.0,
					      lineepp, lineepm, contepp, contepm, lineeplin,
					      conteplin, nph_linelin, nph_contlin, etaph_lin,
					      nph_linegausp, nph_linegausm, weights-1,
					      etaph_gausp, etaph_gausm, nph_contgausp,
					      nph_contgausm, GAUSASZ, coskip_main,
					      coskim_main, coskilin_main, phistep);
		    
		  }
		  
		}
		
		
		if (ndot_blr[zn][k] < 0.0) {
		  
		  printf("In FS, ndotblr should not be negative for ndot_blr[%d][%d] = %e, time_count = %d\n",
			 zn, k, ndot_blr[zn][k], time_count);
		  exitSimulation (start, logfile, -1, sOutputDirectory, nfnfile);
		  
		}

		power_blr[zn][k] = 5.355265e-47 * syn_nu[k] * ndot_blr[zn][k];

                        if (!((timecount_fs == 0) && (zn == (f - 1)))) {

                            if (tau_gg[zn][k] > TAULIMIT) {
			      ndot_blr[zn][k] *= ((1.0 - exp(-tau_gg[zn][k])) / tau_gg[zn][k]);
                            }
                        }
                    }
                // END PARALLEL FOR SECTION
			}
		// END PARALLEL FOR SECTION ON ZN

	  } else {

            for(zn = Nrs; zn < f; zn++) {
	      for (k = 0; k < FREQ_GRID; k++) {
		ndot_blr[zn][k] = 0.0;
	      }
            }

          }

	  
            if (files) {
                for (zn = Nrs; zn < f; zn++) {
                    for (k = 0; k < FREQ_GRID; k++) {

		fprintf(fp_fsblrnfn, "%e %e %e %d\n", syn_nu[k],
			ndot_blr[zn][k], tau_gg[zn][k], zn);
		      }
                    }

		fclose(fp_fsblrnfn);
            }


			// ECBLR FOR REVERSE SHOCK
			//
	  if (LBLR_rs != 0) {
	    
	    //#pragma omp parallel for private(k)
			for (zn = (Nrs - 1); zn >= r; zn--) {

			#pragma omp parallel for
                    for (k = 0; k < FREQ_GRID; k++) {

		if (((accr == FALSE) && (gam_max[zn] < GAMMAXLIM)) || ((accr == TRUE) && (gamma_max_rs < GAMMAXLIM)) || ((accr == TRUE) && (zn != r) && (gam_max[zn] < GAMMAXLIM))) {

		  if ((timecount_rs == 0) && (zn == r)) {

		    ndot_blr[zn][k] = ndotblrth(epsilscatter[k], Gamma_sh, eplabinteg_main,
						-1.0, 0.0, weights-1, etaph_gausp, etaph_gausm,
						etaph_lin, lineNV, gamma_e, n_e_p0_rs, lineepp,
						lineepm, lineeplin, continten_gauspmain,
						continten_gausmmain, continten_linmain,
						lineinten_gauspmain, lineinten_gausmmain,
						lineinten_linmain, GAUSASZ, coskip_main,
						coskim_main, coskilin_main, phistep);

		  } else {

		    ndot_blr[zn][k] = ndotblrth(epsilscatter[k], Gamma_sh, eplabinteg_main,
						-1.0, 0.0, weights-1, etaph_gausp, etaph_gausm,
						etaph_lin, lineNV, gamma_e, n_e_p[zn], lineepp,
						lineepm, lineeplin, continten_gauspmain,
						continten_gausmmain, continten_linmain,
						lineinten_gauspmain, lineinten_gausmmain,
						lineinten_linmain, GAUSASZ, coskip_main,
						coskim_main, coskilin_main, phistep);
		    
		  }
		  
		} else {
		  
		  if ((timecount_rs == 0) && (zn == r)) {
		    
		    ndot_blr[zn][k] = ndotblr(epsilscatter[k],
					      mu_prime, gamma_e, n_e_p0_rs, -1.0, 0.0,
					      lineepp, lineepm, contepp, contepm, lineeplin,
					      conteplin, nph_linelin, nph_contlin, etaph_lin,
					      nph_linegausp, nph_linegausm, weights-1,
					      etaph_gausp, etaph_gausm, nph_contgausp,
					      nph_contgausm, GAUSASZ, coskip_main,
					      coskim_main, coskilin_main, phistep);
		    
		  } else {
		    
		    ndot_blr[zn][k] = ndotblr(epsilscatter[k],
					      mu_prime, gamma_e, n_e_p[zn], -1.0, 0.0,
					      lineepp, lineepm, contepp, contepm, lineeplin,
					      conteplin, nph_linelin, nph_contlin, etaph_lin,
					      nph_linegausp, nph_linegausm, weights-1,
					      etaph_gausp, etaph_gausm, nph_contgausp,
					      nph_contgausm, GAUSASZ, coskip_main,
					      coskim_main, coskilin_main, phistep);
		    
		  }
		  
		}
		
		if (ndot_blr[zn][k] < 0.0) {
		  
		  printf("In RS, ndotblr should not be negative for ndot_blr[%d][%d] = %e, time_count = %d\n",
			 zn, k, ndot_blr[zn][k], time_count);
		  exitSimulation (start, logfile, -1, sOutputDirectory, nfnfile);
		}

		power_blr[zn][k] = 5.355265e-47 * syn_nu[k] * ndot_blr[zn][k];

		if (!((timecount_rs == 0) && (zn == r))) {
		  
		  if (tau_gg[zn][k] > TAULIMIT) {
		    ndot_blr[zn][k] *= ((1.0 - exp(-tau_gg[zn][k])) / tau_gg[zn][k]);

		  }
		}
                    }
                // END PARALLEL FOR SECTION
			}
		// END PARALLEL FOR SECTION ON ZN

	  } else {

            for(zn = (Nrs - 1); zn >= r; zn--) {
	      for (k = 0; k < FREQ_GRID; k++) {
		ndot_blr[zn][k] = 0.0;
	      }
            }

          }


			if (files) {
			    for (zn = (Nrs - 1); zn >= r; zn--) {
			        for (k = 0; k < FREQ_GRID; k++) {

		fprintf(fp_rsblrnfn, "%e %e %e %d\n", syn_nu[k],
			ndot_blr[zn][k], tau_gg[zn][k], zn);
				  }
			        }
	    
			    fclose(fp_rsblrnfn);
			}

		} // END IF (L_BLR != 0.0)


		// Calculate and update the ECDT spectrum at every time step if L_DT
		// is not zero.
		//
		if (L_DT != 0.0) {

			if (files) {

				sprintf(filename12_fs, "%s%d%s", fsdtnfnfile, time_count, ".dat");
				fp_fsdtnfn = fopen(filename12_fs, "a");

				if(!fp_fsdtnfn) {

					printf("\nCould not open the output file fsdt. \n");
					return -1;
				}

				sprintf(filename12_rs, "%s%d%s", rsdtnfnfile, time_count, ".dat");
				fp_rsdtnfn = fopen(filename12_rs, "a");

				if(!fp_rsdtnfn) {

					printf("\nCould not open the output file rsdt. \n");
					return -1;
				}
			}


	  if (LDT_fs != 0) {
	    
			//#pragma omp parallel for private(k)
			for (zn = Nrs; zn < f; zn++) {

                          #pragma omp parallel for
			  for (k = 0; k < FREQ_GRID; k++) {
		
		if (((accf == FALSE) && (gam_max[zn] < GAMMAXDTLIM)) || ((accf == TRUE) && (gamma_max_fs < GAMMAXDTLIM)) || ((accf == TRUE) && (zn != (f - 1)) && (gam_max[zn] < GAMMAXDTLIM))) {
		  
		  if ((timecount_fs == 0) && (zn == (f - 1))) {

		    ndot_dt[zn][k] = ndotdtth(epsilscatter[k], Gamma_sh, mu_prime, phi_main, etadtarr, gamma_e, n_e_p0_fs, phistep, FCOVDT);

		  } else {

		    ndot_dt[zn][k] = ndotdtth(epsilscatter[k], Gamma_sh, mu_prime, phi_main, etadtarr, gamma_e, n_e_p[zn], phistep, FCOVDT);
		  
		  }

		} else {
		  
		  if ((timecount_fs == 0) && (zn == (f - 1))) {
		    
		    ndot_dt[zn][k] = ndotdt(epsilscatter[k], mu_prime, gamma_e, n_e_p0_fs, etadtarr, phi_main, dtep, nph_dt, phistep);
		    
		  } else {
		    
		    ndot_dt[zn][k] = ndotdt(epsilscatter[k], mu_prime, gamma_e, n_e_p[zn], etadtarr, phi_main, dtep, nph_dt, phistep);
		    
		  }
		  
		}
		
                if (ndot_dt[zn][k] < 0.0) {

		  printf("In FS, ndotdt should not be negative for ndot_dt[%d][%d] = %e, time_count = %d\n",
                         zn, k, ndot_dt[zn][k], time_count);
                  exitSimulation (start, logfile, -1, sOutputDirectory, nfnfile);
                }
		
		power_dt[zn][k] = 5.355265e-47 * syn_nu[k] * ndot_dt[zn][k];
		
		if (!((timecount_fs == 0) && (zn == (f - 1)))) {
		  
		  if (tau_gg[zn][k] > TAULIMIT) {
		    
		    ndot_dt[zn][k] *= ((1.0 - exp(-tau_gg[zn][k])) / tau_gg[zn][k]);
		  }
			    }
			  }
			  // END PARALLEL FOR SECTION
			} // END PARALLEL FOR SECTION FOR ZN

	  } else {

            for(zn = Nrs; zn < f; zn++) {
	      for (k = 0; k < FREQ_GRID; k++) {
		ndot_dt[zn][k] = 0.0;
	      }
            }

          }

			
	  if (files) {
	    for (zn = Nrs; zn < f; zn++) {
	      for (k = 0; k < FREQ_GRID; k++) {
		
		fprintf(fp_fsdtnfn, "%e %e %e %d\n", syn_nu[k], ndot_dt[zn][k], tau_gg[zn][k], zn);

	      }
	    }
	    
	    fclose(fp_fsdtnfn);
	  }
	  

	  if (LDT_rs != 0) {

			//#pragma omp parallel for private(k)
			for (zn = (Nrs - 1); zn >= r; zn--) {

                          #pragma omp parallel for
			  for (k = 0; k < FREQ_GRID; k++) {

		if (((accr == FALSE) && (gam_max[zn] < GAMMAXDTLIM)) || ((accr == TRUE) && (gamma_max_rs < GAMMAXDTLIM)) || ((accr == TRUE) && (zn != r) && (gam_max[zn] < GAMMAXDTLIM))) {
		  
		  if ((timecount_rs == 0) && (zn == r)) {

		    ndot_dt[zn][k] = ndotdtth(epsilscatter[k], Gamma_sh, mu_prime, phi_main, etadtarr, gamma_e, n_e_p0_rs, phistep, FCOVDT);

		  } else {

		  ndot_dt[zn][k] = ndotdtth(epsilscatter[k], Gamma_sh, mu_prime, phi_main, etadtarr, gamma_e, n_e_p[zn], phistep, FCOVDT);
		  
		  }

		} else {
		  
		  if ((timecount_rs == 0) && (zn == r)) {
		    
		    ndot_dt[zn][k] = ndotdt(epsilscatter[k], mu_prime, gamma_e, n_e_p0_rs, etadtarr, phi_main, dtep, nph_dt, phistep);
		    
		  } else {
		    
		    ndot_dt[zn][k] = ndotdt(epsilscatter[k], mu_prime, gamma_e, n_e_p[zn], etadtarr, phi_main, dtep, nph_dt, phistep);
		    
		  }
		  
		}
		
		if (ndot_dt[zn][k] < 0.0) {

                  printf("In RS, ndotdt should not be negative for ndot_dt[%d][%d] = %e, time_count = %d\n",
                         zn, k, ndot_dt[zn][k], time_count);
                  exitSimulation (start, logfile, -1, sOutputDirectory, nfnfile);
		}
		
		power_dt[zn][k] = 5.355265e-47 * syn_nu[k] * ndot_dt[zn][k];
		
		if (!((timecount_rs == 0) && (zn == r))) {
		  
		  if (tau_gg[zn][k] > TAULIMIT) {
		    
		    ndot_dt[zn][k] *= ((1.0 - exp(-tau_gg[zn][k])) / tau_gg[zn][k]);
		  }
		}
			  } // END PARALLEL FOR SECTION
			} // END PARALLEL FOR SECTION FOR ZN

	  } else {

            for(zn = (Nrs - 1); zn >= r; zn--) {
	      for (k = 0; k < FREQ_GRID; k++) {
		ndot_dt[zn][k] = 0.0;
	      }
            }

          }


	  if (files) {
	    for (zn = (Nrs - 1); zn >= r; zn--) {
	      for (k = 0; k < FREQ_GRID; k++) {
		
		fprintf(fp_rsdtnfn, "%e %e %e %d\n", syn_nu[k], ndot_dt[zn][k], tau_gg[zn][k], zn);
		
	      }
	    }
	    
	    fclose(fp_rsdtnfn);
	  }
	  
		} // END IF (L_DT != 0.0)


		// Calculate electron population due to pair production
		//
		//#pragma omp parallel for private(j)
		for (zn = Nrs; zn < f; zn++) {
		  
		  if ((zn == (f - 1)) && (timecount_fs == 0))
		    continue;
		  
#pragma omp parallel for
		  for (j = 0; j < EGRID; j++) {
		    
		    ndot_gg[zn][j] = nedotgg(epsil, gamma_e[j], n_phgg[zn], FREQ_GRID);

		  } //End of parallel section
		} //End of parallel section on zn
		
		//#pragma omp parallel for private(j)
		for (zn = (Nrs - 1); zn >= r; zn--) {
		  
		  if ((zn == r) && (timecount_rs == 0))
		    continue;
		  
#pragma omp parallel for
		  for (j = 0; j < EGRID; j++) {

		    ndot_gg[zn][j] = nedotgg(epsil, gamma_e[j], n_phgg[zn], FREQ_GRID);
		    
		  } //End of parallel section
		} //End of parallel section on zn


		/*
		 sprintf(filename3_fs, "%s%d%s", fsnewedenfile, time_count, ".dat");
		 fp_fsnew = fopen(filename3_fs, "a");
		 if(!fp_fsnew) {

		 printf("\nCould not open the output file fsneweden. \n");
		 return -1;
		 }

		 sprintf(filename3_rs, "%s%d%s", rsnewedenfile, time_count, ".dat");
		 fp_rsnew = fopen(filename3_rs, "a");
		 if(!fp_rsnew) {

		 printf("\nCould not open the output file rsneweden. \n");
		 return -1;
		 }
		*/


		// Calculate electron spectrum using tridiagonal system of equations
		// method
		//
		dt_param = 1.0 + (dt / te_esc);
		steps_param = dt / (steps - 1.0);

		for (zn = Nrs; zn < f; zn++) {

			for (j = 0; j < EGRID; j++) {

				if (j == (EGRID - 1)) {

					b_i[zn][j] = dt_param
							+ (steps_param * eloss[zn][j] / gamma_e[j]);
					c_i[zn][j] = -steps_param * eloss[zn][j] / gamma_e[j];

				} else {

					b_i[zn][j] = dt_param
							+ (dt * eloss[zn][j]
							                  / (gamma_e[j + 1] - gamma_e[j]));
					c_i[zn][j] = -dt * eloss[zn][j + 1]
					                             / (gamma_e[j + 1] - gamma_e[j]);

				}

				if ((accf == TRUE) && (zn == (f - 1))) {

				  if (timecount_fs == 0) {
				    r_i[zn][j] = Q_inj_fs[j] * dt;

				  } else {
				    r_i[zn][j] = n_e[zn][j] + ((Q_inj_fs[j] + ndot_gg[zn][j]) * dt);

				  }

				} else {
				  r_i[zn][j] = n_e[zn][j] + (ndot_gg[zn][j] * dt);

				}

				u_i[zn][j] = 0.0;

			}

		}

		for (zn = (Nrs - 1); zn >= r; zn--) {

			for (j = 0; j < EGRID; j++) {

				if (j == (EGRID - 1)) {

					b_i[zn][j] = dt_param
							+ (steps_param * eloss[zn][j] / gamma_e[j]);
					c_i[zn][j] = -steps_param * eloss[zn][j] / gamma_e[j];

				} else {

					b_i[zn][j] = dt_param
							+ (dt * eloss[zn][j]
							                  / (gamma_e[j + 1] - gamma_e[j]));
					c_i[zn][j] = -dt * eloss[zn][j + 1]
					                             / (gamma_e[j + 1] - gamma_e[j]);

				}

				if ((accr == TRUE) && (zn == r)) {

				  if (timecount_rs == 0) {
				    r_i[zn][j] = Q_inj_rs[j] * dt;

				  } else {
				    r_i[zn][j] = n_e[zn][j] + ((Q_inj_rs[j] + ndot_gg[zn][j]) * dt);

				  }

				} else {
				  r_i[zn][j] = n_e[zn][j] + (ndot_gg[zn][j] * dt);

				}

				u_i[zn][j] = 0.0;

			}

		}

		for (zn = Nrs; zn < f; zn++) {
			tridag(a_i, b_i[zn], c_i[zn], r_i[zn], u_i[zn], EGRID);
		}

		for (zn = (Nrs - 1); zn >= r; zn--) {
			tridag(a_i, b_i[zn], c_i[zn], r_i[zn], u_i[zn], EGRID);
		}

		// Result of the tridag routine: new electron densities
		//
		for (zn = Nrs; zn < f; zn++) {

		  for (j = 0; j < EGRID; j++) {

		    if (u_i[zn][j] < 1.0e-40)
		      u_i[zn][j] = 0.0;
		    /* Some small enough limit that below that it is as good as zero. */

		    n_e[zn][j] = u_i[zn][j];
		    n_e_p[zn][j] = 2.0 * n_e[zn][j];

		    if ((accf == TRUE) && (zn == (f - 1)))
		      ne_curnt_fs[j] = n_e[zn][j];

		    /*
		      fprintf(fp_fsnew, "%e %e %e %e %e %e %d\n", r_i[zn][j], a_i[j],
			      b_i[zn][j], c_i[zn][j], u_i[zn][j], ndot_gg[zn][j], zn);
		    */

		  }

		  //fprintf(fp_fsnew, "\n");

		  e_numden[zn] = ne_numtot(gamma_e, n_e[zn], EGRID);

		  /*
		    if (fmod(write_count, 10) == 0) {
		    fp_log = fopen(logfile, "a");
		    if(!fp_log) {
		    printf("\nCould not open the output file log. \n");
		    return -1;
		    }

		    fprintf(fp_log, "time count = %d, e_numden[%d] = %e\n",
		    time_count, zn, e_numden[zn]);

		    fclose(fp_log);

		    }
		  */

		}

		for (zn = (Nrs - 1); zn >= r; zn--) {

		  for (j = 0; j < EGRID; j++) {

		    if (u_i[zn][j] < 1.0e-40)
		      u_i[zn][j] = 0.0;

		    n_e[zn][j] = u_i[zn][j];
		    n_e_p[zn][j] = 2.0 * n_e[zn][j];

		    if ((accr == TRUE) && (zn == r))
		      ne_curnt_rs[j] = n_e[zn][j];

		    /*
		      fprintf(fp_rsnew, "%e %e %e %e %e %e %d\n", r_i[zn][j], a_i[j],
			      b_i[zn][j], c_i[zn][j], u_i[zn][j], ndot_gg[zn][j], zn);
		    */

		  }

		  //fprintf(fp_rsnew, "\n");

		  e_numden[zn] = ne_numtot(gamma_e, n_e[zn], EGRID);

		  /*
		    if (fmod(write_count, 10) == 0) {
		    fp_log = fopen(logfile, "a");
		    if(!fp_log) {
		    printf("\nCould not open the output file log. \n");
		    return -1;
		    }

		    fprintf(fp_log, "time count = %d, e_numden[%d] = %e\n",
		    time_count, zn, e_numden[zn]);

		    fclose(fp_log);

		    }
		  */

		}

		//fclose(fp_fsnew);
		//fclose(fp_rsnew);


		// Vary gamma_min and gamma_max for both FS and RS regions, with time,
		// once the shocks have exited the zone or the system.
		//
		totgammadotmin_fs = 1.0e-10;

		for (zn = Nrs; zn < f; zn++) {

		  if ((accf == FALSE) || ((accf == TRUE) && (zn != (f - 1)))) {

		    gammadot_synmin[zn] = fabs(gammasyndot(gam_min[zn], B_fs));

		    gammadot_sscmin[zn] = fabs(gammasscdot(epsil, n_phgg[zn], gam_min[zn], FREQ_GRID));

		    //printf("In FS, Bflag = %e, gammadot_sscmin[%d] = %e\n", Bflag, zn, gammadot_sscmin[zn]);

		    if ((Ldisk != 0.0) && (LDISK_fs != 0)) {
		      
		      gammadot_ecdmin[zn] = fabs(gammaecdiskdot(gam_min[zn], Gamma_sh, z_c, thetar, raddist));

		    } else {

		      gammadot_ecdmin[zn] = 0.0;

		    }
		    
		    if ((L_BLR != 0.0) && (LBLR_fs != 0)) {
		      
		      gammadot_blrmin[zn] = fabs(gamdotblr(gam_min[zn], -1.0, 0.0, lineepp,
							   lineepm, contepp, contepm, lineeplin,
							   conteplin, nph_linelin, nph_contlin,
							   etaph_lin, nph_linegausp, nph_linegausm,
							   weights-1, etaph_gausp, etaph_gausm,
							   nph_contgausp, nph_contgausm, GAUSASZ));
		      
		    } else {
		      
		      gammadot_blrmin[zn] = 0.0;
		      
		    }
		    
		    if ((L_DT != 0.0) && (LDT_fs != 0)) {
		      
		      gammadot_dtmin[zn] = fabs(gamdotdt(gam_min[zn], etadtarr, dtep, nph_dt));
		      
		    } else {
		      
		      gammadot_dtmin[zn] = 0.0;
		      
		    }
		    
		    totgammadotmin_zn_fs = gammadot_synmin[zn] + gammadot_sscmin[zn] +
		      gammadot_ecdmin[zn] + gammadot_blrmin[zn] + gammadot_dtmin[zn];
		    
		    gam_min[zn] = gam_min[zn] - (dt * totgammadotmin_zn_fs);
		    gam_max[zn] = gam_max[zn] - (dt * totgammadotmax[zn]);

		    if (fmod(write_count, 10) == 0) {
		      fp_log = fopen(logfile, "a");
		      if(!fp_log) {

			printf("\nCould not open the output file log. \n");
			return -1;

		      }

		      fprintf(fp_log, "totgammadotmin_zn_fs = %e, gam_min[%d] = %e, syndotmin = %e, sscdotmin = %e, ecddotmin = %e, blrdotmin = %e, dtdotmin = %e\n",
			      totgammadotmin_zn_fs, zn, gam_min[zn], gammadot_synmin[zn], gammadot_sscmin[zn], gammadot_ecdmin[zn], gammadot_blrmin[zn], gammadot_dtmin[zn]);

	              fprintf(fp_log, "time counter = %d, gam_min[%d] = %e, gam_max[%d] = %e\n",
			      time_count, zn, gam_min[zn], zn, gam_max[zn]);

		      fclose(fp_log);

		    }

		    if (gam_min[zn] < 1.01) gam_min[zn] = 1.01;
		    if (gam_max[zn] < 1.01) gam_max[zn] = 1.01;
		    if (gam_max[zn] < gam_min[zn]) gam_max[zn] = gam_min[zn];

		  }

		  if ((accf == FALSE) && (accr == FALSE)) {

		    if ((gam_min[zn] <= 1.01) || (gam_max[zn] <= 1.01) || (gam_max[zn] <= gam_min[zn])) {

		      printf("\nFS: gam_min[%d] = %e, gam_max[%d] = %e, time_count = %d\n",
			     zn, gam_min[zn], zn, gam_max[zn], time_count);
		      printf("The minimum and/or maximum gamma values are out of range\n");
		      exitSimulation (start, logfile, -1, sOutputDirectory, nfnfile);

		    }

		  }

		}

		for (zn = (Nrs - 1); zn >= r; zn--) {

		  if ((accr == FALSE) || ((accr == TRUE) && (zn != r))) {

		    gammadot_synmin[zn] = fabs(gammasyndot(gam_min[zn], B_rs));

		    gammadot_sscmin[zn] = fabs(gammasscdot(epsil, n_phgg[zn], gam_min[zn], FREQ_GRID));

		    //printf("In RS, Bflag = %e, gammadot_sscmin[%d] = %e\n", Bflag, zn, gammadot_sscmin[zn]);
		    
		    if ((Ldisk != 0.0) && (LDISK_rs != 0)) {
		      
		      gammadot_ecdmin[zn] = fabs(gammaecdiskdot(gam_min[zn], Gamma_sh, z_c, thetar, raddist));
		      
		    } else {
		      
		      gammadot_ecdmin[zn] = 0.0;
		      
		    }
		    
		    if ((L_BLR != 0.0) && (LBLR_rs != 0)) {
		      
		      gammadot_blrmin[zn] = fabs(gamdotblr(gam_min[zn], -1.0, 0.0, lineepp,
							   lineepm, contepp, contepm, lineeplin,
							   conteplin, nph_linelin, nph_contlin,
							   etaph_lin, nph_linegausp, nph_linegausm,
							   weights-1, etaph_gausp, etaph_gausm,
							   nph_contgausp, nph_contgausm, GAUSASZ));
		      
		    } else {
		      
		      gammadot_blrmin[zn] = 0.0;
		      
		    }
		    
		    if ((L_DT != 0.0) && (LDT_rs != 0)) {
		      
		      gammadot_dtmin[zn] = fabs(gamdotdt(gam_min[zn], etadtarr, dtep, nph_dt));
		      
		    } else {
		      
		      gammadot_dtmin[zn] = 0.0;
		      
		    }
		    
		    totgammadotmin_zn_rs = gammadot_synmin[zn] + gammadot_sscmin[zn] +
		      gammadot_ecdmin[zn] + gammadot_blrmin[zn] + gammadot_dtmin[zn];

		    gam_min[zn] = gam_min[zn] - (dt * totgammadotmin_zn_rs);
		    gam_max[zn] = gam_max[zn] - (dt * totgammadotmax[zn]);

		    if (fmod(write_count, 10) == 0) {
		      fp_log = fopen(logfile, "a");
		      if(!fp_log) {

			printf("\nCould not open the output file log. \n");
			return -1;

		      }

		      fprintf(fp_log, "totgammadotmin_zn_rs = %e, gam_min[%d] = %e, syndotmin = %e, sscdotmin = %e, ecddotmin = %e, blrdotmin = %e, dtdotmin = %e\n",
			      totgammadotmin_zn_rs, zn, gam_min[zn], gammadot_synmin[zn], gammadot_sscmin[zn], gammadot_ecdmin[zn], gammadot_blrmin[zn], gammadot_dtmin[zn]);

		      fprintf(fp_log, "time_count = %d, gam_min[%d] = %e, gam_max[%d] = %e\n",
			      time_count, zn, gam_min[zn], zn, gam_max[zn]);

		      fclose(fp_log);

		    }

		    if (gam_min[zn] < 1.01) gam_min[zn] = 1.01;
		    if (gam_max[zn] < 1.01) gam_max[zn] = 1.01;
		    if (gam_max[zn] < gam_min[zn]) gam_max[zn] = gam_min[zn];

		  }

		  if ((accf == FALSE) && (accr == FALSE)) {

		    if ((gam_min[zn] <= 1.01) || (gam_max[zn] <= 1.01) || (gam_max[zn] <= gam_min[zn])) {

		      printf("\nRS: gam_min[%d] = %e, gam_max[%d] = %e, time_count = %d\n",
			     zn, gam_min[zn], zn, gam_max[zn], time_count);
		      printf("The minimum and maximum gamma values are out of range\n");
		      exitSimulation (start, logfile, -1, sOutputDirectory, nfnfile);
		    }

		  }

		}


		/*
		 sprintf(filename6_fs, "%s%d%s", fsphotfile, time_count, ".dat");
		 fp_fsphot = fopen(filename6_fs, "a");
		 if(!fp_fsphot) {

		 printf("\nCould not open the output file fsphot. \n");
		 return -1;
		 }

		 sprintf(filename7_fs, "%s%d%s", fscomphotfile, time_count, ".dat");
		 fp_fscomphot = fopen(filename7_fs, "a");
		 if(!fp_fscomphot) {

		 printf("\nCould not open the output file fscomphot. \n");
		 return -1;
		 }

		 sprintf(filename6_rs, "%s%d%s", rsphotfile, time_count, ".dat");
		 fp_rsphot = fopen(filename6_rs, "a");
		 if(!fp_rsphot) {

		 printf("\nCould not open the output file rsphot. \n");
		 return -1;
		 }

		 sprintf(filename7_rs, "%s%d%s", rscomphotfile, time_count, ".dat");
		 fp_rscomphot = fopen(filename7_rs, "a");
		 if(!fp_rscomphot) {

		 printf("\nCould not open the output file rscomphot. \n");
		 return -1;
		 }
		*/


		// Take feedback from within the forward region and from the reverse region.
		//
		for (zn = Nrs; zn < f; zn++) {

		  for (j = 0; j < FREQ_GRID; j++) {

                    if ((zn == (Ntot - 1)) || ((Ntot == 2) && (zn == 1))) {
		      
		      // To be used for completely randomly oriented field and for SSC and 
		      // pair production calculations as they do not have dependence on the 
		      // orientation of the field:
                      ndot_phgg[zn][j] = ndot_syn_ep[zn][j] + ndot_ssc[zn][j]
                        + ndot_ecd[zn][j] + ndot_blr[zn][j] + ndot_dt[zn][j]
                        + ndot_phfeed_up[zn - 1][j];

		      // To be used for what the observer would see as the final synchrotron
		      // output is dependent on ordered & disordered field components:
		      ndot_ph[zn][j] = ndot_syn_tot[zn][j] + ndot_ssc[zn][j] 
			+ ndot_ecd[zn][j] + ndot_blr[zn][j] + ndot_dt[zn][j] 
			+ ndot_ph_up[zn - 1][j];
		      
		      // To be used for feedback to adjacent zones that includes the anisotropic 
		      // synchrotron radiation going into the adjacent zone and the 
		      // synchrotron radiation due to the disordered component: 
		      if (Bflag != 0.0) {

			ndot_phfeed[zn][j] = ndot_syn_unifeedtot[zn][j] + ndot_ssc[zn][j]
			  + ndot_ecd[zn][j] + ndot_blr[zn][j] + ndot_dt[zn][j] 
			  + ndot_phfeed_up[zn - 1][j];

		      } else ndot_phfeed[zn][j] = ndot_phgg[zn][j];
		      
                    } else {

                      ndot_phgg[zn][j] = ndot_syn_ep[zn][j] + ndot_ssc[zn][j]
                        + ndot_ecd[zn][j] + ndot_blr[zn][j] + ndot_dt[zn][j]
                        + ndot_phfeed_up[zn - 1][j] + ndot_phfeed_do[zn + 1][j];

		      ndot_ph[zn][j] = ndot_syn_tot[zn][j] + ndot_ssc[zn][j]
			+ ndot_ecd[zn][j] + ndot_blr[zn][j] + ndot_dt[zn][j] 
			+ ndot_ph_up[zn - 1][j] + ndot_ph_do[zn + 1][j];
			
		      if (Bflag != 0.0) {
			
			ndot_phfeed[zn][j] = ndot_syn_unifeedtot[zn][j] + ndot_ssc[zn][j] 
			  + ndot_ecd[zn][j] + ndot_blr[zn][j] + ndot_dt[zn][j]
			  + ndot_phfeed_up[zn - 1][j] + ndot_phfeed_do[zn + 1][j];

		      } else ndot_phfeed[zn][j] = ndot_phgg[zn][j];

                    }

		  }

		}


		// Take feedback from within the reverse region and from the forward region.
		//
		for (zn = (Nrs - 1); zn >= r; zn--) {

		  for (j = 0; j < FREQ_GRID; j++) {

                    if ((zn == 0) || ((Ntot == 2) && (zn == 0))) {

                      // To be used for completely randomly oriented field and for SSC and 
		      // pair production calculations as they do not
                      // have dependence on the orientation of the field:
                      ndot_phgg[zn][j] = ndot_syn_ep[zn][j] + ndot_ssc[zn][j]
                        + ndot_ecd[zn][j] + ndot_blr[zn][j] + ndot_dt[zn][j]
                        + ndot_phfeed_do[zn + 1][j];

		      // To be used for what the observer would see as the final synchrotron
		      // output is dependent on ordered & disordered field components:
		      ndot_ph[zn][j] = ndot_syn_tot[zn][j] + ndot_ssc[zn][j] 
			+ ndot_ecd[zn][j] + ndot_blr[zn][j] + ndot_dt[zn][j] 
			+ ndot_ph_do[zn + 1][j];
		      
		      // To be used for feedback to adjacent zones that includes the anisotropic
		      // synchrotron radiation going into the adjacent zone and the
		      // synchrotron radiation due to the disordered component:
		      if (Bflag != 0.0) {
			
			ndot_phfeed[zn][j] = ndot_syn_unifeedtot[zn][j] + ndot_ssc[zn][j] 
			  + ndot_ecd[zn][j] + ndot_blr[zn][j] + ndot_dt[zn][j] 
			  + ndot_phfeed_do[zn + 1][j];
			
		      } else ndot_phfeed[zn][j] = ndot_phgg[zn][j];
		      
                    } else {

                      ndot_phgg[zn][j] = ndot_syn_ep[zn][j] + ndot_ssc[zn][j]
                        + ndot_ecd[zn][j] + ndot_blr[zn][j] + ndot_dt[zn][j]
                        + ndot_phfeed_up[zn - 1][j] + ndot_phfeed_do[zn + 1][j];

		      ndot_ph[zn][j] = ndot_syn_tot[zn][j] + ndot_ssc[zn][j] 
			+ ndot_ecd[zn][j] + ndot_blr[zn][j] + ndot_dt[zn][j] 
			+ ndot_ph_up[zn - 1][j] + ndot_ph_do[zn + 1][j];
		      
		      if (Bflag != 0.0) {
			
			ndot_phfeed[zn][j] = ndot_syn_unifeedtot[zn][j] + ndot_ssc[zn][j] 
			  + ndot_ecd[zn][j] + ndot_blr[zn][j] + ndot_dt[zn][j]
			  + ndot_phfeed_up[zn - 1][j] + ndot_phfeed_do[zn + 1][j];

		      } else ndot_phfeed[zn][j] = ndot_phgg[zn][j];

                    }

		  }

		}


		// Solve the photon continuity equation for the forward region after taking
		// into account the complete feedback.
		//
		for (zn = Nrs; zn < f; zn++) {

		  for (j = 0; j < FREQ_GRID; j++) {
		    
		    n_ph[zn][j] += (ndot_ph[zn][j] * dt);
		    n_phgg[zn][j] += (ndot_phgg[zn][j] * dt);
		    n_phfeed[zn][j] += (ndot_phfeed[zn][j] * dt);

		    n_phsyn[zn][j] += (ndot_syn_tot[zn][j] * dt);
		    n_phsynep[zn][j] += (ndot_syn_ep[zn][j] * dt);

		    if (Bflag != 0.0) {
		      
		      n_phsynfeed[zn][j] += (ndot_syn_unifeedtot[zn][j] * dt);

		    }

		    if ((Bflag == 4.0) || (Bflag == 5.0)) {

		      for (p = 0; p < QUANTIZE; p++) {

			n_phsynord[zn][j][p] += (ndot_syn_ord[zn][j][p] * dt);

		      }

		    } 

		    n_phssc[zn][j] += (ndot_ssc[zn][j] * dt);
		    n_phecd[zn][j] += (ndot_ecd[zn][j] * dt);
		    n_phblr[zn][j] += (ndot_blr[zn][j] * dt);
		    n_phdt[zn][j] += (ndot_dt[zn][j] * dt);

		    if (time_count != 0) {
		      
		      n_phup[zn - 1][j] += (ndot_ph_up[zn - 1][j] * dt);
		      n_phggup[zn - 1][j] += (ndot_phgg_up[zn - 1][j] * dt);
		      n_phfeedup[zn - 1][j] += (ndot_phfeed_up[zn - 1][j] * dt);
		      
		      if (zn != (Ntot - 1)) {
			
			n_phdo[zn + 1][j] += (ndot_ph_do[zn + 1][j] * dt);
			n_phggdo[zn + 1][j] += (ndot_phgg_do[zn + 1][j] * dt);
			n_phfeeddo[zn + 1][j] += (ndot_phfeed_do[zn + 1][j] * dt);
			
		      }

		    }

		    if ((accf == TRUE) && (zn == (f - 1))) {

		      ndot_ph_up[zn][j] = n_ph[zn][j] * pupshfs_param;
		      ndot_ph_do[zn][j] = ndot_ph_up[zn][j];
		      ndot_ph_si[zn][j] = n_ph[zn][j] * psideshfs_param;

		      ndot_phgg_up[zn][j] = n_phgg[zn][j] * pupshfs_param;
		      ndot_phgg_do[zn][j] = ndot_phgg_up[zn][j];
		      ndot_phgg_si[zn][j] = n_phgg[zn][j] * psideshfs_param;

		      ndot_phfeed_up[zn][j] = n_phfeed[zn][j] * pupshfs_param;
		      ndot_phfeed_do[zn][j] = ndot_phfeed_up[zn][j];
		      ndot_phfeed_si[zn][j] = n_phfeed[zn][j] * psideshfs_param;

                      ndot_phsyn_up[zn][j] = n_phsyn[zn][j] * pupshfs_param;
		      ndot_phsyn_do[zn][j] = ndot_phsyn_up[zn][j];
                      ndot_phsyn_si[zn][j] = n_phsyn[zn][j] * psideshfs_param;

		      ndot_phsynep_up[zn][j] = n_phsynep[zn][j] * pupshfs_param;
		      ndot_phsynep_do[zn][j] = ndot_phsynep_up[zn][j];
                      ndot_phsynep_si[zn][j] = n_phsynep[zn][j] * psideshfs_param;

		      if (Bflag != 0.0) {

			ndot_phsynfeed_up[zn][j] = n_phsynfeed[zn][j] * pupshfs_param;
			ndot_phsynfeed_do[zn][j] = ndot_phsynfeed_up[zn][j];
			ndot_phsynfeed_si[zn][j] = n_phsynfeed[zn][j] * psideshfs_param;

		      }

		      if ((Bflag == 4.0) || (Bflag == 5.0)) {

			for (p = 0; p < QUANTIZE; p++) {

			  ndot_synord_up[zn][j][p] = n_phsynord[zn][j][p] * pupshfs_param;
			  ndot_synord_do[zn][j][p] = ndot_synord_up[zn][j][p];
			  ndot_synord_si[zn][j][p] = n_phsynord[zn][j][p] * psideshfs_param;

			}

		      } 

                      ndot_phssc_up[zn][j] = n_phssc[zn][j] * pupshfs_param;
		      ndot_phssc_do[zn][j] = ndot_phssc_up[zn][j];
                      ndot_phssc_si[zn][j] = n_phssc[zn][j] * psideshfs_param;

                      ndot_phecd_up[zn][j] = n_phecd[zn][j] * pupshfs_param;
		      ndot_phecd_do[zn][j] = ndot_phecd_up[zn][j];
                      ndot_phecd_si[zn][j] = n_phecd[zn][j] * psideshfs_param;

                      ndot_phblr_up[zn][j] = n_phblr[zn][j] * pupshfs_param;
		      ndot_phblr_do[zn][j] = ndot_phblr_up[zn][j];
                      ndot_phblr_si[zn][j] = n_phblr[zn][j] * psideshfs_param;

                      ndot_phdt_up[zn][j] = n_phdt[zn][j] * pupshfs_param;
		      ndot_phdt_do[zn][j] = ndot_phdt_up[zn][j];
                      ndot_phdt_si[zn][j] = n_phdt[zn][j] * psideshfs_param;
		      
		      if (time_count != 0) {

			ndot_phup_up[zn - 1][j] = n_phup[zn - 1][j] * pupshfs_param;
			ndot_phup_do[zn - 1][j] = ndot_phup_up[zn - 1][j];
			ndot_phup_si[zn - 1][j] = n_phup[zn - 1][j] * psideshfs_param;

			ndot_phggup_up[zn - 1][j] = n_phggup[zn - 1][j] * pupshfs_param;
			ndot_phggup_do[zn - 1][j] = ndot_phggup_up[zn - 1][j];
			ndot_phggup_si[zn - 1][j] = n_phggup[zn - 1][j] * psideshfs_param;

			ndot_phfeedup_up[zn - 1][j] = n_phfeedup[zn - 1][j] * pupshfs_param;
                        ndot_phfeedup_do[zn - 1][j] = ndot_phfeedup_up[zn - 1][j];
                        ndot_phfeedup_si[zn - 1][j] = n_phfeedup[zn - 1][j] * psideshfs_param;
			
			if (zn != (Ntot - 1)) {
			  
			  ndot_phdo_up[zn + 1][j] = n_phdo[zn + 1][j] * pupshfs_param;
			  ndot_phdo_do[zn + 1][j] = ndot_phdo_up[zn + 1][j];
			  ndot_phdo_si[zn + 1][j] = n_phdo[zn + 1][j] * psideshfs_param;

			  ndot_phggdo_up[zn + 1][j] = n_phggdo[zn + 1][j] * pupshfs_param;
			  ndot_phggdo_do[zn + 1][j] = ndot_phggdo_up[zn + 1][j];
			  ndot_phggdo_si[zn + 1][j] = n_phggdo[zn + 1][j] * psideshfs_param;

			  ndot_phfeeddo_up[zn + 1][j] = n_phfeeddo[zn + 1][j] * pupshfs_param;
                          ndot_phfeeddo_do[zn + 1][j] = ndot_phfeeddo_up[zn + 1][j];
                          ndot_phfeeddo_si[zn + 1][j] = n_phfeeddo[zn + 1][j] * psideshfs_param;
			  
			}

		      }
		      
		    } else {

		      ndot_ph_up[zn][j] = n_ph[zn][j] * pupfs_param;
		      ndot_ph_do[zn][j] = ndot_ph_up[zn][j];
		      ndot_ph_si[zn][j] = n_ph[zn][j] * psidefs_param;

		      ndot_phgg_up[zn][j] = n_phgg[zn][j] * pupfs_param;
		      ndot_phgg_do[zn][j] = ndot_phgg_up[zn][j];
		      ndot_phgg_si[zn][j] = n_phgg[zn][j] * psidefs_param;

		      ndot_phfeed_up[zn][j] = n_phfeed[zn][j] * pupfs_param;
		      ndot_phfeed_do[zn][j] = ndot_phfeed_up[zn][j];
		      ndot_phfeed_si[zn][j] = n_phfeed[zn][j] * psidefs_param;

                      ndot_phsyn_up[zn][j] = n_phsyn[zn][j] * pupfs_param;
		      ndot_phsyn_do[zn][j] = ndot_phsyn_up[zn][j];
                      ndot_phsyn_si[zn][j] = n_phsyn[zn][j] * psidefs_param;

		      ndot_phsynep_up[zn][j] = n_phsynep[zn][j] * pupfs_param;
		      ndot_phsynep_do[zn][j] = ndot_phsynep_up[zn][j];
                      ndot_phsynep_si[zn][j] = n_phsynep[zn][j] * psidefs_param;

		      if (Bflag != 0.0) {
			
			ndot_phsynfeed_up[zn][j] = n_phsynfeed[zn][j] * pupfs_param;
			ndot_phsynfeed_do[zn][j] = ndot_phsynfeed_up[zn][j];
			ndot_phsynfeed_si[zn][j] = n_phsynfeed[zn][j] * psidefs_param; 

		      }

                      if ((Bflag == 4.0) || (Bflag == 5.0)) {

			for (p = 0; p < QUANTIZE; p++) {

                          ndot_synord_up[zn][j][p] = n_phsynord[zn][j][p] * pupfs_param;
                          ndot_synord_do[zn][j][p] = ndot_synord_up[zn][j][p];
                          ndot_synord_si[zn][j][p] = n_phsynord[zn][j][p] * psidefs_param;

			}

                      } 

                      ndot_phssc_up[zn][j] = n_phssc[zn][j] * pupfs_param;
		      ndot_phssc_do[zn][j] = ndot_phssc_up[zn][j];
                      ndot_phssc_si[zn][j] = n_phssc[zn][j] * psidefs_param;

                      ndot_phecd_up[zn][j] = n_phecd[zn][j] * pupfs_param;
		      ndot_phecd_do[zn][j] = ndot_phecd_up[zn][j];
                      ndot_phecd_si[zn][j] = n_phecd[zn][j] * psidefs_param;

                      ndot_phblr_up[zn][j] = n_phblr[zn][j] * pupfs_param;
		      ndot_phblr_do[zn][j] = ndot_phblr_up[zn][j];
                      ndot_phblr_si[zn][j] = n_phblr[zn][j] * psidefs_param;

                      ndot_phdt_up[zn][j] = n_phdt[zn][j] * pupfs_param;
		      ndot_phdt_do[zn][j] = ndot_phdt_up[zn][j];
                      ndot_phdt_si[zn][j] = n_phdt[zn][j] * psidefs_param;

		      if (time_count != 0) {

			ndot_phup_up[zn - 1][j] = n_phup[zn - 1][j] * pupfs_param;
			ndot_phup_do[zn - 1][j] = ndot_phup_up[zn - 1][j];
			ndot_phup_si[zn - 1][j] = n_phup[zn - 1][j] * psidefs_param;

			ndot_phggup_up[zn - 1][j] = n_phggup[zn - 1][j] * pupfs_param;
			ndot_phggup_do[zn - 1][j] = ndot_phggup_up[zn - 1][j];
			ndot_phggup_si[zn - 1][j] = n_phggup[zn - 1][j] * psidefs_param;

			ndot_phfeedup_up[zn - 1][j] = n_phfeedup[zn - 1][j] * pupfs_param;
                        ndot_phfeedup_do[zn - 1][j] = ndot_phfeedup_up[zn - 1][j];
                        ndot_phfeedup_si[zn - 1][j] = n_phfeedup[zn - 1][j] * psidefs_param;
			
			if (zn != (Ntot - 1)) {

			  ndot_phdo_up[zn + 1][j] = n_phdo[zn + 1][j] * pupfs_param;
			  ndot_phdo_do[zn + 1][j] = ndot_phdo_up[zn + 1][j];
			  ndot_phdo_si[zn + 1][j] = n_phdo[zn + 1][j] * psidefs_param;

			  ndot_phggdo_up[zn + 1][j] = n_phggdo[zn + 1][j] * pupfs_param;
			  ndot_phggdo_do[zn + 1][j] = ndot_phggdo_up[zn + 1][j];
			  ndot_phggdo_si[zn + 1][j] = n_phggdo[zn + 1][j] * psidefs_param;

			  ndot_phfeeddo_up[zn + 1][j] = n_phfeeddo[zn + 1][j] * pupfs_param;
                          ndot_phfeeddo_do[zn + 1][j] = ndot_phfeeddo_up[zn + 1][j];
                          ndot_phfeeddo_si[zn + 1][j] = n_phfeeddo[zn + 1][j] * psidefs_param;
			  
			}

		      }

		    }

		    n_ph[zn][j] = n_ph[zn][j] - (((2.0 * ndot_ph_up[zn][j]) + ndot_ph_si[zn][j]) * dt);		    

		    n_phgg[zn][j] = n_phgg[zn][j] - (((2.0 * ndot_phgg_up[zn][j]) + ndot_phgg_si[zn][j]) * dt);

		    n_phfeed[zn][j] = n_phfeed[zn][j] - (((2.0 * ndot_phfeed_up[zn][j]) + ndot_phfeed_si[zn][j]) * dt);

		    n_phsyn[zn][j] = n_phsyn[zn][j] - (((2.0 * ndot_phsyn_up[zn][j]) + ndot_phsyn_si[zn][j]) * dt);
		    
		    n_phsynep[zn][j] = n_phsynep[zn][j] - (((2.0 * ndot_phsynep_up[zn][j]) + ndot_phsynep_si[zn][j]) * dt);
		    
		    if (Bflag != 0.0) {

		      n_phsynfeed[zn][j] = n_phsynfeed[zn][j] - (((2.0 * ndot_phsynfeed_up[zn][j]) + ndot_phsynfeed_si[zn][j]) * dt);

		    }

		    if ((Bflag == 4.0) || (Bflag == 5.0)) {

		      for (p = 0; p < QUANTIZE; p++) {

			n_phsynord[zn][j][p] = n_phsynord[zn][j][p] - (((2.0 * ndot_synord_up[zn][j][p]) + ndot_synord_si[zn][j][p]) * dt);
			//printf("\n In FS, time_count = %d, n_phsynord[%d][%d][%d] = %e\n", time_count, zn, j, p, n_phsynord[zn][j][p]);

		      }

		    } 

		    n_phssc[zn][j] = n_phssc[zn][j] - (((2.0 * ndot_phssc_up[zn][j]) + ndot_phssc_si[zn][j]) * dt);
		    n_phecd[zn][j] = n_phecd[zn][j] - (((2.0 * ndot_phecd_up[zn][j]) + ndot_phecd_si[zn][j]) * dt);
		    n_phblr[zn][j] = n_phblr[zn][j] - (((2.0 * ndot_phblr_up[zn][j]) + ndot_phblr_si[zn][j]) * dt);
		    n_phdt[zn][j] = n_phdt[zn][j] - (((2.0 * ndot_phdt_up[zn][j]) + ndot_phdt_si[zn][j]) * dt);

		    if (time_count != 0) {

		      n_phup[zn - 1][j] = n_phup[zn - 1][j] - (((2.0 * ndot_phup_up[zn - 1][j]) + ndot_phup_si[zn - 1][j]) * dt);
		      n_phggup[zn - 1][j] = n_phggup[zn - 1][j] - (((2.0 * ndot_phggup_up[zn - 1][j]) + ndot_phggup_si[zn - 1][j]) * dt);
		      n_phfeedup[zn - 1][j] = n_phfeedup[zn - 1][j] - (((2.0 * ndot_phfeedup_up[zn - 1][j]) + ndot_phfeedup_si[zn - 1][j]) * dt);
		      
		      if (zn != (Ntot - 1)) {
			n_phdo[zn + 1][j] = n_phdo[zn + 1][j] - (((2.0 * ndot_phdo_up[zn + 1][j]) + ndot_phdo_si[zn + 1][j]) * dt);
			n_phggdo[zn + 1][j] = n_phggdo[zn + 1][j] - (((2.0 * ndot_phggdo_up[zn + 1][j]) + ndot_phggdo_si[zn + 1][j]) * dt);
			n_phfeeddo[zn + 1][j] = n_phfeeddo[zn + 1][j] - (((2.0 * ndot_phfeeddo_up[zn + 1][j]) + ndot_phfeeddo_si[zn + 1][j]) * dt);
		      }

		    }

		    /*
		    fprintf(fp_fsphot, "%e %e %e %e %e %e %d\n", epsil[j], n_ph[zn][j], n_phgg[zn][j], ndot_ph[zn][j],
			    ndot_ph_up[zn][j], ndot_ph_si[zn][j], zn);

		      nuLnu_si[zn][j] = volfs_param * SQR(syn_nu[j]) * ndot_ph_si[zn][j];
		      nuLnu_up[j] = volfs_param * SQR(syn_nu[j]) * ndot_ph_up[zn][j];
		      nuLnu[zn][j] = nuLnu_si[zn][j] + nuLnu_up[j];

		      fprintf(fp_fscomphot, "%e %e %e %e %d\n", syn_nu[j], nuLnu[zn][j], nuLnu_si[zn][j],
		      nuLnu_up[j], zn);
		    */

		  }

		  //fprintf(fp_fsphot, "\n");

		}

		//fclose(fp_fsphot);
		//fclose(fp_fscomphot);

		// Solve the photon continuity equation for the reverse region after taking
		// into account the complete feedback.
		//
		for (zn = (Nrs - 1); zn >= r; zn--) {

		  for (j = 0; j < FREQ_GRID; j++) {

		    n_ph[zn][j] += (ndot_ph[zn][j] * dt);
		    n_phgg[zn][j] += (ndot_phgg[zn][j] * dt);
		    n_phfeed[zn][j] += (ndot_phfeed[zn][j] * dt);

                    n_phsyn[zn][j] += (ndot_syn_tot[zn][j] * dt);
		    n_phsynep[zn][j] += (ndot_syn_ep[zn][j] * dt);

		    if (Bflag != 0.0) {
		      
		      n_phsynfeed[zn][j] += (ndot_syn_unifeedtot[zn][j] * dt);

		    }

                    if ((Bflag == 4.0) || (Bflag == 5.0)) {

                      for (p = 0; p < QUANTIZE; p++) {

			n_phsynord[zn][j][p] += (ndot_syn_ord[zn][j][p] * dt);

		      }

                    } 

                    n_phssc[zn][j] += (ndot_ssc[zn][j] * dt);
                    n_phecd[zn][j] += (ndot_ecd[zn][j] * dt);
                    n_phblr[zn][j] += (ndot_blr[zn][j] * dt);
                    n_phdt[zn][j] += (ndot_dt[zn][j] * dt);

		    if (time_count != 0) {
		      
		      n_phdo[zn + 1][j] += (ndot_ph_do[zn + 1][j] * dt);
		      n_phggdo[zn + 1][j] += (ndot_phgg_do[zn + 1][j] * dt);
		      n_phfeeddo[zn + 1][j] += (ndot_phfeed_do[zn + 1][j] * dt);
		      
		      if (zn != 0) {
			n_phup[zn - 1][j] += (ndot_ph_up[zn - 1][j] * dt);
			n_phggup[zn - 1][j] += (ndot_phgg_up[zn - 1][j] * dt);
			n_phfeedup[zn - 1][j] += (ndot_phfeed_up[zn - 1][j] * dt);
		      }

		    }

		    if ((accr == TRUE) && (zn == r)) {

		      ndot_ph_up[zn][j] = n_ph[zn][j] * pupshrs_param;
		      ndot_ph_do[zn][j] = ndot_ph_up[zn][j];
		      ndot_ph_si[zn][j] = n_ph[zn][j] * psideshrs_param;

		      ndot_phgg_up[zn][j] = n_phgg[zn][j] * pupshrs_param;
		      ndot_phgg_do[zn][j] = ndot_phgg_up[zn][j];
		      ndot_phgg_si[zn][j] = n_phgg[zn][j] * psideshrs_param;

		      ndot_phfeed_up[zn][j] = n_phfeed[zn][j] * pupshrs_param;
		      ndot_phfeed_do[zn][j] = ndot_phfeed_up[zn][j];
		      ndot_phfeed_si[zn][j] = n_phfeed[zn][j] * psideshrs_param;
		      
                      ndot_phsyn_up[zn][j] = n_phsyn[zn][j] * pupshrs_param;
		      ndot_phsyn_do[zn][j] = ndot_phsyn_up[zn][j];
                      ndot_phsyn_si[zn][j] = n_phsyn[zn][j] * psideshrs_param;

		      ndot_phsynep_up[zn][j] = n_phsynep[zn][j] * pupshrs_param;
		      ndot_phsynep_do[zn][j] = ndot_phsynep_up[zn][j];
                      ndot_phsynep_si[zn][j] = n_phsynep[zn][j] * psideshrs_param;

		      if (Bflag != 0.0) {

			ndot_phsynfeed_up[zn][j] = n_phsynfeed[zn][j] * pupshrs_param;
			ndot_phsynfeed_do[zn][j] = ndot_phsynfeed_up[zn][j];
			ndot_phsynfeed_si[zn][j] = n_phsynfeed[zn][j] * psideshrs_param;

		      }

                      if ((Bflag == 4.0) || (Bflag == 5.0)) {

			for (p = 0; p < QUANTIZE; p++) {

                          ndot_synord_up[zn][j][p] = n_phsynord[zn][j][p] * pupshrs_param;
			  ndot_synord_do[zn][j][p] = ndot_synord_up[zn][j][p];
                          ndot_synord_si[zn][j][p] = n_phsynord[zn][j][p] * psideshrs_param;

                        }

                      } 

                      ndot_phssc_up[zn][j] = n_phssc[zn][j] * pupshrs_param;
		      ndot_phssc_do[zn][j] = ndot_phssc_up[zn][j];
                      ndot_phssc_si[zn][j] = n_phssc[zn][j] * psideshrs_param;

                      ndot_phecd_up[zn][j] = n_phecd[zn][j] * pupshrs_param;
		      ndot_phecd_do[zn][j] = ndot_phecd_up[zn][j];
                      ndot_phecd_si[zn][j] = n_phecd[zn][j] * psideshrs_param;

                      ndot_phblr_up[zn][j] = n_phblr[zn][j] * pupshrs_param;
		      ndot_phblr_do[zn][j] = ndot_phblr_up[zn][j];
                      ndot_phblr_si[zn][j] = n_phblr[zn][j] * psideshrs_param;

                      ndot_phdt_up[zn][j] = n_phdt[zn][j] * pupshrs_param;
		      ndot_phdt_do[zn][j] = ndot_phdt_up[zn][j];
                      ndot_phdt_si[zn][j] = n_phdt[zn][j] * psideshrs_param;

		      if (time_count != 0) {

			ndot_phdo_up[zn + 1][j] = n_phdo[zn + 1][j] * pupshrs_param;
			ndot_phdo_do[zn + 1][j] = ndot_phdo_up[zn + 1][j];
			ndot_phdo_si[zn + 1][j] = n_phdo[zn + 1][j] * psideshrs_param;

			ndot_phggdo_up[zn + 1][j] = n_phggdo[zn + 1][j] * pupshrs_param;
			ndot_phggdo_do[zn + 1][j] = ndot_phggdo_up[zn + 1][j];
			ndot_phggdo_si[zn + 1][j] = n_phggdo[zn + 1][j] * psideshrs_param;

			ndot_phfeeddo_up[zn + 1][j] = n_phfeeddo[zn + 1][j] * pupshrs_param;
                        ndot_phfeeddo_do[zn + 1][j] = ndot_phfeeddo_up[zn + 1][j];
                        ndot_phfeeddo_si[zn + 1][j] = n_phfeeddo[zn + 1][j] * psideshrs_param;
			
			if (zn != 0) {
			  
			  ndot_phup_up[zn - 1][j] = n_phup[zn - 1][j] * pupshrs_param;
			  ndot_phup_do[zn - 1][j] = ndot_phup_up[zn - 1][j];
			  ndot_phup_si[zn - 1][j] = n_phup[zn - 1][j] * psideshrs_param;

			  ndot_phggup_up[zn - 1][j] = n_phggup[zn - 1][j] * pupshrs_param;
			  ndot_phggup_do[zn - 1][j] = ndot_phggup_up[zn - 1][j];
			  ndot_phggup_si[zn - 1][j] = n_phggup[zn - 1][j] * psideshrs_param;

			  ndot_phfeedup_up[zn - 1][j] = n_phfeedup[zn - 1][j] * pupshrs_param;
                          ndot_phfeedup_do[zn - 1][j] = ndot_phfeedup_up[zn - 1][j];
                          ndot_phfeedup_si[zn - 1][j] = n_phfeedup[zn - 1][j] * psideshrs_param;
			  
			}

		      }

		    } else {

		      ndot_ph_up[zn][j] = n_ph[zn][j] * puprs_param;
		      ndot_ph_do[zn][j] = ndot_ph_up[zn][j];
		      ndot_ph_si[zn][j] = n_ph[zn][j] * psiders_param;

		      ndot_phgg_up[zn][j] = n_phgg[zn][j] * puprs_param;
		      ndot_phgg_do[zn][j] = ndot_phgg_up[zn][j];
		      ndot_phgg_si[zn][j] = n_phgg[zn][j] * psiders_param;

		      ndot_phfeed_up[zn][j] = n_phfeed[zn][j] * puprs_param;
		      ndot_phfeed_do[zn][j] = ndot_phfeed_up[zn][j];
		      ndot_phfeed_si[zn][j] = n_phfeed[zn][j] * psiders_param;

                      ndot_phsyn_up[zn][j] = n_phsyn[zn][j] * puprs_param;
		      ndot_phsyn_do[zn][j] = ndot_phsyn_up[zn][j];
                      ndot_phsyn_si[zn][j] = n_phsyn[zn][j] * psiders_param;

		      ndot_phsynep_up[zn][j] = n_phsynep[zn][j] * puprs_param;
		      ndot_phsynep_do[zn][j] = ndot_phsynep_up[zn][j];
                      ndot_phsynep_si[zn][j] = n_phsynep[zn][j] * psiders_param;

		      if (Bflag != 0.0) {

			ndot_phsynfeed_up[zn][j] = n_phsynfeed[zn][j] * puprs_param;
			ndot_phsynfeed_do[zn][j] = ndot_phsynfeed_up[zn][j];
			ndot_phsynfeed_si[zn][j] = n_phsynfeed[zn][j] * psiders_param;

		      }

                      if ((Bflag == 4.0) || (Bflag == 5.0)) {

			for (p = 0; p < QUANTIZE; p++) {

                          ndot_synord_up[zn][j][p] = n_phsynord[zn][j][p] * puprs_param;
                          ndot_synord_do[zn][j][p] = ndot_synord_up[zn][j][p];
                          ndot_synord_si[zn][j][p] = n_phsynord[zn][j][p] * psiders_param;

                        }

                      } 

                      ndot_phssc_up[zn][j] = n_phssc[zn][j] * puprs_param;
		      ndot_phssc_do[zn][j] = ndot_phssc_up[zn][j];
                      ndot_phssc_si[zn][j] = n_phssc[zn][j] * psiders_param;

                      ndot_phecd_up[zn][j] = n_phecd[zn][j] * puprs_param;
		      ndot_phecd_do[zn][j] = ndot_phecd_up[zn][j];
                      ndot_phecd_si[zn][j] = n_phecd[zn][j] * psiders_param;

		      ndot_phblr_up[zn][j] = n_phblr[zn][j] * puprs_param;
		      ndot_phblr_do[zn][j] = ndot_phblr_up[zn][j];
                      ndot_phblr_si[zn][j] = n_phblr[zn][j] * psiders_param;

                      ndot_phdt_up[zn][j] = n_phdt[zn][j] * puprs_param;
		      ndot_phdt_do[zn][j] = ndot_phdt_up[zn][j];
                      ndot_phdt_si[zn][j] = n_phdt[zn][j] * psiders_param;

		      if (time_count != 0) {

			ndot_phdo_up[zn + 1][j] = n_phdo[zn + 1][j] * puprs_param;
			ndot_phdo_do[zn + 1][j] = ndot_phdo_up[zn + 1][j];
			ndot_phdo_si[zn + 1][j] = n_phdo[zn + 1][j] * psiders_param;

			ndot_phggdo_up[zn + 1][j] = n_phggdo[zn + 1][j] * puprs_param;
			ndot_phggdo_do[zn + 1][j] = ndot_phggdo_up[zn + 1][j];
			ndot_phggdo_si[zn + 1][j] = n_phggdo[zn + 1][j] * psiders_param;

			ndot_phfeeddo_up[zn + 1][j] = n_phfeeddo[zn + 1][j] * puprs_param;
                        ndot_phfeeddo_do[zn + 1][j] = ndot_phfeeddo_up[zn + 1][j];
                        ndot_phfeeddo_si[zn + 1][j] = n_phfeeddo[zn + 1][j] * psiders_param;
			
			if (zn != 0) {
			  
			  ndot_phup_up[zn - 1][j] = n_phup[zn - 1][j] * puprs_param;
			  ndot_phup_do[zn - 1][j] = ndot_phup_up[zn - 1][j];
			  ndot_phup_si[zn - 1][j] = n_phup[zn - 1][j] * psiders_param;

			  ndot_phggup_up[zn - 1][j] = n_phggup[zn - 1][j] * puprs_param;
			  ndot_phggup_do[zn - 1][j] = ndot_phggup_up[zn - 1][j];
			  ndot_phggup_si[zn - 1][j] = n_phggup[zn - 1][j] * psiders_param;

			  ndot_phfeedup_up[zn - 1][j] = n_phfeedup[zn - 1][j] * puprs_param;
                          ndot_phfeedup_do[zn - 1][j] = ndot_phfeedup_up[zn - 1][j];
                          ndot_phfeedup_si[zn - 1][j] = n_phfeedup[zn - 1][j] * psiders_param;
			  
			}

		      }

		    }

		    n_ph[zn][j] = n_ph[zn][j] - (((2.0 * ndot_ph_up[zn][j]) + ndot_ph_si[zn][j]) * dt);

		    n_phgg[zn][j] = n_phgg[zn][j] - (((2.0 * ndot_phgg_up[zn][j]) + ndot_phgg_si[zn][j]) * dt); 

		    n_phfeed[zn][j] = n_phfeed[zn][j] - (((2.0 * ndot_phfeed_up[zn][j]) + ndot_phfeed_si[zn][j]) * dt);

		    n_phsyn[zn][j] = n_phsyn[zn][j] - (((2.0 * ndot_phsyn_up[zn][j]) + ndot_phsyn_si[zn][j]) * dt);

		    n_phsynep[zn][j] = n_phsynep[zn][j] - (((2.0 * ndot_phsynep_up[zn][j]) + ndot_phsynep_si[zn][j]) * dt);
		    
		    if (Bflag != 0.0) {

		      n_phsynfeed[zn][j] = n_phsynfeed[zn][j] - (((2.0 * ndot_phsynfeed_up[zn][j]) + ndot_phsynfeed_si[zn][j]) * dt);

		    }

                    if ((Bflag == 4.0) || (Bflag == 5.0)) {

                      for (p = 0; p < QUANTIZE; p++) {

                        n_phsynord[zn][j][p] = n_phsynord[zn][j][p] - (((2.0 * ndot_synord_up[zn][j][p]) + ndot_synord_si[zn][j][p]) * dt);
			//printf("\n In RS, time_count = %d, n_phsynord[%d][%d][%d] = %e\n", time_count, zn, j, p, n_phsynord[zn][j][p]);
			
                      }

                    } 

		    n_phssc[zn][j] = n_phssc[zn][j] - (((2.0 * ndot_phssc_up[zn][j]) + ndot_phssc_si[zn][j]) * dt);
		    n_phecd[zn][j] = n_phecd[zn][j] - (((2.0 * ndot_phecd_up[zn][j]) + ndot_phecd_si[zn][j]) * dt);
		    n_phblr[zn][j] = n_phblr[zn][j] - (((2.0 * ndot_phblr_up[zn][j]) + ndot_phblr_si[zn][j]) * dt);
		    n_phdt[zn][j] = n_phdt[zn][j] - (((2.0 * ndot_phdt_up[zn][j]) + ndot_phdt_si[zn][j]) * dt);

		    if (time_count != 0) {
		      
		      n_phdo[zn + 1][j] = n_phdo[zn + 1][j] - (((2.0 * ndot_phdo_up[zn + 1][j]) + ndot_phdo_si[zn + 1][j]) * dt);
		      n_phggdo[zn + 1][j] = n_phggdo[zn + 1][j] - (((2.0 * ndot_phggdo_up[zn + 1][j]) + ndot_phggdo_si[zn + 1][j]) * dt);
		      n_phfeeddo[zn + 1][j] = n_phfeeddo[zn + 1][j] - (((2.0 * ndot_phfeeddo_up[zn + 1][j]) + ndot_phfeeddo_si[zn + 1][j]) * dt);
		      
		      if (zn != 0) {			
			n_phup[zn - 1][j] = n_phup[zn - 1][j] - (((2.0 * ndot_phup_up[zn - 1][j]) + ndot_phup_si[zn - 1][j]) * dt);
			n_phggup[zn - 1][j] = n_phggup[zn - 1][j] - (((2.0 * ndot_phggup_up[zn - 1][j]) + ndot_phggup_si[zn - 1][j]) * dt);
			n_phfeedup[zn - 1][j] = n_phfeedup[zn - 1][j] - (((2.0 * ndot_phfeedup_up[zn - 1][j]) + ndot_phfeedup_si[zn - 1][j]) * dt);

		      }

		    }

		    /*
		    fprintf(fp_rsphot, "%e %e %e %e %e %e %d\n", epsil[j], n_ph[zn][j], n_phgg[zn][j], ndot_ph[zn][j],
			    ndot_ph_up[zn][j], ndot_ph_si[zn][j], zn);

		      nuLnu_si[zn][j] = volrs_param * SQR(syn_nu[j]) * ndot_ph_si[zn][j];
		      nuLnu_do[j] = volrs_param * SQR(syn_nu[j]) * ndot_ph_do[zn][j];
		      nuLnu[zn][j] = nuLnu_si[zn][j] + nuLnu_do[j];

		      fprintf(fp_rscomphot, "%e %e %e %e %d\n", syn_nu[j], nuLnu[zn][j], nuLnu_si[zn][j],
		      nuLnu_do[j], zn);
		    */

		  }

		  //fprintf(fp_rsphot, "\n");

		}

		//fclose(fp_rsphot);
		//fclose(fp_rscomphot);


		// Post-processing starts here in order to calculate the observed spectrum with the
		// correct time-delay taken into account
		//
		// Calculate the full spectrum in the observer's frame recieved from the forward region
		//
		if (mu_prime > 0.0) {

		  for (i = 0; i < FREQ_GRID; i++) {

		    nuFnu_fs_si[i] = nuFnu_fs_up[i] = nuFnu_rs_si[i] = 0.0;
		    nuFnusyn_fs[i] = nuFnussc_fs[i] = nuFnuecd_fs[i] = 0.0;
		    nuFnublr_fs[i] = nuFnudt_fs[i] = 0.0;

		  }

		  distph_fs = (znht_fs * zn_fs) - (beta_fs * c * simtime);
		  if (distph_fs <= 0.0) {

		    distph_fs = 0.0;

		  }

		  if (accf == TRUE) {

		    timeobs_up_fs = (1.0 / D_Z) * (simtime + ((distph_fs * mu_prime) / c));
		    timeobs_upround_fs = myround(timeobs_up_fs);
		    sprintf(filename8, "%s%d%s", nfnfile, timeobs_upround_fs, ".dat");

		    if (fmod(write_count, 10) == 0) {

		      logMsg(logfile, "f", "timeobs_up_fs", f, time_count, timeobs_up_fs, filename8);
		    }

		    /*
		    sendToDetector(f - 1, timeobs_upround_fs, ndot_ph_up[f - 1], voldlshfsparam, TRUE, &nuFnu_max, sOutputDirectory, nfnfile, ndot_phsyn_up[f - 1], ndot_phsyndisord_up[f - 1], ndot_phssc_up[f - 1], ndot_phecd_up[f - 1], ndot_phblr_up[f - 1], ndot_phdt_up[f - 1], ndot_phup_up[f - 2], ndot_ph_dummy[f - 1], mu_prime, time_count);
		    */

		    sendToDetector(f - 1, timeobs_upround_fs, ndot_ph_up[f - 1], voldlshfsparam, TRUE, &nuFnu_max, sOutputDirectory, nfnfile, ndot_phsyn_up[f - 1], ndot_phssc_up[f - 1], ndot_phecd_up[f - 1], ndot_phblr_up[f - 1], ndot_phdt_up[f - 1], ndot_phup_up[f - 2], ndot_ph_dummy[f - 1], mu_prime, time_count);

		    /*
		    sendToDetector(f - 1, timeobs_upround_fs, ndot_phfeed_up[f - 1], voldlshfsparam, TRUE, &nuFnu_max, sOutputDirectory, nfnfile, ndot_phsyn_up[f - 1], ndot_phssc_up[f - 1], ndot_phecd_up[f - 1], ndot_phblr_up[f - 1], ndot_phdt_up[f - 1], ndot_phfeedup_up[f - 2], ndot_ph_dummy[f - 1], mu_prime, time_count);
		    */

		  }		    

		  nuFnu_max = 1.0;

		  for (zn = Nrs; zn < f; zn++) {

		    if ((zn == (f - 1)) && (accf == TRUE)) {

		      timeobs_si_fs = (1.0 / D_Z) * (simtime + ((distph_fs * mu_prime) / c));

		    } else {

		      timeobs_si_fs = (1.0 / D_Z) * (simtime + tdelsi_fs[zn]);

		    }

		    timeobs_siround_fs = myround(timeobs_si_fs);
		    sprintf(filename8, "%s%d%s", nfnfile, timeobs_siround_fs, ".dat");

		    if (fmod(write_count, 10) == 0) {

		      logMsg(logfile, "zn", "timeobs_si_fs", zn, time_count, timeobs_si_fs, filename8);
		    }

		    if ((accf == TRUE) && (zn == (f - 1))) {

		      if ((zn == (Ntot - 1)) || ((Ntot == 2) && (zn == 1))) {

			/*
                        sendToDetector2(zn, timeobs_siround_fs, ndot_ph_si[zn], ndot_ph_up[zn], voldlshfs_param, accf, &nuFnu_max, sOutputDirectory, nfnfile, zn_fs, ndot_phsyn_si[zn], ndot_phsyn_up[zn], ndot_phsyndisord_si[zn], ndot_phsyndisord_up[zn], ndot_phssc_si[zn], ndot_phssc_up[zn], ndot_phecd_si[zn], ndot_phecd_up[zn], ndot_phblr_si[zn], ndot_phblr_up[zn], ndot_phdt_si[zn], ndot_phdt_up[zn], ndot_phup_si[zn - 1], ndot_phup_up[zn - 1], ndot_ph_dummy[zn], ndot_ph_dummy[zn], time_count);
			*/

                        sendToDetector2(zn, timeobs_siround_fs, ndot_ph_si[zn], ndot_ph_up[zn], voldlshfs_param, accf, &nuFnu_max, sOutputDirectory, nfnfile, zn_fs, ndot_phsyn_si[zn], ndot_phsyn_up[zn], ndot_phssc_si[zn], ndot_phssc_up[zn], ndot_phecd_si[zn], ndot_phecd_up[zn], ndot_phblr_si[zn], ndot_phblr_up[zn], ndot_phdt_si[zn], ndot_phdt_up[zn], ndot_phup_si[zn - 1], ndot_phup_up[zn - 1], ndot_ph_dummy[zn], ndot_ph_dummy[zn], time_count);


			/*
			sendToDetector2(zn, timeobs_siround_fs, ndot_ph_si[zn], ndot_phfeed_up[zn], voldlshfs_param, accf, &nuFnu_max, sOutputDirectory, nfnfile, zn_fs, ndot_phsyn_si[zn], ndot_phsyn_up[zn], ndot_phssc_si[zn], ndot_phssc_up[zn], ndot_phecd_si[zn], ndot_phecd_up[zn], ndot_phblr_si[zn], ndot_phblr_up[zn], ndot_phdt_si[zn], ndot_phdt_up[zn], ndot_phfeedup_si[zn - 1], ndot_phfeedup_up[zn - 1], ndot_ph_dummy[zn], ndot_ph_dummy[zn], time_count);
			*/

		      } else {			 

			/*
			sendToDetector2(zn, timeobs_siround_fs, ndot_ph_si[zn], ndot_ph_up[zn], voldlshfs_param, accf, &nuFnu_max, sOutputDirectory, nfnfile, zn_fs, ndot_phsyn_si[zn], ndot_phsyn_up[zn], ndot_phsyndisord_si[zn], ndot_phsyndisord_up[zn], ndot_phssc_si[zn], ndot_phssc_up[zn], ndot_phecd_si[zn], ndot_phecd_up[zn], ndot_phblr_si[zn], ndot_phblr_up[zn], ndot_phdt_si[zn], ndot_phdt_up[zn], ndot_phup_si[zn - 1], ndot_phup_up[zn - 1], ndot_phdo_si[zn + 1], ndot_phdo_up[zn + 1], time_count);
			*/

                        sendToDetector2(zn, timeobs_siround_fs, ndot_ph_si[zn], ndot_ph_up[zn], voldlshfs_param, accf, &nuFnu_max, sOutputDirectory, nfnfile, zn_fs, ndot_phsyn_si[zn], ndot_phsyn_up[zn], ndot_phssc_si[zn], ndot_phssc_up[zn], ndot_phecd_si[zn], ndot_phecd_up[zn], ndot_phblr_si[zn], ndot_phblr_up[zn], ndot_phdt_si[zn], ndot_phdt_up[zn], ndot_phup_si[zn - 1], ndot_phup_up[zn - 1], ndot_phdo_si[zn + 1], ndot_phdo_up[zn + 1], time_count);

			/*
			sendToDetector2(zn, timeobs_siround_fs, ndot_ph_si[zn], ndot_phfeed_up[zn], voldlshfs_param, accf, &nuFnu_max, sOutputDirectory, nfnfile, zn_fs, ndot_phsyn_si[zn], ndot_phsyn_up[zn], ndot_phssc_si[zn], ndot_phssc_up[zn], ndot_phecd_si[zn], ndot_phecd_up[zn], ndot_phblr_si[zn], ndot_phblr_up[zn], ndot_phdt_si[zn], ndot_phdt_up[zn], ndot_phfeedup_si[zn - 1], ndot_phfeedup_up[zn - 1], ndot_phfeeddo_si[zn + 1], ndot_phfeeddo_up[zn + 1], time_count);
			*/
		      }

		    } else {
		     
		      if ((zn == (Ntot - 1)) || ((Ntot == 2) && (zn == 1))) {

			/*
                        sendToDetector2(zn, timeobs_siround_fs, ndot_ph_si[zn], ndot_ph_up[zn], voldlfs_param, accf, &nuFnu_max, sOutputDirectory, nfnfile, zn_fs, ndot_phsyn_si[zn], ndot_phsyn_up[zn], ndot_phsyndisord_si[zn], ndot_phsyndisord_up[zn], ndot_phssc_si[zn], ndot_phssc_up[zn], ndot_phecd_si[zn], ndot_phecd_up[zn], ndot_phblr_si[zn], ndot_phblr_up[zn], ndot_phdt_si[zn], ndot_phdt_up[zn], ndot_phup_si[zn - 1], ndot_phup_up[zn - 1], ndot_ph_dummy[zn], ndot_ph_dummy[zn], time_count);
			*/

                        sendToDetector2(zn, timeobs_siround_fs, ndot_ph_si[zn], ndot_ph_up[zn], voldlfs_param, accf, &nuFnu_max, sOutputDirectory, nfnfile, zn_fs, ndot_phsyn_si[zn], ndot_phsyn_up[zn], ndot_phssc_si[zn], ndot_phssc_up[zn], ndot_phecd_si[zn], ndot_phecd_up[zn], ndot_phblr_si[zn], ndot_phblr_up[zn], ndot_phdt_si[zn], ndot_phdt_up[zn], ndot_phup_si[zn - 1], ndot_phup_up[zn - 1], ndot_ph_dummy[zn], ndot_ph_dummy[zn], time_count);

			/*
			sendToDetector2(zn, timeobs_siround_fs, ndot_ph_si[zn], ndot_phfeed_up[zn], voldlfs_param, accf, &nuFnu_max, sOutputDirectory, nfnfile, zn_fs, ndot_phsyn_si[zn], ndot_phsyn_up[zn], ndot_phssc_si[zn], ndot_phssc_up[zn], ndot_phecd_si[zn], ndot_phecd_up[zn], ndot_phblr_si[zn], ndot_phblr_up[zn], ndot_phdt_si[zn], ndot_phdt_up[zn], ndot_phfeedup_si[zn - 1], ndot_phfeedup_up[zn - 1], ndot_ph_dummy[zn], ndot_ph_dummy[zn], time_count);
			*/

		      } else {

			/*
                        sendToDetector2(zn, timeobs_siround_fs, ndot_ph_si[zn], ndot_ph_up[zn], voldlfs_param, accf, &nuFnu_max, sOutputDirectory, nfnfile, zn_fs, ndot_phsyn_si[zn], ndot_phsyn_up[zn], ndot_phsyndisord_si[zn], ndot_phsyndisord_up[zn], ndot_phssc_si[zn], ndot_phssc_up[zn], ndot_phecd_si[zn], ndot_phecd_up[zn], ndot_phblr_si[zn], ndot_phblr_up[zn], ndot_phdt_si[zn], ndot_phdt_up[zn], ndot_phup_si[zn - 1], ndot_phup_up[zn - 1], ndot_phdo_si[zn + 1], ndot_phdo_up[zn + 1], time_count);
			*/

                        sendToDetector2(zn, timeobs_siround_fs, ndot_ph_si[zn], ndot_ph_up[zn], voldlfs_param, accf, &nuFnu_max, sOutputDirectory, nfnfile, zn_fs, ndot_phsyn_si[zn], ndot_phsyn_up[zn], ndot_phssc_si[zn], ndot_phssc_up[zn], ndot_phecd_si[zn], ndot_phecd_up[zn], ndot_phblr_si[zn], ndot_phblr_up[zn], ndot_phdt_si[zn], ndot_phdt_up[zn], ndot_phup_si[zn - 1], ndot_phup_up[zn - 1], ndot_phdo_si[zn + 1], ndot_phdo_up[zn + 1], time_count);

			/*
			sendToDetector2(zn, timeobs_siround_fs, ndot_ph_si[zn], ndot_phfeed_up[zn], voldlfs_param, accf, &nuFnu_max, sOutputDirectory, nfnfile, zn_fs, ndot_phsyn_si[zn], ndot_phsyn_up[zn], ndot_phssc_si[zn], ndot_phssc_up[zn], ndot_phecd_si[zn], ndot_phecd_up[zn], ndot_phblr_si[zn], ndot_phblr_up[zn], ndot_phdt_si[zn], ndot_phdt_up[zn], ndot_phfeedup_si[zn - 1], ndot_phfeedup_up[zn - 1], ndot_phfeeddo_si[zn + 1], ndot_phfeeddo_up[zn + 1], time_count);
			*/

		      }

		    }

		    for (j = 0; j < FREQ_GRID; j++) {
		      
		      if ((zn == (f - 1)) && (accf == TRUE)) {

			nuFnusyn_fs[j] += voldlshfs_param * SQR(syn_nu[j]) * ndot_syn_ep[zn][j];
			nuFnussc_fs[j] += voldlshfs_param * SQR(syn_nu[j]) * ndot_ssc[zn][j];
			nuFnuecd_fs[j] += voldlshfs_param * SQR(syn_nu[j]) * ndot_ecd[zn][j];
			nuFnublr_fs[j] += voldlshfs_param * SQR(syn_nu[j]) * ndot_blr[zn][j];
			nuFnudt_fs[j] += voldlshfs_param * SQR(syn_nu[j]) * ndot_dt[zn][j];

		      } else {

			nuFnusyn_fs[j] += voldlfs_param * SQR(syn_nu[j]) * ndot_syn_ep[zn][j];
			nuFnussc_fs[j] += voldlfs_param * SQR(syn_nu[j]) * ndot_ssc[zn][j];
			nuFnuecd_fs[j] += voldlfs_param * SQR(syn_nu[j]) * ndot_ecd[zn][j];
			nuFnublr_fs[j] += voldlfs_param * SQR(syn_nu[j]) * ndot_blr[zn][j];
			nuFnudt_fs[j] += voldlfs_param * SQR(syn_nu[j]) * ndot_dt[zn][j];
		      }

		    }

		  }

		  /*
		  sprintf(filename_nfnecd, "%s%e%s", nfnecd_file, z_c, ".dat");
		  fp_nfnecd = fopen(filename_nfnecd, "w");
		  if (!fp_nfnecd) {
		    printf("\nCould not open the output file component nuFnu. \n");
		    return -1;
		  }
		  */

		  nuFnuecd_fsmax = 1.0;
		  nuFnublr_fsmax = 1.0;
		  nuFnudt_fsmax = 1.0;

		  for (j = 0; j < FREQ_GRID; j++) {

		    /*
		    fprintf(fp_nfnecd, "%e %e %e %e %e %e\n", (D_Z * syn_nu[j]), max(1.0e-40, nuFnusyn_fs[j]),
			    max(1.0e-40, nuFnussc_fs[j]), max(1.0e-40, nuFnuecd_fs[j]), max(1.0e-40, nuFnublr_fs[j]),
			    max(1.0e-40, nuFnudt_fs[j]));
		    */

		    nuFnuecd_fsmax = max(nuFnuecd_fs[j], nuFnuecd_fsmax);
		    nuFnublr_fsmax = max(nuFnublr_fs[j], nuFnublr_fsmax);
		    nuFnudt_fsmax = max(nuFnudt_fs[j], nuFnudt_fsmax);

		  }

		  //fclose(fp_nfnecd);

		  if ((z_c >= raddist[RGRID - 1]) && (nuFnuecd_fsmax <= NUFNUMAXLIM)) {
		    LDISK_fs = 0;
		  }

                  if ((z_c >= Rout_BLR) && (nuFnublr_fsmax <= NUFNUMAXLIM)) {
                    LBLR_fs = 0;
                  }

                  if ((z_c > (2.0 * Rout_BLR)) && (nuFnudt_fsmax <= NUFNUMAXLIM)) {
                    LDT_fs = 0;
                  }


		  // Calculate the full spectrum in the observer's frame received from the reverse region
		  //
		  for (i = 0; i < FREQ_GRID; i++) {

		    nuFnusyn_rs[i] = nuFnussc_rs[i] = nuFnuecd_rs[i] = 0.0;
		    nuFnublr_rs[i] = nuFnudt_rs[i] = 0.0;
		  }

		  /*
		  sprintf(filename_nfnblr, "%s%e%s", nfnblr_file, z_c, ".dat");
		  fp_nfnblr = fopen(filename_nfnblr, "w");
		  if (!fp_nfnblr) {
		    printf("\nCould not open the output file component nuFnu. \n");
		    return -1;
		  }
		  */

		  for (zn = (Nrs - 1); zn >= r; zn--) {

		    if ((zn == r) && (accr == TRUE)) {

		      distph_rs = (znht_fs * zn_fs) + (beta_rs * c * simtime);
		      if (distph_rs >= ((znht_rs * zn_rs) + (znht_fs * zn_fs))) {

			distph_rs = (znht_rs * zn_rs) + (znht_fs * zn_fs);

		      }

		      timeobs_si_rs = (1.0 / D_Z) * (simtime + ((distph_rs * mu_prime) / c));

		    } else {

		      timeobs_si_rs = (1.0 / D_Z) * (simtime + tdelsi_rs[zn]);

		    }

		    timeobs_siround_rs = myround(timeobs_si_rs);
		    sprintf(filename8, "%s%d%s", nfnfile, timeobs_siround_rs, ".dat");

		    if (fmod(write_count, 10) == 0) {

		      logMsg(logfile, "zn", "timeobs_si_rs", zn, time_count, timeobs_si_rs, filename8);
		    }

		    if ((accr == TRUE) && (zn == r)) {
		      
		      if ((zn == 0) || ((Ntot == 2) && (zn == 0))) {

			/*
                        sendToDetector(zn, timeobs_siround_rs, ndot_ph_si[zn], voldlshrs_param, accr, &nuFnu_max, sOutputDirectory, nfnfile, ndot_phsyn_si[zn], ndot_phsyndisord_si[zn], ndot_phssc_si[zn], ndot_phecd_si[zn], ndot_phblr_si[zn], ndot_phdt_si[zn], ndot_ph_dummy[zn], ndot_phdo_si[zn + 1], mu_prime, time_count);
			*/

                        sendToDetector(zn, timeobs_siround_rs, ndot_ph_si[zn], voldlshrs_param, accr, &nuFnu_max, sOutputDirectory, nfnfile, ndot_phsyn_si[zn], ndot_phssc_si[zn], ndot_phecd_si[zn], ndot_phblr_si[zn], ndot_phdt_si[zn], ndot_ph_dummy[zn], ndot_phdo_si[zn + 1], mu_prime, time_count);

			/*
			sendToDetector(zn, timeobs_siround_rs, ndot_ph_si[zn], voldlshrs_param, accr, &nuFnu_max, sOutputDirectory, nfnfile, ndot_phsyn_si[zn], ndot_phssc_si[zn], ndot_phecd_si[zn], ndot_phblr_si[zn], ndot_phdt_si[zn], ndot_ph_dummy[zn], ndot_phfeeddo_si[zn + 1], mu_prime, time_count);
			*/

		      } else {

			/*
                        sendToDetector(zn, timeobs_siround_rs, ndot_ph_si[zn], voldlshrs_param, accr, &nuFnu_max, sOutputDirectory, nfnfile, ndot_phsyn_si[zn], ndot_phsyndisord_si[zn], ndot_phssc_si[zn], ndot_phecd_si[zn], ndot_phblr_si[zn], ndot_phdt_si[zn], ndot_phup_si[zn - 1], ndot_phdo_si[zn + 1], mu_prime, time_count);
			*/

                        sendToDetector(zn, timeobs_siround_rs, ndot_ph_si[zn], voldlshrs_param, accr, &nuFnu_max, sOutputDirectory, nfnfile, ndot_phsyn_si[zn], ndot_phssc_si[zn], ndot_phecd_si[zn], ndot_phblr_si[zn], ndot_phdt_si[zn], ndot_phup_si[zn - 1], ndot_phdo_si[zn + 1], mu_prime, time_count);

			/*
			sendToDetector(zn, timeobs_siround_rs, ndot_ph_si[zn], voldlshrs_param, accr, &nuFnu_max, sOutputDirectory, nfnfile, ndot_phsyn_si[zn], ndot_phssc_si[zn], ndot_phecd_si[zn], ndot_phblr_si[zn], ndot_phdt_si[zn], ndot_phfeedup_si[zn - 1], ndot_phfeeddo_si[zn + 1], mu_prime, time_count);
			*/

		      }			
				
		    } else {

		      if ((zn == 0) || ((Ntot == 2) && (zn == 0))) {

			/*
                        sendToDetector(zn, timeobs_siround_rs, ndot_ph_si[zn], voldlrs_param, accr, &nuFnu_max, sOutputDirectory, nfnfile, ndot_phsyn_si[zn], ndot_phsyndisord_si[zn], ndot_phssc_si[zn], ndot_phecd_si[zn], ndot_phblr_si[zn], ndot_phdt_si[zn], ndot_ph_dummy[zn], ndot_phdo_si[zn + 1], mu_prime, time_count);
			*/

                        sendToDetector(zn, timeobs_siround_rs, ndot_ph_si[zn], voldlrs_param, accr, &nuFnu_max, sOutputDirectory, nfnfile, ndot_phsyn_si[zn], ndot_phssc_si[zn], ndot_phecd_si[zn], ndot_phblr_si[zn], ndot_phdt_si[zn], ndot_ph_dummy[zn], ndot_phdo_si[zn + 1], mu_prime, time_count);

			/*
			sendToDetector(zn, timeobs_siround_rs, ndot_ph_si[zn], voldlrs_param, accr, &nuFnu_max, sOutputDirectory, nfnfile, ndot_phsyn_si[zn], ndot_phssc_si[zn], ndot_phecd_si[zn], ndot_phblr_si[zn], ndot_phdt_si[zn], ndot_ph_dummy[zn], ndot_phfeeddo_si[zn + 1], mu_prime, time_count);
			*/

		      } else {

			/*
                        sendToDetector(zn, timeobs_siround_rs, ndot_ph_si[zn], voldlrs_param, accr, &nuFnu_max, sOutputDirectory, nfnfile, ndot_phsyn_si[zn], ndot_phsyndisord_si[zn], ndot_phssc_si[zn], ndot_phecd_si[zn], ndot_phblr_si[zn], ndot_phdt_si[zn], ndot_phup_si[zn - 1], ndot_phdo_si[zn + 1], mu_prime, time_count);
			*/

                        sendToDetector(zn, timeobs_siround_rs, ndot_ph_si[zn], voldlrs_param, accr, &nuFnu_max, sOutputDirectory, nfnfile, ndot_phsyn_si[zn], ndot_phssc_si[zn], ndot_phecd_si[zn], ndot_phblr_si[zn], ndot_phdt_si[zn], ndot_phup_si[zn - 1], ndot_phdo_si[zn + 1], mu_prime, time_count);

			/*
			sendToDetector(zn, timeobs_siround_rs, ndot_ph_si[zn], voldlrs_param, accr, &nuFnu_max, sOutputDirectory, nfnfile, ndot_phsyn_si[zn], ndot_phssc_si[zn], ndot_phecd_si[zn], ndot_phblr_si[zn], ndot_phdt_si[zn], ndot_phfeedup_si[zn - 1], ndot_phfeeddo_si[zn + 1], mu_prime, time_count);
			*/
		      }

		    }

		    for (j = 0; j < FREQ_GRID; j++) {

		      if ((zn == r) && (accr == TRUE)) {
			
			nuFnusyn_rs[j] += voldlshrs_param * SQR(syn_nu[j]) * ndot_syn_ep[zn][j];
			nuFnussc_rs[j] += voldlshrs_param * SQR(syn_nu[j]) * ndot_ssc[zn][j];
			nuFnuecd_rs[j] += voldlshrs_param * SQR(syn_nu[j]) * ndot_ecd[zn][j];
			nuFnublr_rs[j] += voldlshrs_param * SQR(syn_nu[j]) * ndot_blr[zn][j];
			nuFnudt_rs[j] += voldlshrs_param * SQR(syn_nu[j]) * ndot_dt[zn][j];

		      } else {
			
			nuFnusyn_rs[j] += voldlrs_param * SQR(syn_nu[j]) * ndot_syn_ep[zn][j];
			nuFnussc_rs[j] += voldlrs_param * SQR(syn_nu[j]) * ndot_ssc[zn][j];
			nuFnuecd_rs[j] += voldlrs_param * SQR(syn_nu[j]) * ndot_ecd[zn][j];
			nuFnublr_rs[j] += voldlrs_param * SQR(syn_nu[j]) * ndot_blr[zn][j];
			nuFnudt_rs[j] += voldlrs_param * SQR(syn_nu[j]) * ndot_dt[zn][j];
		      }

		    }

		  }

                  nuFnuecd_rsmax = 1.0;
                  nuFnublr_rsmax = 1.0;
                  nuFnudt_rsmax = 1.0;

		  for (j = 0; j < FREQ_GRID; j++) {

		    /*
		    fprintf(fp_nfnblr, "%e %e %e %e %e %e\n", (D_Z * syn_nu[j]), max(1.0e-40, nuFnusyn_rs[j]),
			    max(1.0e-40, nuFnussc_rs[j]), max(1.0e-40, nuFnuecd_rs[j]),
			    max(1.0e-40, nuFnublr_rs[j]), max(1.0e-40, nuFnudt_rs[j]));
		    */

                    nuFnuecd_rsmax = max(nuFnuecd_rs[j], nuFnuecd_rsmax);
                    nuFnublr_rsmax = max(nuFnublr_rs[j], nuFnublr_rsmax);
                    nuFnudt_rsmax = max(nuFnudt_rs[j], nuFnudt_rsmax);

		  }

		  //fclose(fp_nfnblr);

                  if ((z_c >= raddist[RGRID - 1]) && (nuFnuecd_rsmax <= NUFNUMAXLIM)) {
                    LDISK_rs = 0;
		}

                  if ((z_c >= Rout_BLR) && (nuFnublr_rsmax <= NUFNUMAXLIM)) {
                    LBLR_rs = 0;
                  }

                  if ((z_c > (2.0 * Rout_BLR)) && (nuFnudt_rsmax <= NUFNUMAXLIM)) {
                    LDT_rs = 0;
                  }
		  
		}

		else {

		  // Calculate the full spectrum in the observer's frame recieved from the reverse region
		  // when mu < beta
		  //
		  //printf("In the mu < beta part for this simulation\n");

		  for (i = 0; i < FREQ_GRID; i++) {

		    nuFnu_rs_si[i] = nuFnu_rs_do[i] = nuFnu_fs_si[i] = 0.0;

		  }

		  distph_rs = (znht_rs * zn_rs) - (beta_rs * c * simtime);
		  if (distph_rs <= 0.0) {

		    distph_rs = 0.0;

		  }

		  if (accr == TRUE) {

		    timeobs_do_rs = (1.0 / D_Z) * (simtime + ((distph_rs * mu_prime_abs) / c));
		    timeobs_doround_rs = myround(timeobs_do_rs);
		    sprintf(filename8, "%s%d%s", nfnfile, timeobs_doround_rs, ".dat");

		    if (fmod(write_count, 10) == 0) {

		      logMsg(logfile, "r", "timeobs_do_rs", r, time_count, timeobs_do_rs, filename8);
		    }

		    /*
                    sendToDetector(r, timeobs_doround_rs, ndot_ph_do[r], voldlshrsparam, TRUE, &nuFnu_max, sOutputDirectory, nfnfile, ndot_phsyn_do[r], ndot_phsyndisord_do[r], ndot_phssc_do[r], ndot_phecd_do[r], ndot_phblr_do[r], ndot_phdt_do[r], ndot_ph_dummy[r], ndot_phdo_do[r + 1], mu_prime, time_count);
		    */

                    sendToDetector(r, timeobs_doround_rs, ndot_ph_do[r], voldlshrsparam, TRUE, &nuFnu_max, sOutputDirectory, nfnfile, ndot_phsyn_do[r], ndot_phssc_do[r], ndot_phecd_do[r], ndot_phblr_do[r], ndot_phdt_do[r], ndot_ph_dummy[r], ndot_phdo_do[r + 1], mu_prime, time_count);

		    /*
		    sendToDetector(r, timeobs_doround_rs, ndot_phfeed_do[r], voldlshrsparam, TRUE, &nuFnu_max, sOutputDirectory, nfnfile, ndot_phsyn_do[r], ndot_phssc_do[r], ndot_phecd_do[r], ndot_phblr_do[r], ndot_phdt_do[r], ndot_ph_dummy[r], ndot_phfeeddo_do[r + 1], mu_prime, time_count);
		    */
		  }

		  nuFnu_max = 1.0;

		  for (zn = (Nrs - 1); zn >= r; zn--) {

		    if ((zn == r) && (accr == TRUE)) {

		      timeobs_si_rs = (1.0 / D_Z) * (simtime + ((distph_rs * mu_prime_abs) / c));

		    } else {

		      timeobs_si_rs = (1.0 / D_Z) * (simtime + tdelsi_rs[zn]);

		    }

		    timeobs_siround_rs = myround(timeobs_si_rs);
		    sprintf(filename8, "%s%d%s", nfnfile, timeobs_siround_rs, ".dat");

		    if (fmod(write_count, 10) == 0) {

		      logMsg(logfile, "zn", "timeobs_si_rs", zn, time_count, timeobs_si_rs, filename8);
		    }

		    if ((accr == TRUE) && (zn == r)) {

		      if ((zn == 0) || ((Ntot == 2) && (zn == 0))) {

			/*
                        sendToDetector2(zn, timeobs_siround_rs, ndot_ph_si[zn], ndot_ph_do[zn], voldlshrs_param, accr, &nuFnu_max, sOutputDirectory, nfnfile, zn_rs, ndot_phsyn_si[zn], ndot_phsyn_do[zn], ndot_phsyndisord_si[zn], ndot_phsyndisord_do[zn], ndot_phssc_si[zn], ndot_phssc_do[zn], ndot_phecd_si[zn], ndot_phecd_do[zn], ndot_phblr_si[zn], ndot_phblr_do[zn], ndot_phdt_si[zn], ndot_phdt_do[zn], ndot_ph_dummy[zn], ndot_ph_dummy[zn], ndot_phdo_si[zn + 1], ndot_phdo_do[zn + 1], time_count);
			*/

                        sendToDetector2(zn, timeobs_siround_rs, ndot_ph_si[zn], ndot_ph_do[zn], voldlshrs_param, accr, &nuFnu_max, sOutputDirectory, nfnfile, zn_rs, ndot_phsyn_si[zn], ndot_phsyn_do[zn], ndot_phssc_si[zn], ndot_phssc_do[zn], ndot_phecd_si[zn], ndot_phecd_do[zn], ndot_phblr_si[zn], ndot_phblr_do[zn], ndot_phdt_si[zn], ndot_phdt_do[zn], ndot_ph_dummy[zn], ndot_ph_dummy[zn], ndot_phdo_si[zn + 1], ndot_phdo_do[zn + 1], time_count);

			/*
			sendToDetector2(zn, timeobs_siround_rs, ndot_ph_si[zn], ndot_phfeed_do[zn], voldlshrs_param, accr, &nuFnu_max, sOutputDirectory, nfnfile, zn_rs, ndot_phsyn_si[zn], ndot_phsyn_do[zn], ndot_phssc_si[zn], ndot_phssc_do[zn], ndot_phecd_si[zn], ndot_phecd_do[zn], ndot_phblr_si[zn], ndot_phblr_do[zn], ndot_phdt_si[zn], ndot_phdt_do[zn], ndot_ph_dummy[zn], ndot_ph_dummy[zn], ndot_phfeeddo_si[zn + 1], ndot_phfeeddo_do[zn + 1], time_count);
			*/

                      } else {

			/*
                        sendToDetector2(zn, timeobs_siround_rs, ndot_ph_si[zn], ndot_ph_do[zn], voldlshrs_param, accr, &nuFnu_max, sOutputDirectory, nfnfile, zn_rs, ndot_phsyn_si[zn], ndot_phsyn_do[zn], ndot_phsyndisord_si[zn], ndot_phsyndisord_do[zn], ndot_phssc_si[zn], ndot_phssc_do[zn], ndot_phecd_si[zn], ndot_phecd_do[zn], ndot_phblr_si[zn], ndot_phblr_do[zn], ndot_phdt_si[zn], ndot_phdt_do[zn], ndot_phup_si[zn - 1], ndot_phup_do[zn - 1], ndot_phdo_si[zn + 1], ndot_phdo_do[zn + 1], time_count);
			*/

                        sendToDetector2(zn, timeobs_siround_rs, ndot_ph_si[zn], ndot_ph_do[zn], voldlshrs_param, accr, &nuFnu_max, sOutputDirectory, nfnfile, zn_rs, ndot_phsyn_si[zn], ndot_phsyn_do[zn], ndot_phssc_si[zn], ndot_phssc_do[zn], ndot_phecd_si[zn], ndot_phecd_do[zn], ndot_phblr_si[zn], ndot_phblr_do[zn], ndot_phdt_si[zn], ndot_phdt_do[zn], ndot_phup_si[zn - 1], ndot_phup_do[zn - 1], ndot_phdo_si[zn + 1], ndot_phdo_do[zn + 1], time_count);

			/*
			sendToDetector2(zn, timeobs_siround_rs, ndot_ph_si[zn], ndot_phfeed_do[zn], voldlshrs_param, accr, &nuFnu_max, sOutputDirectory, nfnfile, zn_rs, ndot_phsyn_si[zn], ndot_phsyn_do[zn], ndot_phssc_si[zn], ndot_phssc_do[zn], ndot_phecd_si[zn], ndot_phecd_do[zn], ndot_phblr_si[zn], ndot_phblr_do[zn], ndot_phdt_si[zn], ndot_phdt_do[zn], ndot_phfeedup_si[zn - 1], ndot_phfeedup_do[zn - 1], ndot_phfeeddo_si[zn + 1], ndot_phfeeddo_do[zn + 1], time_count);
			*/
                      }

		    } else {

                      if ((zn == 0) || ((Ntot == 2) && (zn == 0))) {

			/*
                        sendToDetector2(zn, timeobs_siround_rs, ndot_ph_si[zn], ndot_ph_do[zn], voldlrs_param, accr, &nuFnu_max, sOutputDirectory, nfnfile, zn_rs, ndot_phsyn_si[zn], ndot_phsyn_do[zn], ndot_phsyndisord_si[zn], ndot_phsyndisord_do[zn], ndot_phssc_si[zn], ndot_phssc_do[zn], ndot_phecd_si[zn], ndot_phecd_do[zn], ndot_phblr_si[zn], ndot_phblr_do[zn], ndot_phdt_si[zn], ndot_phdt_do[zn], ndot_ph_dummy[zn], ndot_ph_dummy[zn], ndot_phdo_si[zn + 1], ndot_phdo_do[zn + 1], time_count);
			*/

                        sendToDetector2(zn, timeobs_siround_rs, ndot_ph_si[zn], ndot_ph_do[zn], voldlrs_param, accr, &nuFnu_max, sOutputDirectory, nfnfile, zn_rs, ndot_phsyn_si[zn], ndot_phsyn_do[zn], ndot_phssc_si[zn], ndot_phssc_do[zn], ndot_phecd_si[zn], ndot_phecd_do[zn], ndot_phblr_si[zn], ndot_phblr_do[zn], ndot_phdt_si[zn], ndot_phdt_do[zn], ndot_ph_dummy[zn], ndot_ph_dummy[zn], ndot_phdo_si[zn + 1], ndot_phdo_do[zn + 1], time_count);

			/*
			sendToDetector2(zn, timeobs_siround_rs, ndot_ph_si[zn], ndot_phfeed_do[zn], voldlrs_param, accr, &nuFnu_max, sOutputDirectory, nfnfile, zn_rs, ndot_phsyn_si[zn], ndot_phsyn_do[zn], ndot_phssc_si[zn], ndot_phssc_do[zn], ndot_phecd_si[zn], ndot_phecd_do[zn], ndot_phblr_si[zn], ndot_phblr_do[zn], ndot_phdt_si[zn], ndot_phdt_do[zn], ndot_ph_dummy[zn], ndot_ph_dummy[zn], ndot_phfeeddo_si[zn + 1], ndot_phfeeddo_do[zn + 1], time_count);
			*/

                      } else {

			/*
                        sendToDetector2(zn, timeobs_siround_rs, ndot_ph_si[zn], ndot_ph_do[zn], voldlrs_param, accr, &nuFnu_max, sOutputDirectory, nfnfile, zn_rs, ndot_phsyn_si[zn], ndot_phsyn_do[zn], ndot_phsyndisord_si[zn], ndot_phsyndisord_do[zn], ndot_phssc_si[zn], ndot_phssc_do[zn], ndot_phecd_si[zn], ndot_phecd_do[zn], ndot_phblr_si[zn], ndot_phblr_do[zn], ndot_phdt_si[zn], ndot_phdt_do[zn], ndot_phup_si[zn - 1], ndot_phup_do[zn - 1], ndot_phdo_si[zn + 1], ndot_phdo_do[zn + 1], time_count);
			*/

                        sendToDetector2(zn, timeobs_siround_rs, ndot_ph_si[zn], ndot_ph_do[zn], voldlrs_param, accr, &nuFnu_max, sOutputDirectory, nfnfile, zn_rs, ndot_phsyn_si[zn], ndot_phsyn_do[zn], ndot_phssc_si[zn], ndot_phssc_do[zn], ndot_phecd_si[zn], ndot_phecd_do[zn], ndot_phblr_si[zn], ndot_phblr_do[zn], ndot_phdt_si[zn], ndot_phdt_do[zn], ndot_phup_si[zn - 1], ndot_phup_do[zn - 1], ndot_phdo_si[zn + 1], ndot_phdo_do[zn + 1], time_count);

			/*
			sendToDetector2(zn, timeobs_siround_rs, ndot_ph_si[zn], ndot_phfeed_do[zn], voldlrs_param, accr, &nuFnu_max, sOutputDirectory, nfnfile, zn_rs, ndot_phsyn_si[zn], ndot_phsyn_do[zn], ndot_phssc_si[zn], ndot_phssc_do[zn], ndot_phecd_si[zn], ndot_phecd_do[zn], ndot_phblr_si[zn], ndot_phblr_do[zn], ndot_phdt_si[zn], ndot_phdt_do[zn], ndot_phfeedup_si[zn - 1], ndot_phfeedup_do[zn - 1], ndot_phfeeddo_si[zn + 1], ndot_phfeeddo_do[zn + 1], time_count);
			*/
                      }

		    }

		  }


		  // Calculate the full spectrum in the observer's frame received from the forward region
		  //
		  for (zn = Nrs; zn < f; zn++) {

		    if ((zn == (f - 1)) && (accf == TRUE)) {

		      distph_fs = (znht_rs * zn_rs) + (beta_fs * c * simtime);
		      if (distph_fs >= ((znht_rs * zn_rs) + (znht_fs * zn_fs))) {

			distph_fs = (znht_rs * zn_rs) + (znht_fs * zn_fs);

		      }

		      timeobs_si_fs = (1.0 / D_Z)
			* (simtime + ((distph_fs * mu_prime_abs) / c));

		    } else {

		      timeobs_si_fs = (1.0 / D_Z) * (simtime + tdelsi_fs[zn]);

		    }

		    timeobs_siround_fs = myround(timeobs_si_fs);
		    sprintf(filename8, "%s%d%s", nfnfile, timeobs_siround_fs, ".dat");

		    if (fmod(write_count, 10) == 0) {

		      logMsg(logfile, "zn", "timeobs_si_fs", zn, time_count, timeobs_si_fs, filename8);
		    }
		    
		    
		    if ((accf == TRUE) && (zn == (f - 1))) {
		      
		      if ((zn == (Ntot - 1)) || ((Ntot == 2) && (zn == 1))) {

			/*
                        sendToDetector(zn, timeobs_siround_fs, ndot_ph_si[zn], voldlshfs_param, accf, &nuFnu_max, sOutputDirectory, nfnfile, ndot_phsyn_si[zn], ndot_phsyndisord_si[zn], ndot_phssc_si[zn], ndot_phecd_si[zn], ndot_phblr_si[zn], ndot_phdt_si[zn], ndot_phup_si[zn - 1], ndot_ph_dummy[zn], mu_prime, time_count);
			*/

			sendToDetector(zn, timeobs_siround_fs, ndot_ph_si[zn], voldlshfs_param, accf, &nuFnu_max, sOutputDirectory, nfnfile, ndot_phsyn_si[zn], ndot_phssc_si[zn], ndot_phecd_si[zn], ndot_phblr_si[zn], ndot_phdt_si[zn], ndot_phup_si[zn - 1], ndot_ph_dummy[zn], mu_prime, time_count);

			/*
			sendToDetector(zn, timeobs_siround_fs, ndot_ph_si[zn], voldlshfs_param, accf, &nuFnu_max, sOutputDirectory, nfnfile, ndot_phsyn_si[zn], ndot_phssc_si[zn], ndot_phecd_si[zn], ndot_phblr_si[zn], ndot_phdt_si[zn], ndot_phfeedup_si[zn - 1], ndot_ph_dummy[zn], mu_prime, time_count);
			*/

		      } else {

			/*
                        sendToDetector(zn, timeobs_siround_fs, ndot_ph_si[zn], voldlshfs_param, accf, &nuFnu_max, sOutputDirectory, nfnfile, ndot_phsyn_si[zn], ndot_phsyndisord_si[zn], ndot_phssc_si[zn], ndot_phecd_si[zn], ndot_phblr_si[zn], ndot_phdt_si[zn], ndot_phup_si[zn - 1], ndot_phdo_si[zn + 1], mu_prime, time_count);
			*/

			sendToDetector(zn, timeobs_siround_fs, ndot_ph_si[zn], voldlshfs_param, accf, &nuFnu_max, sOutputDirectory, nfnfile, ndot_phsyn_si[zn], ndot_phssc_si[zn], ndot_phecd_si[zn], ndot_phblr_si[zn], ndot_phdt_si[zn], ndot_phup_si[zn - 1], ndot_phdo_si[zn + 1], mu_prime, time_count);

			/*
			sendToDetector(zn, timeobs_siround_fs, ndot_ph_si[zn], voldlshfs_param, accf, &nuFnu_max, sOutputDirectory, nfnfile, ndot_phsyn_si[zn], ndot_phssc_si[zn], ndot_phecd_si[zn], ndot_phblr_si[zn], ndot_phdt_si[zn], ndot_phfeedup_si[zn - 1], ndot_phfeeddo_si[zn + 1], mu_prime, time_count);
			*/

		      }
		      
		    } else {
		      
		      if ((zn == (Ntot - 1)) || ((Ntot == 2) && (zn == 1))) {

			/*
                        sendToDetector(zn, timeobs_siround_fs, ndot_ph_si[zn], voldlfs_param, accf, &nuFnu_max, sOutputDirectory, nfnfile, ndot_phsyn_si[zn], ndot_phsyndisord_si[zn], ndot_phssc_si[zn], ndot_phecd_si[zn], ndot_phblr_si[zn], ndot_phdt_si[zn], ndot_phup_si[zn - 1], ndot_ph_dummy[zn], mu_prime, time_count);
			*/

			sendToDetector(zn, timeobs_siround_fs, ndot_ph_si[zn], voldlfs_param, accf, &nuFnu_max, sOutputDirectory, nfnfile, ndot_phsyn_si[zn], ndot_phssc_si[zn], ndot_phecd_si[zn], ndot_phblr_si[zn], ndot_phdt_si[zn], ndot_phup_si[zn - 1], ndot_ph_dummy[zn], mu_prime, time_count);

			/*
			sendToDetector(zn, timeobs_siround_fs, ndot_ph_si[zn], voldlfs_param, accf, &nuFnu_max, sOutputDirectory, nfnfile, ndot_phsyn_si[zn], ndot_phssc_si[zn], ndot_phecd_si[zn], ndot_phblr_si[zn], ndot_phdt_si[zn], ndot_phfeedup_si[zn - 1], ndot_ph_dummy[zn], mu_prime, time_count);
			*/

		      } else {

			/*
                        sendToDetector(zn, timeobs_siround_fs, ndot_ph_si[zn], voldlfs_param, accf, &nuFnu_max, sOutputDirectory, nfnfile, ndot_phsyn_si[zn], ndot_phsyndisord_si[zn], ndot_phssc_si[zn], ndot_phecd_si[zn], ndot_phblr_si[zn], ndot_phdt_si[zn], ndot_phup_si[zn - 1], ndot_phdo_si[zn + 1], mu_prime, time_count);
			*/

                        sendToDetector(zn, timeobs_siround_fs, ndot_ph_si[zn], voldlfs_param, accf, &nuFnu_max, sOutputDirectory, nfnfile, ndot_phsyn_si[zn], ndot_phssc_si[zn], ndot_phecd_si[zn], ndot_phblr_si[zn], ndot_phdt_si[zn], ndot_phup_si[zn - 1], ndot_phdo_si[zn + 1], mu_prime, time_count);

			/*
			sendToDetector(zn, timeobs_siround_fs, ndot_ph_si[zn], voldlfs_param, accf, &nuFnu_max, sOutputDirectory, nfnfile, ndot_phsyn_si[zn], ndot_phssc_si[zn], ndot_phecd_si[zn], ndot_phblr_si[zn], ndot_phdt_si[zn], ndot_phfeedup_si[zn - 1], ndot_phfeeddo_si[zn + 1], mu_prime, time_count);
			*/
			
		      }

		    }

		  }

		}


		z_c += (Gamma_sh * Beta_sh * c * dt);

		if (fmod(write_count, 10) == 0) {
			fp_log = fopen(logfile, "a");
			if (!fp_log) {
				printf("\nCould not open the output file log. \n");
				return -1;
			}

			fprintf(fp_log, "time count = %d, z_c = %e, dt = %e\n", time_count,
					z_c, dt);

			fclose(fp_log);

		}

		time_count++;
		write_count++;
		timecount_fs++;
		timecount_rs++;
		Izindex++;

		if ((accf == FALSE) && (accr == FALSE)) {

			if (nuFnu_max <= NUFNUMAXLIM) {

		    FLUX = FALSE;
		    printf("The maximum nuFnu is %e, and hence flux_fs = %d\n",
			   nuFnu_max, FLUX);
				printf("Acceleration in both regions is %d & %d for time counter = %d\n",
			   accf, accr, time_count);
		    
		  }
		  
		}

		end = time(NULL);
		walltime = end - start;

		// IF PROGRAM IS GETTING CLOSE TO THE 24 HOUR MAX RUNTIME THEN SAVE THE FILES AND EXIT
		//
		if( (restrict_wall_time == 1) && (walltime >= (max_wall_time - file_write_time)) ) {

		  exitSimulation (start, logfile, 0, sOutputDirectory, nfnfile);
		}

	} // End of time loop.

	exitSimulation (start, logfile, 0, sOutputDirectory, nfnfile);

} // End of main program.
