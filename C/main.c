#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../mfdcca.h"

double data[2][MAX_DATA];

/*********************************************************************************/
int load(char* filename)
{
	int size = 0;
	double val1, val2;
	char cnom[2][256];
	FILE* fp;

	fp = fopen(filename, "r+");						/* open file */
	fscanf(fp, "%s\t%s\n", cnom[0], cnom[1]);
	while (fscanf(fp, "%lf\t%lf\n", &val1, &val2) == 2) {			/* loop until end of file */
		data[0][size] = val1;							/* store data */
		data[1][size] = val2;							/* store data */
		size++;
	}

	fclose(fp);										/* close file */

	return size;
}

/*********************************************************************************/

void main() {

  int i, w, m, n;
  double r;
  double q, qmin, qmax, dq;
  int window = 250, jump = 1, total, first, offset = 0, nlin=0;
  FILE *h;
  char fname[256], fnom[] = "BIOMIAL_MULTIFRACTAL_65536_0.3_0.4";
  //  char fname[256], fnom[]="sug_eta_returns";

  int nq, integrate;
  double eps = 0.001, perc;
  int dcca_version = MFDXA;
//MFDXA, ABS, MFCCA, PS, MS, PB, MB, PP, PM, MP, MM


  qmin = -10.000001; qmax = 10.0; dq = 0.1;

  sprintf(fname, "../data/%s.dat", fnom);
  total = load(fname);
  cfg.npts = total;

  printf("%d data read from %s\n",total,fname);

  cfg.minbox = 4;
  cfg.maxbox = total / 4;
  if (cfg.maxbox < cfg.minbox)
	  cfg.maxbox = cfg.minbox;

  cfg.boxratio = pow(2.0, 1.0 / 8.0);
  rscale(&cfg);		// prepares scale, allocates and fills in x values
  cfg.nfit = 1;		// nfit==1 means linear fit, nfit==2 means 2nd degree polynomial, and nfit==3 is 3rd degree
  integrate = 1;	// series repersnts icremets, integrate it to form a profile
  cfg.goback = 1;	// if non-overlapping windows, go backwards as well

  setup(&cfg);

  //cfg.sw = 1;	// sliding window (more statistics, but slower...)

  prepare_data(data[0], total);
  prepare_data(data[1], total);
  calc_residues(data[0], data[1], total, &cfg);
  perc=calc_mfdcca(dcca_version, qmin, qmax, dq);

  printf("used %.2f%% data pairs\n", perc);

  // print out the results to screen and file
  sprintf(fname, "%s_mfdfa%s.txt", fnom, cfg.sw ? "_SW":"");
  h = fopen(fname, "w+");
  printf("q\th(q)\ttau(q)\talpha\tf\n");
  fprintf(h, "q\th(q)\ttau(q)\talpha\tf\n");
  nq = 0;
  for (q = qmin; q < qmax; q += dq) {
	  printf("%f\t%f\t%f\t%f\t%f\n", q, H[nq], tau[nq], alpha[nq], f[nq]);
	  fprintf(h, "%f\t%f\t%f\t%f\t%f\n", q, H[nq], tau[nq], alpha[nq], f[nq]);
	  fflush(h);
	  nq++;
  }
  fclose(h);
}
