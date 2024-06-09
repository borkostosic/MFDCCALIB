#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_DATA	1048576	// maximum data points 2^20
#define MAX_BOX	200			// maximum 200 points on logarithmic scale...
#define MAXQ	201			// max q resolution -10,...,10 dq=0.1
#define MAX_FIT	3
// different MFDCCA versions
#define MFDXA	1	// Original MF-DXA W.-X. Zhou, Phys. Rev. E 77, 066211 (2008).
#define ABS		2	// MF-DXA with absolute values of fluctuation products
#define MFCCA	3	// MFCCA Phys. Rev. E 89, 023305
#define PS		4	// Plus sum
#define MS		5	// Minus sum
#define PB		6	// Plus Box
#define MB		7	// Minus box
#define PP		8	// Plus plus
#define PM		9	// Plus minus
#define MP		10	// Minus plus
#define MM		11	// Minus minus

double mse[MAX_BOX];
double logst[MAX_BOX], logft[MAXQ][MAX_BOX], logft2[MAXQ][MAX_BOX];
double H[MAXQ], H2[MAXQ], tau[MAXQ], alpha[MAXQ], f[MAXQ];

typedef struct {			// configurateion structure, holds all parameters
	int npts;
	int minbox;				// minimum box size
	int maxbox;				// maximum box size
	double boxratio;		// multiplicative factor for box size
	int rs[MAX_BOX];		// box size array 
	double x[MAX_FIT * 2][MAX_DATA];					// absicssa for fitting
	int nfit;				// order of the regression fit
	int nr;					// number of box sizes 
	int sw;					// sliding window flag
	int goback;				// go backwards
}DFA_CONFIG;

DFA_CONFIG cfg = { 0 };		// configuration structure
double *res1, *res2;	// residues

// function declarations
void setup(void);
int rscale(DFA_CONFIG* cfg);
double polyfit(double x[][MAX_DATA], double* y, int boxsize, int nfit);
double fit_poly1(double x[][MAX_DATA], double* y, int n, double* a, double* chisq);
double fit_poly2(double x[][MAX_DATA], double* y, int n, double* a, double* chisq);
double fit_poly3(double x[][MAX_DATA], double* y, int n, double* a, double* chisq);
double do_calc_slope(double* x, double* y, int nn);

/**********************************************************************************/

double calc_mse(int ver, double q)
{
int i, j, k, inc, boxsize, sgn;
int totsize = 0, count, nb;
double *s1, *s2, temp, msetemp;
int totbox = 0, totpass = 0;

	s1 = res1; s2 = res2;

	for (i = 0; i < cfg.nr; i++) {
		count = 0;
		mse[i] = 0;
		boxsize = cfg.rs[i];
		inc = cfg.sw ? 1 : boxsize;
		for (j = 0; j <= cfg.npts - boxsize; j += inc) {
			totsize += boxsize;
			msetemp = nb = 0;
			for (k = 0; k < boxsize; k++) {
				totpass++;
				switch (ver) {
				case MFDXA:					// Original MF-DXA W.-X. Zhou, Phys. Rev. E 77, 066211 (2008).
				case MFCCA:					// MFCCA Phys. Rev. E 89, 023305
				case PS:					// Plus sum
				case MS:					// Minus sum
					temp = s1[k] * s2[k];
					msetemp += temp;
					nb++;
					break;
				case PB:					// Plus Box
					temp = s1[k] * s2[k];
					if (temp > 0) {
						msetemp += temp;
						nb++;
					}
					break;
				case MB:					// Minus box
					temp = s1[k] * s2[k];
					if (temp < 0) {
						msetemp += -temp;
						nb++;
					}
					break;
				case PP:					// PP
					if ((s1[k]>0)&&(s2[k]>0)) {
						temp = s1[k] * s2[k];
						msetemp += temp;
						nb++;
					}
					break;
				case PM:					// PM
					if ((s1[k] > 0) && (s2[k] < 0)) {
						temp = s1[k] * s2[k];
						msetemp += -temp;
						nb++;
					}
					break;
				case MP:					// MP
					if ((s1[k] < 0) && (s2[k] > 0)) {
						temp = s1[k] * s2[k];
						msetemp += -temp;
						nb++;
					}
					break;
				case MM:					// MM
					if ((s1[k] < 0) && (s2[k] < 0)) {
						temp = s1[k] * s2[k];
						msetemp += temp;
						nb++;
					}
					break;
				case ABS:					// ABS
					temp = fabs(s1[k] * s2[k]);
					msetemp += temp;
					nb++;
					break;
				default:
					break;
				}
			}
			msetemp /= nb ? nb : 1; 

			switch (ver) {
			case PB:				// Plus Box
			case MB:				// Minus box
			case PP:				// PP
			case PM:				// PM
			case MP:				// MP
			case MM:				// MM
				if (!msetemp)
					break;
			case MFDXA:				// Original MF-DXA W.-X. Zhou, Phys. Rev. E 77, 066211 (2008).
			case ABS:				// ABS
					mse[i] += q ? pow(msetemp, q / 2) : log(msetemp) / 2;
					count++;
					totbox += nb;
				break;
			case PS:				// Plus sum
				if (msetemp > 0) {
					mse[i] += q ? pow(msetemp, q / 2) : log(msetemp) / 2;
					count++;
					totbox += nb;
				}
				break;
			case MS:				// Minus sum
				if (msetemp < 0) {
					msetemp = -msetemp;
					mse[i] += q ? pow(msetemp, q / 2) : log(msetemp) / 2;
					count++;
					totbox += nb;
				}
				break;
			case MFCCA:				// MFCCA Phys. Rev. E 89, 023305
				sgn = msetemp >= 0 ? 1 : -1;
				msetemp = fabs(msetemp);
				mse[i] += q ? sgn*pow(msetemp, q / 2) : sgn*log(msetemp) / 2;
				count++;
				totbox += nb;
				break;
			default:
				break;
			}
			s1 += boxsize;
			s2 += boxsize;
		}
		mse[i] /= count ? count : 1;
	}
	return(totbox * 100.0 / totpass);
}

/***********************************************************************/
void prepare_data(double *data, int total)
{
int i, c;
double average;

average = 0;
for (i = 0; i < total; i++)
	average += data[i];
average /= total;
data[0] -= average;
for (i = 1; i < total; i++)
	data[i] = data[i-1] + (data[i] - average);
}

/********** find residues ********************************/
void calc_residues(double *data1, double *data2, long npts, DFA_CONFIG *cfg)
{
	long i, boxsize, inc, j, k, bsize, np, nm;
	double stat, mseorig, msetemp, mtp, mtm, temp;
	int totsize = 0;
	double *s1, *s2;
	int nfit = cfg->nfit;
	int* rs = cfg->rs;
	int nr = cfg->nr;
	int sw = cfg->sw;

	for (i = 0; i < nr; i++) {
		boxsize = rs[i];
		inc = sw ? 1 : boxsize;
		for (j = 0; j <= npts - boxsize; j += inc) {
			totsize += boxsize;
		}
	}

	if (res1)
		free(res1);
	s1 = res1 = (double *)malloc(totsize * sizeof(double));
	if (res2)
		free(res2);
	s2 = res2 = (double*)malloc(totsize * sizeof(double));

	for (i = 0; i < nr; i++) {
		boxsize = rs[i];
		bsize = boxsize * sizeof(double);
		inc = sw ? 1 : boxsize;
		if (i == 41)
			i = i;
		for (j = 0; j <= npts - boxsize; j += inc) {
			memcpy(s1, &data1[j], bsize);
			polyfit(cfg->x, s1, boxsize, nfit);
			s1+= boxsize;
			memcpy(s2, &data2[j], bsize);
			polyfit(cfg->x, s2, boxsize, nfit);
			s2+= boxsize;
		}
	}
}

int rslen;	/* length of rs[] */

/************************************************************************************************
modified from C.K.Peng's original code to use the DFA_CONFIG structure and 0 offset for rs array:
rscale() allocates and fills rs[], the array of box sizes used by mfdfa()
below.  The box sizes range from (exactly) minbox to (approximately) maxbox,
and are arranged in a geometric series such that the ratio between
consecutive box sizes is (approximately) boxratio.  The return value is
the number of box sizes in rs[].
************************************************************************************************/
int rscale(DFA_CONFIG* cfg)
{
	int i, j, ir, n, tot;
	long rw, rslen;

	/* Determine how many scales are needed. */
	rslen = (int)(log10(cfg->maxbox / (double)cfg->minbox) / log10(cfg->boxratio) + 1.5);
	/* Thanks to Peter Domitrovich for pointing out that a previous version
	of the above calculation undercounted the number of scales in some situations. */
	for (ir = 1, n = 1, cfg->rs[0] = cfg->minbox; n <= rslen && cfg->rs[n - 1] < cfg->maxbox; ir++)
		if ((rw = cfg->minbox * pow(cfg->boxratio, ir) + 0.5) > cfg->rs[n - 1])
			cfg->rs[n++] = rw;
	if (cfg->rs[n - 1] > cfg->maxbox) --n;
	cfg->nr = n;
	tot = cfg->rs[n - 1];

	return (n);
}

/*************** prepare for fitting *****************************/
void setup(DFA_CONFIG* cfg)
{
	int i, j; 
	for (i = 0; i < MAX_DATA; i++) {
		cfg->x[0][i] = i + 1;
		cfg->x[1][i] = cfg->x[0][i] * cfg->x[0][i];
		for (j = 2; j < MAX_FIT * 2; j++)
			cfg->x[j][i] = cfg->x[j - 1][i] * cfg->x[0][i];
	}
}

/*************** wrapper function for polynomial fit *****************************/
double polyfit(double x[][MAX_DATA], double* y, int boxsize, int nfit)
{
	double mse, chisq, a[5];
	switch (nfit) {
	case 1:
		mse = fit_poly1(x, y, boxsize, a, &chisq);	// linear fit
		break;
	case 2:
		mse = fit_poly2(x, y, boxsize, a, &chisq);	// second order polynomil
		break;
	case 3:
		mse = fit_poly3(x, y, boxsize, a, &chisq);	// third order polynomial
		break;
	default:
		return -1.0;
	}
	return mse;
}

/*************** linear fit *************************************/
double fit_poly1(double x[][MAX_DATA], double* y, int n, double* a, double* chisq)
{
	double sx = 0, sx2 = 0, sy = 0, sxy = 0, temp, mse;
	int i;

	for (i = 0; i < n; i++) {
		sx += x[0][i];
		sy += y[i];
		sx2 += x[1][i];
		sxy += x[0][i] * y[i];
	}
	a[1] = (n * sxy - sx * sy) / (n * sx2 - sx * sx);
	a[0] = (sy - (a[1]) * sx) / n;

	*chisq = 0;
	for (i = 0; i < n; i++) {
		temp = a[1] * x[0][i] + a[0] - y[i];
		y[i] = temp;
		*chisq += temp * temp;
	}
	mse = *chisq;
	*chisq /= n;
	return mse;
}

/************** 2nd order polynomial fit ************************************/
double fit_poly2(double x[][MAX_DATA], double* y, int n, double* a, double* chisq)
{
	double sx, sy, sx2, sxy, sx3, sx4, syx2, temp, xx, mse;
	int i;

	sx = sy = sx2 = sxy = sx3 = sx4 = syx2 = 0;
	for (i = 0; i < n; i++) {
		sx += x[0][i];
		sy += y[i];
		sx2 += x[1][i];
		sxy += x[0][i] * y[i];
		syx2 += y[i] * x[1][i];
		sx3 += x[2][i];
		sx4 += x[3][i];
	}

	temp = (n * sx3 * sx3 - n * sx2 * sx4 + sx2 * sx2 * sx2 + sx * sx * sx4 - 2.0 * sx3 * sx * sx2);
	if (temp == 0.0)
		return 0;
	a[1] = (-sx * sx2 * syx2 + sx * sy * sx4 + sx2 * sx2 * sxy + n * sx3 * syx2 - sx2 * sx3 * sy - n * sxy * sx4) / temp;
	a[0] = -(-sy * sx3 * sx3 + sy * sx2 * sx4 + sx * sx3 * syx2 - sx * sxy * sx4 - sx2 * sx2 * syx2 + sx2 * sx3 * sxy) / temp;
	a[2] = -(n * sx2 * syx2 + sx * sx2 * sxy - sx2 * sx2 * sy - n * sx3 * sxy - sx * sx * syx2 + sx3 * sx * sy) / temp;

	*chisq = 0;
	for (i = 0; i < n; i++) {
		temp = a[2] * x[1][i] + a[1] * x[0][i] + a[0] - y[i];
		y[i] = temp;
		*chisq += temp * temp;
	}
	mse = *chisq;
	*chisq /= n;
	return mse;
}

/************** 3rd order polynomial fit ************************************/
double fit_poly3(double x[][MAX_DATA], double* y, int n, double* a, double* chisq)
{
	double sx = 0, sy = 0, sx2 = 0, syx = 0, sx3 = 0, sx4 = 0, syx2 = 0, temp = 0, xx = 0, mse;
	double sx5 = 0, sx6 = 0, syx3 = 0, xxx = 0, xxxx = 0;
	int i;

	for (i = 0; i < n; i++) {
		sx += x[0][i];
		sy += y[i];
		sx2 += x[1][i];
		sx3 += x[2][i];
		sx4 += x[3][i];
		sx5 += x[4][i];
		sx6 += x[5][i];
		syx += x[0][i] * y[i];
		syx2 += y[i] * x[1][i];
		syx3 += y[i] * x[2][i];
	}

	temp = (sx4 * n * sx2 * sx6 + sx * sx * sx5 * sx5 - n * sx3 *
		sx3 * sx6 + 2.0 * sx5 * sx2 * sx2 * sx3 + 2.0 * sx * sx3 * sx2 * sx6 + 2.0 * sx5 * n * sx3 * sx4 - sx2 * sx2 * sx2 * sx6 - sx * sx * sx4 * sx6 + sx2 * sx2
		* sx4 * sx4 - 2.0 * sx * sx5 * sx3 * sx3 - sx4 * sx4 * sx4 * n - 3.0 * sx3 * sx3 * sx2 * sx4 - sx5 * sx5 * n * sx2 - 2.0 * sx * sx5 * sx2 * sx4 + 2.0 *
		sx * sx4 * sx4 * sx3 + sx3 * sx3 * sx3 * sx3);

	a[0] = (-sy * sx4 * sx4 * sx4 - sx2 * sx4 * sx5 * syx + sx2 * sx3 * syx * sx6 - sx3 * sx3 * sx4 * syx2 - sx5 * sx3 * sx3 * syx + sx3 * sx3 * sx3 * syx3 + 2.0 * sy * sx5 * sx3 * sx4 - sx * syx * sx4 * sx6 + sx * sx3 * syx2 * sx6 - sx * sx4 * sx5 * syx2 + sx * syx * sx5 * sx5 + sx2 * sx4 * sx4 *
		syx2 + sx3 * sx4 * sx4 * syx - sy * sx3 * sx3 * sx6 + sx * sx4 * sx4 * syx3 - sx * sx3 * sx5 * syx3 - sy * sx5 * sx5 * sx2 - sx2 * sx2 * syx2 * sx6 - 2.0 *
		sx3 * sx4 * sx2 * syx3 + sx2 * sx2 * sx5 * syx3 + sy * sx4 * sx2 * sx6 + sx3 * sx5 * sx2 * syx2) / temp;

	a[1] = (sx * sy * sx5 * sx5 + sx * sx2 * syx2 * sx6 - sx * sx5 * sx2 * syx3 + sx * sx4 * sx3 * syx3 - sx * sx3 * sx5 * syx2 - sx * sy * sx4 * sx6 + sx3 * sx4 * sx4 * sy - n * sx4 * sx4 * syx3 - sx2 * sx2 * syx * sx6 - sx2 * sx3 * sx3 * syx3 + n * sx3 * sx5 * syx3 - sx2 * sx4 * sx5 * sy + n * syx * sx4
		* sx6 + 2.0 * sx2 * sx5 * sx3 * syx - n * sx3 * syx2 * sx6 + n * sx4 * sx5 * syx2 + sx2 * sx2 * sx4 * syx3 - sx3 * sx3 * sx5 * sy - n * syx * sx5 * sx5 -
		sx3 * sx3 * sx4 * syx + sx3 * sx3 * sx3 * syx2 - sx2 * sx3 * sx4 * syx2 + sx2 * sx3 * sy * sx6) / temp;

	a[2] = (-sx * sx * syx2 * sx6 + sx * sx * sx5 * syx3 - sx * sx4 * sx2 * syx3 - sx * sx4 * sx5 * sy + 2.0 * sx * sx3 * sx4 * syx2 - sx * sx3 * sx3 * syx3 + sx * sx2 * syx * sx6 - sx * sx5 * sx3 * syx + sx * sx3 * sy * sx6 - sx4 * sx4 * n * syx2 + sx4 * sx4 * sx2 * sy - sx4 * sx2 * sx3 * syx - sx4 * sx3
		* sx3 * sy + sx4 * n * sx3 * syx3 + sx4 * sx5 * n * syx + n * sx2 * syx2 * sx6 - sx2 * sx2 * sy * sx6 + sx5 * sx3 * sx2 * sy - sx3 * sx3 * sx2 * syx2 + sx2 *
		sx2 * sx3 * syx3 - n * sx3 * syx * sx6 + sx3 * sx3 * sx3 * syx - sx5 * n * sx2 * syx3) / temp;

	a[3] = -1 / (sx4 * n * sx2 * sx6 + sx * sx * sx5 * sx5 - n * sx3 * sx3 * sx6 + 2.0 * sx5 * sx2 * sx2 * sx3 + 2.0 * sx * sx3 * sx2 * sx6 + 2.0 * sx5 * n * sx3 * sx4 - sx2 * sx2 * sx2 * sx6 - sx * sx * sx4 * sx6 + sx2 * sx2 * sx4 * sx4 - 2.0 * sx * sx5 * sx3 * sx3 - sx4 * sx4 * sx4 * n - 3.0 *
		sx3 * sx3 * sx2 * sx4 - sx5 * sx5 * n * sx2 - 2.0 * sx * sx5 * sx2 * sx4 + 2.0 * sx * sx4 * sx4 * sx3 + sx3 * sx3 * sx3 * sx3) * (-sx * sx * sx5 * syx2 +
			n * sx3 * sx3 * syx3 - n * sx3 * sx4 * syx2 + sx4 * sx4 * n * syx - sx * sx4 * sx4 * sy - sx4 * n * sx2 * syx3 + sx * sx3 * sx3 * syx2 + sx2 * sx2 * sx2 *
			syx3 - sx3 * sx3 * sx3 * sy - sx5 * sx2 * sx2 * sy + sx5 * n * sx2 * syx2 + sx3 * sx3 * sx2 * syx - sx2 * sx2 * sx3 * syx2 - sx2 * sx2 * sx4 * syx - sx * sx4
			* sx3 * syx + sx * sx * sx4 * syx3 + sx * sx5 * sx2 * syx + 2.0 * sx3 * sx2 * sx4 * sy - 2.0 * sx * sx3 * sx2 * syx3 - sx5 * n * sx3 * syx + sx * sx5 * sx3 *
			sy + sx * sx2 * sx4 * syx2);


	*chisq = 0;
	for (i = 0; i < n; i++) {
		temp = a[3] * x[2][i] + a[2] * x[1][i] + a[1] * x[0][i] + a[0] - y[i];
		y[i] = temp;
		*chisq += temp * temp;
	}
	mse = *chisq;
	*chisq /= n;
	return mse;
}

/******************************************************************************************/

double do_calc_slope(double* x, double* y, int nn)
{
	double sx = 0, sx2 = 0, sy = 0, sxy = 0, temp, sstot = 0, ymean = 0, slope;
	int i, n;

	n = 0;
	for (i = 0; i < nn; i++) {
		if (isnan(x[i]))
			continue;
		if (isnan(y[i]))
			continue;
		if (!isfinite(y[i]))
			continue;
		sx += x[i];
		sy += y[i];
		sx2 += x[i] * x[i];
		sxy += x[i] * y[i];
		n++;
	}
	slope = (n * sxy - sx * sy) / (n * sx2 - sx * sx);
	return slope;
}

/******************************************************************/

double calc_mfdcca(int dcca_version, double minq, double maxq, double dq) {	// MFDFA
	int i, nq;
	double q;
	char string[256];
	double eps = 0.001, perc=0;

	for (i = 0; i < cfg.nr; i++)
		logst[i] = log10((double)cfg.rs[i]);

	nq = 0;
	for (q = minq; q <= maxq; q += dq) {
		perc+=calc_mse(dcca_version, q);
		for (i = 0; i < cfg.nr; i++){
			if (q)
				logft[nq][i] = log10(mse[i]) / q;
			else
				logft[nq][i] = log10(exp(mse[i]));
		}

		calc_mse(dcca_version, q + eps);
		for (i = 0; i < cfg.nr; i++) {
			if (q + eps)
				logft2[nq][i] = log10(mse[i]) / (q + eps);
			else
				logft2[nq][i] = log10(exp(mse[i]));
		}
		nq++;
	}

	for (i = 0; i < nq; i++) {
		H[i] = do_calc_slope(logst, logft[i], cfg.nr);
		H2[i] = do_calc_slope(logst, logft2[i], cfg.nr);
		q = minq + i * dq;
		tau[i] = q * H[i] - 1.0;
		alpha[i] = H[i] + q * (H2[i] - H[i]) / eps;
		f[i] = q * alpha[i] - tau[i];
	}
	return(perc / nq);
}	//end MFDCCA