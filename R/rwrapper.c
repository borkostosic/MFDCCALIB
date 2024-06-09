#include "..\mfdcca.h"

/********* wrapper function; only pointers can be passed from R **************/

void mfdcca(int *ver, double *data1, double *data2, int* total,
	double *HR, double* tauR, double* fR, double* alphaR,
	double* dmse, int *rs, int *nrs, double* perc,
	double *qmin, double *qmax, double *dq,
	int *minbox, int *maxbox, double *boxratio,
	int *nfit, int *sw)
{
	int i, nq, nr, integrate=1;
	double q, eps = 0.001;
	int dcca_version;

	cfg.npts = *total;
	cfg.minbox = *minbox ? *minbox : 4;
	cfg.maxbox = *maxbox ? *maxbox : *total / 4;
	cfg.boxratio = *boxratio ? *boxratio : pow(2.0, 1.0 / 8.0);
	cfg.nfit = *nfit ? *nfit : 1;		// nfit==1 linear fit, nfit==2 2nd degree poly, etc.
	cfg.sw = *sw ? 1 : 0;
	
	if (*nrs) {				// user supplied scale
		int tot, n;
		cfg.nr = n= *nrs;
		for (nr = 0; nr < cfg.nr; nr++) {
			cfg.rs[nr] = rs[nr];			// copy to cfg structure
		}
	}
	else
		rscale(&cfg);	// prepares scale

	setup(&cfg);		// fills in x values and powers for fitting
	
	dcca_version = *ver;
	prepare_data(data1, *total);
	if(data1!=data2)			// check if this is MFDCCA of MFDFA if data1 and data2 are the same
		prepare_data(data2, *total);
	calc_residues(data1, data2, *total, &cfg);
	*perc=calc_mfdcca(dcca_version, *qmin, *qmax, *dq);

//	mfdcca(&cfg, *total, data, integrate,
//		*qmin, *qmax, *dq, eps, H, tau, alpha, f);	//call the MFDFA algorithm

	if (!*nrs) {
		*nrs = cfg.nr;
		for (nr = 0; nr < cfg.nr; nr++) {
			rs[nr] = cfg.rs[nr];
		}
	}

	nq = 0;
	for (q = *qmin; q < *qmax; q += *dq) {
		for (nr = 0; nr < cfg.nr; nr++) {
			if(!isnan(logft[nq][nr]))
				dmse[nq * MAX_BOX + nr] = pow(10.0,logft[nq][nr]);
		}
		HR[nq] = H[nq];
		tauR[nq] = tau[nq];
		alphaR[nq] = alpha[nq];
		fR[nq] = f[nq];
		nq++;
	}



	return cfg.nr;
}

