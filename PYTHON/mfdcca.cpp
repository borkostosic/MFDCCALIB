#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <math.h>
#include "../mfdcca.h"

namespace py = pybind11;

double calc(py::array_t<double> data1, py::array_t<double> data2, 
	int total, int ver, 
	py::array_t<double> HP, py::array_t<double> tauP,  
	py::array_t<double> fP, py::array_t<double> alphaP,
	py::array_t<double> dmseP,
	py::array_t<int> rs, int nrs,
	double qmin, double qmax, double dq, 
	int minbox, int maxbox, double boxratio, int nfit, int sw) {
    double* data1_ptr = static_cast<double*>(data1.request().ptr);
    double* data2_ptr = static_cast<double*>(data2.request().ptr);
    double* H_ptr = static_cast<double*>(HP.request().ptr);
    double* tau_ptr = static_cast<double*>(tauP.request().ptr);
    double* f_ptr = static_cast<double*>(fP.request().ptr);
    double* alpha_ptr = static_cast<double*>(alphaP.request().ptr);
    int* rs_ptr = static_cast<int*>(rs.request().ptr);
	double* dmse_ptr = static_cast<double*>(dmseP.request().ptr);
    double eps = 0.001;
    int integrate = 1;	// series repersents icremets, integrate it to form a profile
    int nq, nr;
    double q, perc;
	int dcca_version = ver;

	cfg.npts = total;
	cfg.minbox = minbox ? minbox : 4;
	cfg.maxbox = maxbox ? maxbox : total / 4;
	cfg.boxratio = boxratio ? boxratio : pow(2.0, 1.0 / 8.0);
	cfg.nfit = nfit ? nfit : 1;		// nfit==1 linear fit, nfit==2 2nd degree poly, etc.
	cfg.sw = sw ? 1 : 0;

	if (nrs) {				// user supplied scale
		int tot, n, i;
		cfg.nr = n = nrs;
		for (nr = 0; nr < cfg.nr; nr++) {
			cfg.rs[nr] = rs_ptr[nr];			// copy to cfg structure
		}
		tot = cfg.rs[n - 1];
	}
	else
		rscale(&cfg);	// prepares scale, allocates and fills in x value

	setup(&cfg);		// fills in x values and powers for fitting

	prepare_data(data1_ptr, total);
	if(data1_ptr!=data2_ptr)	// check if this is MFDCCA or MFDFA if two series are the same
		prepare_data(data2_ptr, total);
	calc_residues(data1_ptr, data2_ptr, total, &cfg);
	perc = calc_mfdcca(dcca_version, qmin, qmax, dq);
	
	for (nr = 0; nr<cfg.nr; nr++) {
		rs_ptr[nr] = cfg.rs[nr];
	}
    
	nq = 0;
    for (q = qmin; q < qmax; q += dq) {
	  for(nr=0; nr<cfg.nr; nr++){
		  if (!isnan(logft[nq][nr]))
			  dmse_ptr[nq * MAX_BOX + nr] = pow(10.0, logft[nq][nr]);
		  else
			  dmse_ptr[nq * MAX_BOX + nr] = 0;
	  }
	  H_ptr[nq] = H[nq];
	  tau_ptr[nq] = tau[nq];
	  alpha_ptr[nq] = alpha[nq];
	  f_ptr[nq] = f[nq];
	  nq++;
    }

    return perc;

}

PYBIND11_MODULE(TORCH_EXTENSION_NAME, m)
{
    m.def("calc", &calc, "calc");
}