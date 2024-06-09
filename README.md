
# MFDCCALIB
## Introduction
This code repository presents a library implementing a set of novel algorithms proposed to deal with negative cross-covariance among two series 
in the scope of Multifractal detrended cross-correlation analysis
[[1]](https://doi.org/xxx),
together with
three existing algorithms for MFDCCA found in the literature.

## Setup
The library is written in C, contained in the ```mfdcca.h``` header with examples of integration in ```C```, ```R```, and ```Python```, found in the respective directories. Below are instructions for running examples in each environment.
* ```C```: Run the Microsoft Visual Studio project.
* ```Python```: Compile the dynamic library  ```python setup.py build```, copy it over ```cp build/lib/mfdcca*```, run the example ```python mfdcca.py``` (precompiled Windows dynamic link library ```mfdcca.cp311-win_amd64.pyd``` is also provided)
* ```R```: Run the example ```mfdcca.R``` using the precompiled Windows libraries (Visual Studio project and code for compiling the dlls is also provided for portability).

## Library
The library exposes a single function that computes the multifractal spectrum. The API of that function is defined as follows and has similar arguments for the ```R``` and ```Python``` wrappers.
```
double calc_mfdcca(int dcca_version, double minq, double maxq, double dq);
```
where different MFDCCA versions are defined as (see [[1]](https://doi.org/xxx))

```
#define MFDXA	1	// Original MF-DXA W.-X. Zhou, Phys. Rev. E 77, 066211 (2008).
#define ABS	2	// MF-DXA with absolute values of fluctuation products
#define MFCCA	3	// MFCCA Phys. Rev. E 89, 023305
#define PS	4	// Plus sum
#define MS	5	// Minus sum
#define PB	6	// Plus Box
#define MB	7	// Minus box
#define PP	8	// Plus plus
#define PM	9	// Plus minus
#define MP	10	// Minus plus
#define MM	11	// Minus minus
```
and  ```qmin, qmax, dq``` represent the scaling parameter range

The configuration structure defined in ```mfdccaa.h``` header as:
```
typedef struct {			// configurateion structure, holds all parameters
	int npts;			// number of data point pairs
	int minbox;			// minimum box size
	int maxbox;			// maximum box size
	double boxratio;		// multiplicative factor for box size
	int rs[MAX_BOX];		// box size array 
	double x[MAX_FIT * 2][MAX_DATA];// absicssa for fitting
	int nfit;			// order of the regression fit
	int nr;				// number of box sizes 
	int sw;				// sliding window flag
	int goback;			// go backwards
}DFA_CONFIG;
```
contains the necessary parameters and is filled in by the wrapper code for C, R and Python together with the input data, while the output arrays
```
double H[MAXQ], tau[MAXQ], alpha[MAXQ], f[MAXQ];
```
are exposed as globals in the ```mfdccaa.h``` header.

## Results
An example of application of the original MF-DXA algorithm for the Binomial multifractal model [[2]](https://doi.org/10.1016/S0378-4371(02)01383-3), is shown below for two sequences of 2^20=1048576 numbers each, for p=0.3 and p=0.4, included under ```data/``` 

<img width="" alt="" src="./Fig1c.png">

The red curves correspond to the average of the MFDFA theoretical curves for the two series.

## Citation
If you use this work in academic research, citating the following reference would be appreciated:

```
@software{borkostosic2024MFDFA,
  author = {Stosic, Borko},
  title = {Multifractal detrended fluctuation analysis software},
  url = {https://github.com/borkostosic/mfdfa},
  version = {1.0.0},
  year = {2024},
}
```

## Contact
Borko Stosic (borkostosic@gmail.com)

 
