#include "math.h"
#include <omp.h>
#include "blas.h"
#include <stdlib.h>
#include "mex.h"
#include "GSHLibrary.h"

#define RANK_MAX  64
#define sign(i) 		(((i)%2) ? -1.0:1.0)
#define fsign(i) 		(((i) < 0) ? -1.0:1.0)

#define sq2pi		2.506628274631000502415765284811
#define sq2			1.4142135623730950488016887242097
#define TWOPI		6.283185307179586476925286766559
#define PI			3.14159265358979323846264338327950

// All of these are the same as the CS.id in mtex so we can use the mtex symmetries directly
#define SYM_CUBIC          45
#define SYM_HEXAGONAL      40
#define SYM_ORTHORHOMBIC   16
#define SYM_TRICLINIC       2
#define SYM_TETRAGONAL     32
#define SYM_TRIGONAL       24


#define kahansum(LIM, FUN, SUM) {double kahanc, kahany, kahant; kahanc = 0.0;SUM = 0.0;for (size_t k = 0; k < LIM; k++) {kahany = FUN - kahanc;kahant = SUM + kahany;kahanc = (kahant - SUM) - kahany;SUM = kahant;}}

//#define RDEBUG

//structure that contains real and imaginary parts of a T and their derivatives
typedef struct {
   double re, im;
   double dre[3];
   double dim[3];
   double ddre[6];
   double ddim[6];
} gshe_complex;

// Function Header Definitions

//complex math
gshe_complex gshe_init(size_t DF);
gshe_complex cx_add(size_t DF, gshe_complex x, gshe_complex y);
gshe_complex cx_sub(size_t DF, gshe_complex x, gshe_complex y);
gshe_complex cx_conj(size_t DF, gshe_complex z);
gshe_complex sc_add(size_t DF, gshe_complex x, double y);
gshe_complex sc_mult(size_t DF, gshe_complex x, double y);

//index of databases
//size_t almnsp_ind(size_t l, size_t m, size_t n);


//Symmetry Assignments
size_t getSymCombo(size_t CS, size_t SS);
size_t sLinEq(size_t iSymCombo, size_t symm, size_t l);
void NumLinEq(size_t iSymCombo, size_t symm, size_t maxL, size_t* NumLin);

//mu -> m & nu -> n
double normFactorN(size_t symm, size_t l, size_t nu, int* n);
double normFactorM(size_t symm, size_t mu, int* m);
int TriMap(size_t nu);

// Assigns Calculation of T
gshe_complex symmetricT(size_t DF, size_t iSymCombo, size_t CS, size_t SS, size_t l, size_t mu, size_t nu, double sphi1[], double cphi1[], double sPhi[], double cPhi[], double sphi2[], double cphi2[]);

// Combines base calculations of T to form symmetrised coefficients
gshe_complex tricTricT(size_t DF, size_t l, size_t mu, size_t nu, double sphi1[], double cphi1[], double sPhi[], double cPhi[], double sphi2[], double cphi2[]);
gshe_complex lowerTricT(size_t DF, size_t CS, size_t l, size_t mu, size_t nu, double sphi1[], double cphi1[], double sPhi[], double cPhi[], double sphi2[], double cphi2[]);
gshe_complex lowLowerT(size_t DF, size_t CS, size_t SS, size_t l, size_t mu, size_t nu, double sphi1[], double cphi1[], double sPhi[], double cPhi[], double sphi2[], double cphi2[]);
gshe_complex cubLowerT(size_t DF, size_t SS, size_t l, size_t mu, size_t nu, double sphi1[], double cphi1[], double sPhi[], double cPhi[], double sphi2[], double cphi2[]);
gshe_complex cubCubicT(size_t DF, size_t l, size_t mu, size_t nu, double sphi1[], double cphi1[], double sPhi[], double cPhi[], double sphi2[], double cphi2[]);
gshe_complex cubTricT(size_t DF, size_t l, size_t mu, size_t nu, double sphi1[], double cphi1[], double sPhi[], double cPhi[], double sphi2[], double cphi2[]);

// Basic Calculations of T and its combinations
gshe_complex noSymT(size_t DF, size_t l, int m, int n, double sphi1[], double cphi1[], double sPhi[], double cPhi[], double sphi2[], double cphi2[]);
gshe_complex Slmn(size_t DF, size_t l, int m, int n, double sphi1[], double cphi1[], double sPhi[], double cPhi[], double sphi2[], double cphi2[]);
gshe_complex Rlmn(size_t DF, size_t l, int m, int n, double sphi1[], double cphi1[], double sPhi[], double cPhi[], double sphi2[], double cphi2[]);

size_t calcTLEN(size_t maxL, size_t CS, size_t SS);
size_t calcTLENK(size_t maxL, size_t CS, size_t SS, double* kernel);
double calcTDISCL(size_t maxL, size_t CS);

void calcF(double* mF, double* phi, double* kernel, size_t CS, size_t SS, size_t maxL, size_t N);
void calcMF(double* mF, double* phi, double* alpha, double* kernel, size_t CS, size_t SS, size_t maxL, size_t N);
void calcMFK(double* mF, double* phi, double* alpha, double* kernel, size_t CS, size_t SS, size_t maxL, size_t N);
void calcMFR(double* mF, double* phi, double* alpha, size_t CS, size_t SS, size_t maxL, size_t N);

void calcPHIS(double* E, double* J, double* H, double* phi, double* alpha, double* mF, double* kernel, size_t CS, size_t SS, size_t maxL, size_t N);
void calcPHIS2(double* E, double* J, double* H, double* phi, double* alpha, double* mF, double* kernel, size_t CS, size_t SS, size_t maxL, size_t N);
void calcFULL(double* E, double* J, double* H, double* phi, double* mF, double* kernel, size_t CS, size_t SS, size_t maxL, size_t N);
void calcFULL2(double* E, double* J, double* H, double* phi, double* mF, double* kernel, size_t CS, size_t SS, size_t maxL, size_t N);
void calcE(double* E, double* phi, double* mF, double* kernel, size_t CS, size_t SS, size_t maxL, size_t N);

void calcSTD(double* E, double* J, double* H, double* phi, size_t N, double K);
void calcESTD(double* E, double* phi, size_t N, double K);

void calcDIFFPHI(double* E, double* J, double* H, double* phi, double* alpha, double* mF, double* kernel, const size_t CS, const size_t SS, const size_t maxL, const size_t N, const size_t nvars, const size_t Z, double r);
void calcEDIFFPHI(double* E, double* phi, double* alpha, double* mF, const size_t CS, const size_t SS, size_t const maxL, size_t const N, const size_t nvars, const size_t Z, double r);

void calcMAXPHI(double* E, double* J, double* H, double* phi, double* alpha, double* F, double* kernel, const size_t CS, const size_t SS, size_t const maxL, size_t const N, const size_t nvars, const size_t Z, double r);
void calcEMAXPHI(double* E, double* phi, double* alpha, double* F, double* kernel, const size_t CS, const size_t SS, size_t const maxL, size_t const N, const size_t nvars, const size_t Z, double r);
