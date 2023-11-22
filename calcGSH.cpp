/*=================================================================
 * calcDT.c - Calculate the symmetrised Fourier Coefficients and their derivatives in matlab
 *
 *
 * Input:   phi1, Phi, phi2, CS, SS, maxL
 * Output:  F, L, DF, DDF
 *
 *
 *=================================================================*/


 //#define RDEBUG

 //structure that contains real and imaginary parts of a T and their derivatives
#include "mex.h"
#include "GSHFun.h"

// CASE 1 : mF = calcGSH(1, phi, alpha, CS, SS, maxL)
// 
// CASE 2 : [E, J, H] = calcGSH(2, phi, mF, CS, SS, maxL)
//         where phi is [phi1(1); Phi(1); phi2(1), alpha(1); phi1(2); Phi(2); phi2(2), alpha(2); ...]
// 
// CASE 3 : [E, J, H] = calcGSH(3, phi, alpha, mF, CS, SS, maxL)
//         where phi is [phi1(1); Phi(1); phi2(1); phi1(2); Phi(2); phi2(2); ...]
//         only varies phi
// 
// CASE 4 : [E, J, H] = calcGSH(4, phi, mF, CS, SS, maxL, K)
//         where phi is [phi1(1); Phi(1); phi2(1), alpha(1); phi1(2); Phi(2); phi2(2), alpha(2); ...]
//         E = E + Es*K
// 
// CASE 5 : E = calcGSH(5, phi, mF, CS, SS, maxL)
//         where phi is [phi1(1); Phi(1); phi2(1), alpha(1); phi1(2); Phi(2); phi2(2), alpha(2); ...]
//         only calculates E
// 
// CASE 6 : E = calcGSH(6, phi, mF, CS, SS, maxL, K)
//         where phi is [phi1(1); Phi(1); phi2(1), alpha(1); phi1(2); Phi(2); phi2(2), alpha(2); ...]
//         only calculates E
//         E = E + Es*K
//
// CASE 7 : [E J H] = calcGSH(7, phi, alpha, Fz, CS, SS, maxL, r, kernel)
//         where phi is [phi1(1); Phi(1); phi2(1), alpha(1); phi1(2); Phi(2); phi2(2), alpha(2); ...]
//         this is a part of the different texture calculator described in the paper
// 
// CASE 8 : [E] = calcGSH(8, phi, alpha, Fz, CS, SS, maxL, r, kernel)
//         where phi is [phi1(1); Phi(1); phi2(1), alpha(1); phi1(2); Phi(2); phi2(2), alpha(2); ...]
//         this is a part of the different texture calculator
//
// CASE 9 : F = calcGSH(9, phi, CS, SS, maxL, kernel)
//         where phi is [phi1(1); Phi(1); phi2(1); phi1(2); Phi(2); phi2(2); ...]
//         calculates F for each phi
// 
// CASE 10: [E J H] = calcGSH(10, phi, alpha, Fz, CS, SS, maxL, r, kernel)
//         where phi is [phi1(1); Phi(1); phi2(1); phi1(2); Phi(2); phi2(2); ...]
//         calculates error gradient hessian for maximum pole figure intensity of 6, 6 is hardcoded :)
// 
// CASE 11: mF = calcGSH(11, phi, alpha, CS, SS, maxL, kernel)
//          where phi is [phi1(1); Phi(1); phi2(1); phi1(2); Phi(2); phi2(2); ...]
//          calculates mean GSH coeffs but only on l where kernel[l] > 0
//          used as part of the different texture routine
// 
// CASE 12: FS = calcGSH(12, phi, alpha, CS, SS, maxL, kernel)
//          where phi is [phi1(1); Phi(1); phi2(1); phi1(2); Phi(2); phi2(2); ...]
//          for l = 0,maxL; FS[l] += |mF[l]|;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
   //output
   double* E, * J, * H;
   //input
   double* mF, * phi, * alpha, * kernel;
   size_t CS, SS, maxL, calcFlag;
   //quasi input
   size_t N, L, nvars, Z;
   double K, TDIS, * ES, r;

   calcFlag = (size_t)*mxGetPr(prhs[0]);

   if (mxGetNumberOfElements(prhs[0]) > 1) calcFlag = 0;
   //mexPrintf("calcFlag = %i\n", (int)calcFlag);
   switch (calcFlag) {
   case 0:
      //When you forget calcFlag
      phi = mxGetPr(prhs[0]);
      alpha = mxGetPr(prhs[1]);

      CS = (size_t)*mxGetPr(prhs[2]);
      SS = (size_t)*mxGetPr(prhs[3]);
      maxL = (size_t)*mxGetPr(prhs[4]);

      if (nrhs > 5) kernel = mxGetPr(prhs[5]);
      else {
         kernel = (double*)calloc(maxL + 1, sizeof(double));
         for (int i = 0; i < (maxL + 1); i++) kernel[i] = 1.0;
      }

      N = mxGetNumberOfElements(prhs[1]);
      L = calcTLEN(maxL, CS, SS);
      /* Create a 0-by-0 mxArray; you will allocate the memory dynamically */
      plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);

      mF = (double*)mxMalloc(L * sizeof(double));

      calcMF(mF, phi, alpha, kernel, CS, SS, maxL, N);

      mxSetPr(plhs[0], mF);
      mxSetM(plhs[0], L);
      mxSetN(plhs[0], 1);
      break;
   case 1:
      phi = mxGetPr(prhs[1]);
      alpha = mxGetPr(prhs[2]);

      CS = (size_t)*mxGetPr(prhs[3]);
      SS = (size_t)*mxGetPr(prhs[4]);
      maxL = (size_t)*mxGetPr(prhs[5]);

      if (nrhs > 6) kernel = mxGetPr(prhs[6]);
      else {
         kernel = (double*)calloc(maxL + 1, sizeof(double));
         for (int i = 0; i < (maxL + 1); i++) kernel[i] = 1.0;
      }

      N = mxGetNumberOfElements(prhs[2]);
      L = calcTLEN(maxL, CS, SS);
      /* Create a 0-by-0 mxArray; you will allocate the memory dynamically */
      plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);

      mF = (double*)mxMalloc(L * sizeof(double));

      calcMF(mF, phi, alpha, kernel, CS, SS, maxL, N);

      mxSetPr(plhs[0], mF);
      mxSetM(plhs[0], L);
      mxSetN(plhs[0], 1);

      if (!(nrhs > 6)) free(kernel);
      break;
   case 2:
      phi = mxGetPr(prhs[1]);
      mF = mxGetPr(prhs[2]);

      CS = (size_t)*mxGetPr(prhs[3]);
      SS = (size_t)*mxGetPr(prhs[4]);
      maxL = (size_t)*mxGetPr(prhs[5]);

      if (nrhs > 6) kernel = mxGetPr(prhs[6]);
      else {
         kernel = (double*)calloc(maxL + 1, sizeof(double));
         for (int i = 0; i <= maxL; i++) kernel[i] = 1.0;
      }

      nvars = mxGetNumberOfElements(prhs[1]);
      N = nvars / 4;

      /* Create a 0-by-0 mxArray; you will allocate the memory dynamically */
      plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
      plhs[1] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
      plhs[2] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);

      E = (double*)mxMalloc(sizeof(double));
      J = (double*)mxMalloc(nvars * sizeof(double));
      H = (double*)mxMalloc(nvars * nvars * sizeof(double));

      calcFULL(E, J, H, phi, mF, kernel, CS, SS, maxL, N);


      mxSetPr(plhs[0], E);
      mxSetM(plhs[0], 1);
      mxSetN(plhs[0], 1);

      mxSetPr(plhs[1], J);
      mxSetM(plhs[1], nvars);
      mxSetN(plhs[1], 1);

      mxSetPr(plhs[2], H);
      mxSetM(plhs[2], nvars);
      mxSetN(plhs[2], nvars);

      if (!(nrhs > 6)) free(kernel);
      break;
   case 3:
      phi = mxGetPr(prhs[1]);
      alpha = mxGetPr(prhs[2]);
      mF = mxGetPr(prhs[3]);

      CS = (size_t)*mxGetPr(prhs[4]);
      SS = (size_t)*mxGetPr(prhs[5]);
      maxL = (size_t)*mxGetPr(prhs[6]);

      if (nrhs > 7) kernel = mxGetPr(prhs[7]);
      else {
         kernel = (double*)calloc(maxL + 1, sizeof(double));
         for (int i = 0; i <= maxL; i++) kernel[i] = 1.0;
      }

      nvars = mxGetNumberOfElements(prhs[1]);
      N = nvars / 3;

      /* Create a 0-by-0 mxArray; you will allocate the memory dynamically */
      plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
      plhs[1] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
      plhs[2] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);

      E = (double*)mxMalloc(sizeof(double));
      J = (double*)mxMalloc(nvars * sizeof(double));
      H = (double*)mxMalloc(nvars * nvars * sizeof(double));

      calcPHIS(E, J, H, phi, alpha, mF, kernel, CS, SS, maxL, N);


      mxSetPr(plhs[0], E);
      mxSetM(plhs[0], 1);
      mxSetN(plhs[0], 1);

      mxSetPr(plhs[1], J);
      mxSetM(plhs[1], nvars);
      mxSetN(plhs[1], 1);

      mxSetPr(plhs[2], H);
      mxSetM(plhs[2], nvars);
      mxSetN(plhs[2], nvars);

      if (!(nrhs > 7)) free(kernel);
      break;
   case 4:
      phi = mxGetPr(prhs[1]);
      mF = mxGetPr(prhs[2]);

      CS = (size_t)*mxGetPr(prhs[3]);
      SS = (size_t)*mxGetPr(prhs[4]);
      maxL = (size_t)*mxGetPr(prhs[5]);
      K = (double)*mxGetPr(prhs[6]);

      if (nrhs > 7) kernel = mxGetPr(prhs[7]);
      else {
         kernel = (double*)calloc(maxL + 1, sizeof(double));
         for (int i = 0; i <= maxL; i++) kernel[i] = 1.0;
      }

      nvars = mxGetNumberOfElements(prhs[1]);
      N = nvars / 4;

      /* Create a 0-by-0 mxArray; you will allocate the memory dynamically */
      plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
      plhs[1] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
      plhs[2] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);

      if (nlhs > 3) plhs[3] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);

      E = (double*)mxMalloc(sizeof(double));
      J = (double*)mxMalloc(nvars * sizeof(double));
      H = (double*)mxMalloc(nvars * nvars * sizeof(double));

      if (nlhs > 3) ES = (double*)mxMalloc(sizeof(double));

      calcFULL(E, J, H, phi, mF, kernel, CS, SS, maxL, N);

      //TDIS = calcTDISCL(maxL, CS);
      ES[0] = E[0];
      //mexPrintf("   %16.8e\n", E[0]*TDIS);
      calcSTD(E, J, H, phi, N, K);


      mxSetPr(plhs[0], E);
      mxSetM(plhs[0], 1);
      mxSetN(plhs[0], 1);

      mxSetPr(plhs[1], J);
      mxSetM(plhs[1], nvars);
      mxSetN(plhs[1], 1);

      mxSetPr(plhs[2], H);
      mxSetM(plhs[2], nvars);
      mxSetN(plhs[2], nvars);
      if (nlhs > 3) {
         mxSetPr(plhs[3], ES);
         mxSetM(plhs[3], 1);
         mxSetN(plhs[3], 1);
      }
      if (!(nrhs > 7)) free(kernel);
      break;
   case 5:
      phi = mxGetPr(prhs[1]);
      mF = mxGetPr(prhs[2]);

      CS = (size_t)*mxGetPr(prhs[3]);
      SS = (size_t)*mxGetPr(prhs[4]);
      maxL = (size_t)*mxGetPr(prhs[5]);
      
      if (nrhs > 6) kernel = mxGetPr(prhs[6]);
      else {
         kernel = (double*)calloc(maxL + 1, sizeof(double));
         for (int i = 0; i <= maxL; i++) kernel[i] = 1.0;
      }


      nvars = mxGetNumberOfElements(prhs[1]);
      N = nvars / 4;

      /* Create a 0-by-0 mxArray; you will allocate the memory dynamically */
      plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);

      E = (double*)mxMalloc(sizeof(double));

      calcE(E, phi, mF, kernel, CS, SS, maxL, N);


      mxSetPr(plhs[0], E);
      mxSetM(plhs[0], 1);
      mxSetN(plhs[0], 1);
      
      if (!(nrhs > 6)) free(kernel);
      break;
   case 6:
      phi = mxGetPr(prhs[1]);
      mF = mxGetPr(prhs[2]);

      CS = (size_t)*mxGetPr(prhs[3]);
      SS = (size_t)*mxGetPr(prhs[4]);
      maxL = (size_t)*mxGetPr(prhs[5]);
      K = (double)*mxGetPr(prhs[6]);
      
      if (nrhs > 7) kernel = mxGetPr(prhs[7]);
      else {
         kernel = (double*)calloc(maxL + 1, sizeof(double));
         for (int i = 0; i <= maxL; i++) kernel[i] = 1.0;
      }

      nvars = mxGetNumberOfElements(prhs[1]);
      N = nvars / 4;

      /* Create a 0-by-0 mxArray; you will allocate the memory dynamically */
      plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);

      E = (double*)mxMalloc(sizeof(double));

      calcE(E, phi, mF, kernel, CS, SS, maxL, N);
      calcESTD(E, phi, N, K);


      mxSetPr(plhs[0], E);
      mxSetM(plhs[0], 1);
      mxSetN(plhs[0], 1);
      
      if (!(nrhs > 7)) free(kernel);
      break;
   case 7:
      phi = mxGetPr(prhs[1]);
      alpha = mxGetPr(prhs[2]);
      mF = mxGetPr(prhs[3]);

      CS = (size_t)*mxGetPr(prhs[4]);
      SS = (size_t)*mxGetPr(prhs[5]);
      maxL = (size_t)*mxGetPr(prhs[6]);
      r = (double)*mxGetPr(prhs[7]);

      if (nrhs > 8) kernel = mxGetPr(prhs[8]);
      else {
         kernel = (double*)calloc(maxL + 1, sizeof(double));
         for (int i = 0; i < (maxL + 1); i++) kernel[i] = 1.0;
      }

      if (maxL > RANK_MAX) maxL = RANK_MAX;
      L = calcTLENK(maxL, CS, SS, kernel);

      nvars = mxGetNumberOfElements(prhs[1]);
      Z = mxGetNumberOfElements(prhs[3]); // scale L later when we have it
      N = nvars / 3;
      Z = Z / L;

      /* Create a 0-by-0 mxArray; you will allocate the memory dynamically */
      plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
      plhs[1] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
      plhs[2] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);

      E = (double*)mxMalloc(sizeof(double));
      J = (double*)mxMalloc(nvars * sizeof(double));
      H = (double*)mxMalloc(nvars * nvars * sizeof(double));

      calcDIFFPHI(E, J, H, phi, alpha, mF, kernel, CS, SS, maxL, N, nvars, Z, r);


      mxSetPr(plhs[0], E);
      mxSetM(plhs[0], 1);
      mxSetN(plhs[0], 1);

      mxSetPr(plhs[1], J);
      mxSetM(plhs[1], nvars);
      mxSetN(plhs[1], 1);

      mxSetPr(plhs[2], H);
      mxSetM(plhs[2], nvars);
      mxSetN(plhs[2], nvars);

      if (!(nrhs > 8)) free(kernel);
      break;
   case 8:
      phi = mxGetPr(prhs[1]);
      alpha = mxGetPr(prhs[2]);
      mF = mxGetPr(prhs[3]);

      CS = (size_t)*mxGetPr(prhs[4]);
      SS = (size_t)*mxGetPr(prhs[5]);
      maxL = (size_t)*mxGetPr(prhs[6]);
      r = (double)*mxGetPr(prhs[7]);

      if (maxL > RANK_MAX) maxL = RANK_MAX;
      L = calcTLEN(maxL, CS, SS);

      nvars = mxGetNumberOfElements(prhs[1]);
      Z = mxGetNumberOfElements(prhs[3]);

      N = nvars / 3;
      Z = Z / L;



      /* Create a 0-by-0 mxArray; you will allocate the memory dynamically */
      plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
      plhs[1] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
      plhs[2] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);

      E = (double*)mxMalloc(sizeof(double));

      calcEDIFFPHI(E, phi, alpha, mF, CS, SS, maxL, N, nvars, Z, r);

      mxSetPr(plhs[0], E);
      mxSetM(plhs[0], 1);
      mxSetN(plhs[0], 1);
      break;
   case 9:
      phi = mxGetPr(prhs[1]);

      CS = (size_t)*mxGetPr(prhs[2]);
      SS = (size_t)*mxGetPr(prhs[3]);
      maxL = (size_t)*mxGetPr(prhs[4]);

      if (nrhs > 5) kernel = mxGetPr(prhs[5]);
      else {
         kernel = (double*)calloc(maxL + 1, sizeof(double));
         for (int i = 0; i < (maxL + 1); i++) kernel[i] = 1.0;
      }

      N = mxGetNumberOfElements(prhs[1]);
      L = calcTLEN(maxL, CS, SS);
      N = N / 3;

      /* Create a 0-by-0 mxArray; you will allocate the memory dynamically */
      plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);

      mF = (double*)mxMalloc(L * N * sizeof(double));

      calcF(mF, phi, kernel, CS, SS, maxL, N);

      mxSetPr(plhs[0], mF);
      mxSetM(plhs[0], L);
      mxSetN(plhs[0], N);

      if (!(nrhs > 5)) free(kernel);
      break;
   case 10:
      phi = mxGetPr(prhs[1]);
      alpha = mxGetPr(prhs[2]);
      mF = mxGetPr(prhs[3]);

      CS = (size_t)*mxGetPr(prhs[4]);
      SS = (size_t)*mxGetPr(prhs[5]);
      maxL = (size_t)*mxGetPr(prhs[6]);
      r = (double)*mxGetPr(prhs[7]);

      if (nrhs > 8) kernel = mxGetPr(prhs[8]);
      else {
         kernel = (double*)calloc(maxL + 1, sizeof(double));
         for (int i = 0; i < (maxL + 1); i++) kernel[i] = 1.0;
      }

      if (maxL > RANK_MAX) maxL = RANK_MAX;
      L = calcTLEN(maxL, CS, SS);

      nvars = mxGetNumberOfElements(prhs[1]);
      Z = mxGetNumberOfElements(prhs[3]); // scale L later when we have it
      N = nvars / 3;
      Z = Z / L;

      if (nlhs == 3) {

         /* Create a 0-by-0 mxArray; you will allocate the memory dynamically */
         plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
         plhs[1] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
         plhs[2] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);

         E = (double*)mxMalloc(sizeof(double));
         J = (double*)mxMalloc(nvars * sizeof(double));
         H = (double*)mxMalloc(nvars * nvars * sizeof(double));

         calcMAXPHI(E, J, H, phi, alpha, mF, kernel, CS, SS, maxL, N, nvars, Z, r);


         mxSetPr(plhs[0], E);
         mxSetM(plhs[0], 1);
         mxSetN(plhs[0], 1);

         mxSetPr(plhs[1], J);
         mxSetM(plhs[1], nvars);
         mxSetN(plhs[1], 1);

         mxSetPr(plhs[2], H);
         mxSetM(plhs[2], nvars);
         mxSetN(plhs[2], nvars);
      }
      else {
         /* Create a 0-by-0 mxArray; you will allocate the memory dynamically */
         plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);

         E = (double*)mxMalloc(sizeof(double));

         calcEMAXPHI(E, phi, alpha, mF, kernel, CS, SS, maxL, N, nvars, Z, r);


         mxSetPr(plhs[0], E);
         mxSetM(plhs[0], 1);
         mxSetN(plhs[0], 1);
      }

      if (!(nrhs > 8)) free(kernel);
      break;
    case 11:
      phi = mxGetPr(prhs[1]);
      alpha = mxGetPr(prhs[2]);

      CS = (size_t)*mxGetPr(prhs[3]);
      SS = (size_t)*mxGetPr(prhs[4]);
      maxL = (size_t)*mxGetPr(prhs[5]);
      kernel = mxGetPr(prhs[6]);

      N = mxGetNumberOfElements(prhs[2]);
      L = calcTLENK(maxL, CS, SS, kernel);
      /* Create a 0-by-0 mxArray; you will allocate the memory dynamically */
      plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);

      mF = (double*)mxMalloc(L * sizeof(double));

      calcMFK(mF, phi, alpha, kernel, CS, SS, maxL, N);

      mxSetPr(plhs[0], mF);
      mxSetM(plhs[0], L);
      mxSetN(plhs[0], 1);
      break;
   case 12:
      phi = mxGetPr(prhs[1]);
      alpha = mxGetPr(prhs[2]);

      CS = (size_t)*mxGetPr(prhs[3]);
      SS = (size_t)*mxGetPr(prhs[4]);
      maxL = (size_t)*mxGetPr(prhs[5]);
      

      N = mxGetNumberOfElements(prhs[2]);
      L = calcTLENK(maxL, CS, SS, kernel);
      /* Create a 0-by-0 mxArray; you will allocate the memory dynamically */
      plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);

      mF = (double*)mxMalloc((maxL+1) * sizeof(double));

      calcMFR(mF, phi, alpha, CS, SS, maxL, N);

      mxSetPr(plhs[0], mF);
      mxSetM(plhs[0], (maxL+1));
      mxSetN(plhs[0], 1);
      break;
      case 13:
         phi = mxGetPr(prhs[1]);
         mF = mxGetPr(prhs[2]);

         CS = (size_t)*mxGetPr(prhs[3]);
         SS = (size_t)*mxGetPr(prhs[4]);
         maxL = (size_t)*mxGetPr(prhs[5]);

         if (nrhs > 6) kernel = mxGetPr(prhs[6]);
         else {
            kernel = (double*)calloc(maxL + 1, sizeof(double));
            for (int i = 0; i <= maxL; i++) kernel[i] = 1.0;
         }

         nvars = mxGetNumberOfElements(prhs[1]);
         N = nvars / 4;

         /* Create a 0-by-0 mxArray; you will allocate the memory dynamically */
         plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
         plhs[1] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
         plhs[2] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);

         E = (double*)mxMalloc(sizeof(double));
         J = (double*)mxMalloc(nvars * sizeof(double));
         H = (double*)mxMalloc(nvars * nvars * sizeof(double));

         calcFULL2(E, J, H, phi, mF, kernel, CS, SS, maxL, N);


         mxSetPr(plhs[0], E);
         mxSetM(plhs[0], 1);
         mxSetN(plhs[0], 1);

         mxSetPr(plhs[1], J);
         mxSetM(plhs[1], nvars);
         mxSetN(plhs[1], 1);

         mxSetPr(plhs[2], H);
         mxSetM(plhs[2], nvars);
         mxSetN(plhs[2], nvars);

         if (!(nrhs > 6)) free(kernel);
         break;
         
         case 14:
         phi = mxGetPr(prhs[1]);
      alpha = mxGetPr(prhs[2]);
      mF = mxGetPr(prhs[3]);

      CS = (size_t)*mxGetPr(prhs[4]);
      SS = (size_t)*mxGetPr(prhs[5]);
      maxL = (size_t)*mxGetPr(prhs[6]);

      if (nrhs > 7) kernel = mxGetPr(prhs[7]);
      else {
         kernel = (double*)calloc(maxL + 1, sizeof(double));
         for (int i = 0; i <= maxL; i++) kernel[i] = 1.0;
      }

      nvars = mxGetNumberOfElements(prhs[1]);
      N = nvars / 3;

      /* Create a 0-by-0 mxArray; you will allocate the memory dynamically */
      plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
      plhs[1] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
      plhs[2] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);

      E = (double*)mxMalloc(sizeof(double));
      J = (double*)mxMalloc(nvars * sizeof(double));
      H = (double*)mxMalloc(nvars * nvars * sizeof(double));

      calcPHIS2(E, J, H, phi, alpha, mF, kernel, CS, SS, maxL, N);


      mxSetPr(plhs[0], E);
      mxSetM(plhs[0], 1);
      mxSetN(plhs[0], 1);

      mxSetPr(plhs[1], J);
      mxSetM(plhs[1], nvars);
      mxSetN(plhs[1], 1);

      mxSetPr(plhs[2], H);
      mxSetM(plhs[2], nvars);
      mxSetN(plhs[2], nvars);

      if (!(nrhs > 7)) free(kernel);
      break;
   }




   return;
}