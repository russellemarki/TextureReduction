#include "GSHFun.h"

// Functions
gshe_complex gshe_init(size_t DF) {
   // Initializes Tlmn to zero
   gshe_complex z;
   z.re = 0.0;
   z.im = 0.0;
   if (DF > 0) {
      z.dre[0] = 0.0;
      z.dre[1] = 0.0;
      z.dre[2] = 0.0;

      z.dim[0] = 0.0;
      z.dim[1] = 0.0;
      z.dim[2] = 0.0;

      if (DF > 1) {
         for (size_t i = 0; i < 6; i++) {
            z.ddre[i] = 0.0;
            z.ddim[i] = 0.0;
         }
      }
   }
   return z;
}

gshe_complex cx_add(size_t DF, gshe_complex x, gshe_complex y) {
   // adds two complex numbers together
   gshe_complex z;
   z.re = x.re + y.re; z.im = x.im + y.im;
   if (DF > 0) {
      z.dre[0] = x.dre[0] + y.dre[0];
      z.dre[1] = x.dre[1] + y.dre[1];
      z.dre[2] = x.dre[2] + y.dre[2];

      z.dim[0] = x.dim[0] + y.dim[0];
      z.dim[1] = x.dim[1] + y.dim[1];
      z.dim[2] = x.dim[2] + y.dim[2];
      if (DF > 1) {
         for (size_t i = 0; i < 6; i++) {
            z.ddre[i] = x.ddre[i] + y.ddre[i];
            z.ddim[i] = x.ddim[i] + y.ddim[i];
         }
      }
   }
   return z;
}

gshe_complex cx_sub(size_t DF, gshe_complex x, gshe_complex y) {
   //subtracts complex x from complex y
   gshe_complex z;
   z.re = x.re - y.re; z.im = x.im - y.im;
   if (DF > 0) {
      z.dre[0] = x.dre[0] - y.dre[0];
      z.dre[1] = x.dre[1] - y.dre[1];
      z.dre[2] = x.dre[2] - y.dre[2];

      z.dim[0] = x.dim[0] - y.dim[0];
      z.dim[1] = x.dim[1] - y.dim[1];
      z.dim[2] = x.dim[2] - y.dim[2];
      if (DF > 1) {
         for (size_t i = 0; i < 6; i++) {
            z.ddre[i] = x.ddre[i] - y.ddre[i];
            z.ddim[i] = x.ddim[i] - y.ddim[i];
         }
      }
   }
   return z;
}

gshe_complex cx_conj(size_t DF, gshe_complex z) {
   //finds the complex conjugate
   z.im *= -1.0;
   if (DF > 0) {
      z.dim[0] *= -1.0;
      z.dim[1] *= -1.0;
      z.dim[2] *= -1.0;

      if (DF > 1) {
         z.ddim[0] *= -1.0;
         z.ddim[1] *= -1.0;
         z.ddim[2] *= -1.0;
         z.ddim[3] *= -1.0;
         z.ddim[4] *= -1.0;
         z.ddim[5] *= -1.0;
      }
   }


   return z;
}

gshe_complex sc_add(size_t DF, gshe_complex x, double y) {
   //multiply gshe_complex x by real y
   x.re += y;
   x.im += y;
   if (DF > 0) {
      x.dre[0] += y;
      x.dre[1] += y;
      x.dre[2] += y;

      x.dim[0] += y;
      x.dim[1] += y;
      x.dim[2] += y;

      if (DF > 1) {
         x.ddre[0] += y;
         x.ddre[1] += y;
         x.ddre[2] += y;
         x.ddre[3] += y;
         x.ddre[4] += y;
         x.ddre[5] += y;

         x.ddim[0] += y;
         x.ddim[1] += y;
         x.ddim[2] += y;
         x.ddim[3] += y;
         x.ddim[4] += y;
         x.ddim[5] += y;
      }
   }
   return x;
}

gshe_complex sc_mult(size_t DF, gshe_complex x, double y) {
   //multiply gshe_complex x by real y
   x.re *= y;
   x.im *= y;
   if (DF > 0) {
      x.dre[0] *= y;
      x.dre[1] *= y;
      x.dre[2] *= y;

      x.dim[0] *= y;
      x.dim[1] *= y;
      x.dim[2] *= y;

      if (DF > 1) {
         x.ddre[0] *= y;
         x.ddre[1] *= y;
         x.ddre[2] *= y;
         x.ddre[3] *= y;
         x.ddre[4] *= y;
         x.ddre[5] *= y;

         x.ddim[0] *= y;
         x.ddim[1] *= y;
         x.ddim[2] *= y;
         x.ddim[3] *= y;
         x.ddim[4] *= y;
         x.ddim[5] *= y;
      }
   }
   return x;
}


size_t getSymCombo(size_t CS, size_t SS) {

   size_t CSb, SSb, iSymCombo = 7;
   size_t BASE_SYM_TRICLINIC = 1, BASE_SYM_CUBIC = 2, BASE_SYM_LOWER = 3;

   if (CS == SYM_CUBIC)        CSb = 2;
   else {
      if (CS == SYM_TRICLINIC)    CSb = 1;
      else { CSb = 3; }
   }

   if (SS == SYM_CUBIC)        SSb = 2;
   else {
      if (SS == SYM_TRICLINIC)    SSb = 1;
      else { SSb = 3; }
   }

   // This is the best way I could think to do this without rearranging the symmetry ids

   if (CSb == BASE_SYM_TRICLINIC && SSb == BASE_SYM_LOWER) 		iSymCombo = 1;
   if (CSb == BASE_SYM_LOWER && SSb == BASE_SYM_TRICLINIC)	iSymCombo = 2;
   if (CSb == BASE_SYM_LOWER && SSb == BASE_SYM_LOWER)		iSymCombo = 3;
   if (CSb == BASE_SYM_CUBIC && SSb == BASE_SYM_LOWER)		iSymCombo = 4;
   if (CSb == BASE_SYM_CUBIC && SSb == BASE_SYM_CUBIC)		iSymCombo = 5;
   if (CSb == BASE_SYM_CUBIC && SSb == BASE_SYM_TRICLINIC)	iSymCombo = 6;
   if (CSb == BASE_SYM_TRICLINIC && SSb == BASE_SYM_TRICLINIC)	iSymCombo = 7;

   return iSymCombo;
}

size_t sLinEq(size_t iSymCombo, size_t symm, size_t l) {
   size_t SizeLin;
   double x3, x4;
   // symmetry is matched with crystalSymmetry('symm').id
   SizeLin = 0;
   if (iSymCombo == 7) {
      // Special Case for Triclinic Triclinic
      SizeLin = 2 * l + 1;//
      return SizeLin;
   }


   switch (symm) {
   case SYM_ORTHORHOMBIC:
      if (l % 2) SizeLin = 0;//(l - 1) / 2;
      else SizeLin = l / 2 + 1;
      break;
   case SYM_TETRAGONAL:
      if (l % 2) SizeLin = (l + 4) / 4 - 1;
      else SizeLin = (l + 4) / 4;
      break;
   case SYM_CUBIC:
      if (l % 2) SizeLin = 0;//cubLin(l)
      else SizeLin = cubLin(l);
      break;
   case SYM_HEXAGONAL:
      if (l % 2) SizeLin = 0;//l / 6.0;
      else SizeLin = l / 6.0 + 1.0;
      break;

   case SYM_TRICLINIC:
      SizeLin = l + 1;//2*
      break;
   }

   return SizeLin;
}

void NumLinEq(size_t iSymCombo, size_t symm, size_t maxL, size_t* NumLin) {
   size_t l;
   // symmetry is matched with crystalSymmetry('symm').id

   if (iSymCombo == 7) {
      // Special Case for Triclinic Triclinic
      for (l = 0; l <= maxL; ++l) {
         NumLin[l] = 2 * l + 1;//
      }
      return;
   }


   switch (symm) {
   case SYM_ORTHORHOMBIC:
      for (l = 0; l <= maxL; ++l) {
         if (l % 2) NumLin[l] = 0;//(l - 1) / 2;
         else NumLin[l] = l / 2 + 1;
      }
      break;
   case SYM_TETRAGONAL:
      for (l = 0; l <= maxL; ++l) {
         if (l % 2) NumLin[l] = (l + 4) / 4 - 1;
         else NumLin[l] = (l + 4) / 4;
      }
      break;
   case SYM_CUBIC:
      for (l = 0; l <= maxL; ++l) {
         if (l % 2) NumLin[l] = 0;//cubLin(l)
         else NumLin[l] = cubLin(l);
      }
      break;
   case SYM_HEXAGONAL:
      for (l = 0; l <= maxL; ++l) {
         if (l % 2) NumLin[l] = 0;//l / 6.0;
         else NumLin[l] = l / 6.0 + 1.0;
      }
      break;

   case SYM_TRICLINIC:
      for (l = 0; l <= maxL; ++l) {
         NumLin[l] = l + 1;//2*
      }
      break;
   }
}

double normFactorN(size_t symm, size_t l, size_t nu, int* n) {
   // Normalization factors for Lower specimen symmetries
   double En;
   int Z = 0;
   switch (symm) {
   case SYM_ORTHORHOMBIC:
      Z = 2;  break;
   case SYM_HEXAGONAL:
      Z = 6;  break;
   case SYM_TETRAGONAL:
      Z = 4;  break;
   case SYM_TRIGONAL:
      Z = 3;  break;
   }

   if (l % 2) {
      n[0] = (Z * nu);
   }
   else {
      n[0] = (Z * (nu - 1));
   }

   if (n[0] == 0) En = 1.0;
   else En = sq2;

   return En;
}

double normFactorM(size_t symm, size_t mu, int* m) {
   // Normalization factors for Lower crystal symmetries
   double Em;
   int Z = 0;
   switch (symm) {
   case SYM_HEXAGONAL:
      Z = 6;  break;
   case SYM_TETRAGONAL:
      Z = 4;  break;
   case SYM_TRIGONAL:
      Z = 3;  break;
   case SYM_ORTHORHOMBIC:
      Z = 2;  break;
   case SYM_TRICLINIC:
      Z = 1;  break;
   }
   m[0] = (Z * (mu - 1));

   if (m[0] == 0.0) Em = 1.0;
   else Em = sq2;

   return Em;
}

int TriMap(size_t nu) {
   // maps mu or nu back to m or n
   int n;
   if (nu % 2) n = (int)((nu - 1) / 2);
   else n = -((int)(nu / 2));
   return n;
}


// ------------------------------------------------------------------------------
// The following functions are to calculate the T functions.
// See my notebook (Yale 1991) for reference pages 55-57.
// ------------------------------------------------------------------------------
gshe_complex symmetricT(size_t DF, size_t iSymCombo, size_t CS, size_t SS, size_t l, size_t mu, size_t nu, double sphi1[], double cphi1[], double sPhi[], double cPhi[], double sphi2[], double cphi2[]) {

   gshe_complex T;
   if (l == 0) {
      T = gshe_init(DF);
      T.re = 1.0;
      return T;
   }

   switch (iSymCombo) {
      /*case 1:
       * Tlmn = TriLowerT(l,mu,nu,phi1,PHI,phi2); break;*/
   case 2:
      T = lowerTricT(DF, CS, l, mu, nu, sphi1, cphi1, sPhi, cPhi, sphi2, cphi2); return T;
   case 3:
      T = lowLowerT(DF, CS, SS, l, mu, nu, sphi1, cphi1, sPhi, cPhi, sphi2, cphi2); return T;
   case 4:
      T = cubLowerT(DF, SS, l, mu, nu, sphi1, cphi1, sPhi, cPhi, sphi2, cphi2); return T;
   case 5:
      T = cubCubicT(DF, l, mu, nu, sphi1, cphi1, sPhi, cPhi, sphi2, cphi2); return T;
   case 6:
      T = cubTricT(DF, l, mu, nu, sphi1, cphi1, sPhi, cPhi, sphi2, cphi2); return T;
   case 7:
      T = tricTricT(DF, l, mu, nu, sphi1, cphi1, sPhi, cPhi, sphi2, cphi2); return T;
   }


}

gshe_complex lowerTricT(size_t DF, size_t CS, size_t l, size_t mu, size_t nu, double sphi1[], double cphi1[], double sPhi[], double cPhi[], double sphi2[], double cphi2[]) {
   gshe_complex Tlmn;
   double Em;
   int m, n, Z = 1;

   //Em = normFactorN(CS, l, mu, &m);

   // This is for the new map of triclinic
   n = (nu - 1);

   switch (CS) {
   case SYM_ORTHORHOMBIC:
      Z = 2;  break;
   case SYM_HEXAGONAL:
      Z = 6;  break;
   case SYM_TETRAGONAL:
      Z = 4;  break;
   case SYM_TRIGONAL:
      Z = 3;  break;
   }

   if (l % 2) {
      m = (Z * mu);
   }
   else {
      m = (Z * (mu - 1));
   }

   //if (m == 0) Em = 1.0;
   //else Em = sq2;
   //if (n != 0.0) Em *= sq2;

   if (m != 0.0 && n != 0.0) Em = 2.0;
   else if (m != 0.0 || n != 0.0) Em = sq2;
   else Em = 1.0;


   Tlmn = Rlmn(DF, l, m, n, sphi1, cphi1, sPhi, cPhi, sphi2, cphi2);
   Tlmn = sc_mult(DF, Tlmn, Em);

   return Tlmn;
}

gshe_complex lowLowerT(size_t DF, size_t CS, size_t SS, size_t l, size_t mu, size_t nu, double sphi1[], double cphi1[], double sPhi[], double cPhi[], double sphi2[], double cphi2[]) {

   gshe_complex Tlmn;
   double		Em, En;
   int m, n;


   Em = normFactorM(CS, mu, &m);
   En = normFactorN(SS, l, nu, &n);

   Tlmn = Slmn(DF, l, m, n, sphi1, cphi1, sPhi, cPhi, sphi2, cphi2);

   Tlmn = sc_mult(DF, Tlmn, Em * En);
   return Tlmn;
}

gshe_complex cubLowerT(size_t DF, size_t SS, size_t l, size_t mu, size_t nu, double sphi1[], double cphi1[], double sPhi[], double cPhi[], double sphi2[], double cphi2[]) {

   gshe_complex Tlmn, S;
   double En;
   int m, n;
   En = normFactorN(SS, l, nu, &n);
   Tlmn = gshe_init(DF);

   En *= sq2pi;

   for (m = 0; m <= l; m += 4) {
      //need a double for Slmn

      S = Slmn(DF, l, m, n, sphi1, cphi1, sPhi, cPhi, sphi2, cphi2);
      S = sc_mult(DF, S, blmn(l, m / 4, mu));

      Tlmn = cx_add(DF, Tlmn, S);
   }

   Tlmn = sc_mult(DF, Tlmn, En);
   return Tlmn;

}

gshe_complex cubCubicT(size_t DF, size_t l, size_t mu, size_t nu, double sphi1[], double cphi1[], double sPhi[], double cPhi[], double sphi2[], double cphi2[]) {

   gshe_complex Tlmn, S;
   size_t	m, n;
   Tlmn = gshe_init(DF);
   for (m = 0; m <= l; m += 4) for (n = 0; n <= l; n += 4) {

      S = Slmn(DF, l, m, n, sphi1, cphi1, sPhi, cPhi, sphi2, cphi2);
      S = sc_mult(DF, S, blmn(l, m / 4, mu) * blmn(l, n / 4, nu));
      Tlmn = cx_add(DF, Tlmn, S);
   }
   Tlmn = sc_mult(DF, Tlmn, TWOPI);
   return Tlmn;
}

gshe_complex cubTricT(size_t DF, size_t l, size_t mu, size_t nu, double sphi1[], double cphi1[], double sPhi[], double cPhi[], double sphi2[], double cphi2[]) {
   double En;
   int m, n;
   gshe_complex Tlmunu, T2;
   En = normFactorM(SYM_TRICLINIC, nu, &n);

   Tlmunu = noSymT(DF, l, 0, n, sphi1, cphi1, sPhi, cPhi, sphi2, cphi2);

   Tlmunu = sc_mult(DF, Tlmunu, blmn(l, (size_t)0, mu));

   for (m = 4; m <= l; m += 4)
   {
      T2 = Rlmn(DF, l, m, n, sphi1, cphi1, sPhi, cPhi, sphi2, cphi2);

      T2 = sc_mult(DF, T2, blmn(l, (size_t)(m) / 4, mu));

      Tlmunu = cx_add(DF, Tlmunu, T2);
   }

   Tlmunu = sc_mult(DF, Tlmunu, En * sq2pi);

   return Tlmunu;
}

gshe_complex tricTricT(size_t DF, size_t l, size_t mu, size_t nu, double sphi1[], double cphi1[], double sPhi[], double cPhi[], double sphi2[], double cphi2[]) {
   gshe_complex Tlmn;
   int m, n;
   m = TriMap(mu);
   n = TriMap(nu);

   Tlmn = noSymT(DF, l, m, n, sphi1, cphi1, sPhi, cPhi, sphi2, cphi2);
   return Tlmn;
}

gshe_complex noSymT(size_t DF, size_t l, int m, int n, double sphi1[], double cphi1[], double sPhi[], double cPhi[], double sphi2[], double cphi2[]) {
   gshe_complex Tlmn;
   double	Plmn = 0.0;
   size_t	 s, mp, np;
   double cc, ss, cs, sc;
   double cc1, ss1, cs1, sc1, cc2, ss2, cs2, sc2;
   double cc11, ss11, cs11, sc11, cc22, ss22, cs22, sc22, cc12, ss12, cs12, sc12;
   double almnsv[RANK_MAX + 1];


   mp = abs(m); np = abs(n);

   Tlmn = gshe_init(DF);

   for (s = 0; s <= l; s++) almnsv[s] = almnsp(l, mp, np, s);

   cc = cphi1[np] * cphi2[mp];//cos(n*phi1)*cos(m*phi2);
   ss = sphi1[np] * sphi2[mp] * fsign(m * n);//sin(n*phi1)*sin(m*phi2);
   cs = cphi1[np] * sphi2[mp] * fsign(m);//cos(n*phi1)*sin(m*phi2);
   sc = sphi1[np] * cphi2[mp] * fsign(n);//sin(n*phi1)*cos(m*phi2);

   if (DF > 0)
   {
      cc1 = -n * sc;
      ss1 = n * cs;
      cs1 = -n * ss;
      sc1 = n * cc;

      cc2 = -m * cs;
      ss2 = m * sc;
      cs2 = m * cc;
      sc2 = -m * ss;
      if (DF > 1)
      {
         double c_mult = -n * n;

         cc11 = c_mult * cc;
         ss11 = c_mult * ss;
         cs11 = c_mult * cs;
         sc11 = c_mult * sc;

         c_mult = -m * m;

         cc22 = c_mult * cc;
         ss22 = c_mult * ss;
         cs22 = c_mult * cs;
         sc22 = c_mult * sc;

         c_mult = m * n;

         cc12 = c_mult * ss;
         ss12 = c_mult * cc;
         cs12 = -c_mult * sc;
         sc12 = -c_mult * cs;
      }
   }



   //by hari
   if (!((mp + np) % 2))
   {  // --- m+n even --- Plmn is purely real
      if (m * n < 0.0) { for (s = 0; s <= l; ++s) Plmn += sign(l + s) * almnsv[s] * cPhi[s]; }

      else { for (s = 0; s <= l; ++s) Plmn += almnsv[s] * cPhi[s]; }

      Tlmn.re = (cc - ss) * Plmn;  Tlmn.im = (cs + sc) * Plmn;

      if (DF > 0) {
         double dPlmn = 0.0;
         if (m * n < 0) { for (s = 1; s <= l; ++s) dPlmn += -sign(l + s) * almnsv[s] * sPhi[s] * s; }

         else { for (s = 1; s <= l; ++s) dPlmn += -almnsv[s] * sPhi[s] * s; }

         Tlmn.dre[0] = (cc1 - ss1) * Plmn;
         Tlmn.dre[1] = (cc - ss) * dPlmn;
         Tlmn.dre[2] = (cc2 - ss2) * Plmn;

         Tlmn.dim[0] = (cs1 + sc1) * Plmn;
         Tlmn.dim[1] = (cs + sc) * dPlmn;
         Tlmn.dim[2] = (cs2 + sc2) * Plmn;

         if (DF > 1) {
            double ddPlmn = 0.0;
            if (m * n < 0.0) {
               for (s = 0; s <= l; ++s) {
                  ddPlmn -= sign(l + s) * almnsv[s] * cPhi[s] * (double)(s * s);
               }
            }

            else {
               for (s = 0; s <= l; ++s) {
                  ddPlmn -= almnsv[s] * cPhi[s] * (double)(s * s);
               }
            }

            Tlmn.ddre[0] = (cc11 - ss11) * Plmn;
            Tlmn.ddre[1] = (cc1 - ss1) * dPlmn;
            Tlmn.ddre[2] = (cc12 - ss12) * Plmn;
            Tlmn.ddre[3] = (cc - ss) * ddPlmn;
            Tlmn.ddre[4] = (cc2 - ss2) * dPlmn;
            Tlmn.ddre[5] = (cc22 - ss22) * Plmn;

            Tlmn.ddim[0] = (cs11 + sc11) * Plmn;
            Tlmn.ddim[1] = (cs1 + sc1) * dPlmn;
            Tlmn.ddim[2] = (cs12 + sc12) * Plmn;
            Tlmn.ddim[3] = (cs + sc) * ddPlmn;
            Tlmn.ddim[4] = (cs2 + sc2) * dPlmn;
            Tlmn.ddim[5] = (cs22 + sc22) * Plmn;
         }
      }
   }
   else { //--- m+n odd --- Plmn is purely imaginary
      if (m * n < 0) for (s = 1; s <= l; ++s) { Plmn += sign(l + s) * almnsv[s] * sPhi[s]; }

      else { for (s = 1; s <= l; ++s) { Plmn += almnsv[s] * sPhi[s]; } }

      Tlmn.re = -(cs + sc) * Plmn;  Tlmn.im = (cc - ss) * Plmn;

      if (DF > 0) {
         double dPlmn = 0.0;

         if (m * n < 0) { for (s = 1; s <= l; ++s) dPlmn += sign(l + s) * almnsv[s] * cPhi[s] * s; }

         else { for (s = 1; s <= l; ++s) dPlmn += almnsv[s] * cPhi[s] * s; }

         Tlmn.dre[0] = -(cs1 + sc1) * Plmn;
         Tlmn.dre[1] = -(cs + sc) * dPlmn;
         Tlmn.dre[2] = -(cs2 + sc2) * Plmn;

         Tlmn.dim[0] = (cc1 - ss1) * Plmn;
         Tlmn.dim[1] = (cc - ss) * dPlmn;
         Tlmn.dim[2] = (cc2 - ss2) * Plmn;

         if (DF > 1) {
            double ddPlmn = 0.0;

            if (m * n < 0.0) {
               for (s = 0; s <= l; ++s) {
                  ddPlmn -= sign(l + s) * almnsv[s] * sPhi[s] * (s * s);
               }
            }

            else {
               for (s = 0; s <= l; ++s) {
                  ddPlmn -= almnsv[s] * sPhi[s] * (s * s);
               }
            }

            Tlmn.ddre[0] = -(cs11 + sc11) * Plmn;
            Tlmn.ddre[1] = -(cs1 + sc1) * dPlmn;
            Tlmn.ddre[2] = -(cs12 + sc12) * Plmn;
            Tlmn.ddre[3] = -(cs + sc) * ddPlmn;
            Tlmn.ddre[4] = -(cs2 + sc2) * dPlmn;
            Tlmn.ddre[5] = -(cs22 + sc22) * Plmn;

            Tlmn.ddim[0] = (cc11 - ss11) * Plmn;
            Tlmn.ddim[1] = (cc1 - ss1) * dPlmn;
            Tlmn.ddim[2] = (cc12 - ss12) * Plmn;
            Tlmn.ddim[3] = (cc - ss) * ddPlmn;
            Tlmn.ddim[4] = (cc2 - ss2) * dPlmn;
            Tlmn.ddim[5] = (cc22 - ss22) * Plmn;
         }
      }
   }
   return Tlmn;
}

gshe_complex Slmn(size_t DF, size_t l, int m, int n, double sphi1[], double cphi1[], double sPhi[], double cPhi[], double sphi2[], double cphi2[]) {

   double	sum1 = 0.0, sum2 = 0.0;
   double cc, ss, cs, sc;
   double cc1, ss1, cc2, ss2;
   double cc11, ss11, cc22, ss22, cc12, ss12;
   size_t np, mp;
   gshe_complex	rS;
   size_t 	s;
   double almnsv[RANK_MAX + 1];


   mp = abs(m); np = abs(n);

   for (s = 0; s <= l; s++) almnsv[s] = almnsp(l, mp, np, s);

   rS = gshe_init(DF);

   //by hari
   cc = cphi1[np] * cphi2[mp];//cos(n*phi1)*cos(m*phi2);
   ss = sphi1[np] * sphi2[mp] * fsign(m * n);//sin(n*phi1)*sin(m*phi2);


   if (DF > 0)
   {
      cs = cphi1[np] * sphi2[mp] * fsign(m);//cos(n*phi1)*sin(m*phi2);
      sc = sphi1[np] * cphi2[mp] * fsign(n);//sin(n*phi1)*cos(m*phi2);

      cc1 = -n * sc;
      ss1 = n * cs;

      cc2 = -m * cs;
      ss2 = m * sc;
      if (DF > 1)
      {
         double c_mult = -n * n;

         cc11 = c_mult * cc;
         ss11 = c_mult * ss;

         c_mult = -m * m;

         cc22 = c_mult * cc;
         ss22 = c_mult * ss;

         c_mult = m * n;

         cc12 = c_mult * ss;
         ss12 = c_mult * cc;
      }
   }

   for (s = 0; s <= l; s += 2) { sum1 += almnsv[s] * cPhi[s]; }
   for (s = 1; s <= l; s += 2) { sum2 += almnsv[s] * cPhi[s]; }

   rS.re = sum1 * cc - sum2 * ss;

   if (DF > 0) {
      double dsum1 = 0.0, dsum2 = 0.0;

      for (s = 0; s <= l; s += 2) { dsum1 += -almnsv[s] * sPhi[s] * s; }
      for (s = 1; s <= l; s += 2) { dsum2 += -almnsv[s] * sPhi[s] * s; }

      rS.dre[0] = sum1 * cc1 - sum2 * ss1;
      rS.dre[1] = dsum1 * cc - dsum2 * ss;
      rS.dre[2] = sum1 * cc2 - sum2 * ss2;

      if (DF > 1) {
         double ddsum1 = 0.0, ddsum2 = 0.0;
         for (s = 0; s <= l; s += 2) { ddsum1 += -almnsv[s] * cPhi[s] * s * s; }
         for (s = 1; s <= l; s += 2) { ddsum2 += -almnsv[s] * cPhi[s] * s * s; }
         rS.ddre[0] = sum1 * cc11 - sum2 * ss11;
         rS.ddre[1] = dsum1 * cc1 - dsum2 * ss1;
         rS.ddre[2] = sum1 * cc12 - sum2 * ss12;
         rS.ddre[3] = ddsum1 * cc - ddsum2 * ss;
         rS.ddre[4] = dsum1 * cc2 - dsum2 * ss2;
         rS.ddre[5] = sum1 * cc22 - sum2 * ss22;
      }


   }
   return rS;
}

gshe_complex Rlmn(size_t DF, size_t l, int m, int n, double sphi1[], double cphi1[], double sPhi[], double cPhi[], double sphi2[], double cphi2[]) {

   gshe_complex R;
   double  sum1, sum2, lf, cc, ss, cs, sc;
   size_t 	s, mp, np;
   double cc1, ss1, cs1, sc1, cc2, ss2, cs2, sc2;
   double cc11, ss11, cs11, sc11, cc22, ss22, cs22, sc22, cc12, ss12, cs12, sc12;
   double almnsv[RANK_MAX + 1];



   R = gshe_init(DF);

   mp = abs(m); np = abs(n);

   for (s = 0; s <= l; s++) almnsv[s] = almnsp(l, mp, np, s);

   cc = cphi1[np] * cphi2[mp];//cos(n*phi1)*cos(m*phi2);
   ss = sphi1[np] * sphi2[mp] * fsign(m * n);//sin(n*phi1)*sin(m*phi2);
   cs = cphi1[np] * sphi2[mp] * fsign(m);//cos(n*phi1)*sin(m*phi2);
   sc = sphi1[np] * cphi2[mp] * fsign(n);//sin(n*phi1)*cos(m*phi2);

   if (DF > 0)
   {
      cc1 = -n * sc;
      ss1 = n * cs;
      cs1 = -n * ss;
      sc1 = n * cc;

      cc2 = -m * cs;
      ss2 = m * sc;
      cs2 = m * cc;
      sc2 = -m * ss;
      if (DF > 1)
      {
         double c_mult = -n * n;

         cc11 = c_mult * cc;
         ss11 = c_mult * ss;
         cs11 = c_mult * cs;
         sc11 = c_mult * sc;

         c_mult = -m * m;

         cc22 = c_mult * cc;
         ss22 = c_mult * ss;
         cs22 = c_mult * cs;
         sc22 = c_mult * sc;

         c_mult = m * n;

         cc12 = c_mult * ss;
         ss12 = c_mult * cc;
         cs12 = -c_mult * sc;
         sc12 = -c_mult * cs;
      }
   }

   if (np % 2) {  /* --- n odd --- */
      sum1 = 0.0;
      sum2 = 0.0;
      /* for (s=0; s<=l; s+=2) sum1 += a[l][m][n][s]*sin(s*PHI); */
      for (s = 1; s <= l; s += 2) sum2 += almnsv[s] * sPhi[s];
      for (s = 2; s <= l; s += 2) sum1 += almnsv[s] * sPhi[s];


      if (n < 0.0) {
         lf = sign(l);
         sum1 *= lf; sum2 *= (-lf);
      }


      R.re = -sc * sum1 - cs * sum2;  R.im = cc * sum1 - ss * sum2;

      if (DF > 0) {
         double dsum1 = 0.0, dsum2 = 0.0;
         for (s = 2; s <= l; s += 2) dsum1 += almnsv[s] * cPhi[s] * ((double)s);
         for (s = 1; s <= l; s += 2) dsum2 += almnsv[s] * cPhi[s] * ((double)s);

         if (n < 0.0) {
            dsum1 *= lf; dsum2 *= (-lf);
         }

         R.dre[0] = -sc1 * sum1 - cs1 * sum2;
         R.dre[1] = -sc * dsum1 - cs * dsum2;
         R.dre[2] = -sc2 * sum1 - cs2 * sum2;
         R.dim[0] = cc1 * sum1 - ss1 * sum2;
         R.dim[1] = cc * dsum1 - ss * dsum2;
         R.dim[2] = cc2 * sum1 - ss2 * sum2;
         if (DF > 1) {
            double ddsum1 = 0.0, ddsum2 = 0.0;
            for (s = 2; s <= l; s += 2) ddsum1 += -almnsv[s] * sPhi[s] * ((double)(s * s));
            for (s = 1; s <= l; s += 2) ddsum2 += -almnsv[s] * sPhi[s] * ((double)(s * s));

            if (n < 0.0) {
               ddsum1 *= lf; ddsum2 *= (-lf);
            }

            R.ddre[0] = -sc11 * sum1 - cs11 * sum2;
            R.ddre[1] = -sc1 * dsum1 - cs1 * dsum2;
            R.ddre[2] = -sc12 * sum1 - cs12 * sum2;
            R.ddre[3] = -sc * ddsum1 - cs * ddsum2;
            R.ddre[4] = -sc2 * dsum1 - cs2 * dsum2;
            R.ddre[5] = -sc22 * sum1 - cs22 * sum2;
            R.ddim[0] = cc11 * sum1 - ss11 * sum2;
            R.ddim[1] = cc1 * dsum1 - ss1 * dsum2;
            R.ddim[2] = cc12 * sum1 - ss12 * sum2;
            R.ddim[3] = cc * ddsum1 - ss * ddsum2;
            R.ddim[4] = cc2 * dsum1 - ss2 * dsum2;
            R.ddim[5] = cc22 * sum1 - ss22 * sum2;
         }

      }

   }
   else { /* --- n even --- */

      sum1 = 0.0;
      sum2 = 0.0;

      for (s = 0; s <= l; s += 2) sum1 += almnsv[s] * cPhi[s];
      for (s = 1; s <= l; s += 2) sum2 += almnsv[s] * cPhi[s];

      if (n < 0.0) {
         lf = sign(l);
         sum1 *= lf; sum2 *= (-lf);
      }
      R.re = cc * sum1 - ss * sum2;  R.im = sc * sum1 + cs * sum2;

      if (DF > 0) {
         double dsum1 = 0.0, dsum2 = 0.0;
         for (s = 2; s <= l; s += 2) dsum1 += -almnsv[s] * sPhi[s] * ((double)s);
         for (s = 1; s <= l; s += 2) dsum2 += -almnsv[s] * sPhi[s] * ((double)s);

         if (n < 0.0) {
            dsum1 *= lf; dsum2 *= (-lf);
         }

         R.dre[0] = cc1 * sum1 - ss1 * sum2;
         R.dre[1] = cc * dsum1 - ss * dsum2;
         R.dre[2] = cc2 * sum1 - ss2 * sum2;

         R.dim[0] = sc1 * sum1 + cs1 * sum2;
         R.dim[1] = sc * dsum1 + cs * dsum2;
         R.dim[2] = sc2 * sum1 + cs2 * sum2;

         if (DF > 1) {
            double ddsum1 = 0.0, ddsum2 = 0.0;
            for (s = 2; s <= l; s += 2) ddsum1 += -almnsv[s] * cPhi[s] * ((double)(s * s));
            for (s = 1; s <= l; s += 2) ddsum2 += -almnsv[s] * cPhi[s] * ((double)(s * s));

            if (n < 0.0) {
               ddsum1 *= lf; ddsum2 *= (-lf);
            }

            R.ddre[0] = cc11 * sum1 - ss11 * sum2;
            R.ddre[1] = cc1 * dsum1 - ss1 * dsum2;
            R.ddre[2] = cc12 * sum1 - ss12 * sum2;
            R.ddre[3] = cc * ddsum1 - ss * ddsum2;
            R.ddre[4] = cc2 * dsum1 - ss2 * dsum2;
            R.ddre[5] = cc22 * sum1 - ss22 * sum2;

            R.ddim[0] = sc11 * sum1 + cs11 * sum2;
            R.ddim[1] = sc1 * dsum1 + cs1 * dsum2;
            R.ddim[2] = sc12 * sum1 + cs12 * sum2;
            R.ddim[3] = sc * ddsum1 + cs * ddsum2;
            R.ddim[4] = sc2 * dsum1 + cs2 * dsum2;
            R.ddim[5] = sc22 * sum1 + cs22 * sum2;
         }

      }
   }

   return R;
}

size_t calcTLEN(size_t maxL, size_t CS, size_t SS)
{
   size_t L, incr, iSymCombo;

   iSymCombo = getSymCombo(CS, SS);

   if (iSymCombo == 4 || iSymCombo == 3) incr = 1;
   else incr = 2;


   L = 0;

   for (size_t l = 0; l <= maxL; l++) {
      L += incr * sLinEq(iSymCombo, CS, l) * sLinEq(iSymCombo, SS, l);
   }

   return L;
}

size_t calcTLENK(size_t maxL, size_t CS, size_t SS, double* kernel)
{
   size_t L, incr, iSymCombo;

   iSymCombo = getSymCombo(CS, SS);

   if (iSymCombo == 4 || iSymCombo == 3) incr = 1;
   else incr = 2;


   L = 0;

   for (size_t l = 0; l <= maxL; l++) {
      if (kernel[l] != 0.0) {
         L += incr * sLinEq(iSymCombo, CS, l) * sLinEq(iSymCombo, SS, l);
      }
   }

   return L;
}

double calcTDISCL(size_t maxL, size_t CS) {
   double  TDIS;

   TDIS = 1.00;
   switch (CS) {
   case SYM_TRICLINIC:
      TDIS = sqrt(2) / (maxL + 1);
      break;
   case SYM_ORTHORHOMBIC:
      TDIS = sqrt(2) / (maxL + 1);
      break;
   case SYM_HEXAGONAL:
      TDIS = sqrt(6) / (maxL + 1);
      break;
   case SYM_CUBIC:
      TDIS = sqrt(12) / (maxL + 1);
      break;
   }

   return TDIS;
}

void calcF(double* mF, double* phi, double* kernel, size_t CS, size_t SS, size_t maxL, size_t N) {
   size_t DF, iSymCombo, incr;
   size_t L;
   double sum = 0.0;
   //Symmetry Assignments
   size_t* mLinEq, * nLinEq;

   if (maxL > RANK_MAX) maxL = RANK_MAX;

   mLinEq = (size_t*)calloc(maxL + 1, sizeof(size_t));
   nLinEq = (size_t*)calloc(maxL + 1, sizeof(size_t));

   iSymCombo = getSymCombo(CS, SS);

   NumLinEq(iSymCombo, CS, maxL, mLinEq);
   NumLinEq(iSymCombo, SS, maxL, nLinEq);

   if (iSymCombo == 4 || iSymCombo == 3) {
      incr = 1;
   }
   else {
      incr = 2;
   }

   L = 0;

   for (size_t l = 0; l <= maxL; l++) {
      for (size_t mu = 1; mu <= mLinEq[l]; mu++) {
         for (size_t nu = 1; nu <= nLinEq[l]; nu++) {
            L += incr;
         }
      }
   }

   DF = 0;

#pragma omp parallel default(none) firstprivate(DF, iSymCombo, CS, SS, maxL, L, N, incr, mLinEq, nLinEq) shared(mF, phi, kernel)
   {
      size_t  j, l, mu, nu;

      double* sphi1 = (double*)calloc(maxL + 1, sizeof(double));
      double* cphi1 = (double*)calloc(maxL + 1, sizeof(double));
      double* sPhi = (double*)calloc(maxL + 1, sizeof(double));
      double* cPhi = (double*)calloc(maxL + 1, sizeof(double));
      double* sphi2 = (double*)calloc(maxL + 1, sizeof(double));
      double* cphi2 = (double*)calloc(maxL + 1, sizeof(double));

#pragma omp for// schedule(static)
      for (size_t i = 0; i < N; i++) {
         //pre-calculate all the sines and cosines we will need
         for (j = 0; j <= maxL; j++) {
            sphi1[j] = sin((double)j * phi[3 * i]);
            cphi1[j] = cos((double)j * phi[3 * i]);
            sPhi[j] = sin((double)j * phi[1 + 3 * i]);
            cPhi[j] = cos((double)j * phi[1 + 3 * i]);
            sphi2[j] = sin((double)j * phi[2 + 3 * i]);
            cphi2[j] = cos((double)j * phi[2 + 3 * i]);
         }

         // cnt counts the index of the fourier coeff
         size_t cnt;
         cnt = 0;
         for (l = 0; l <= maxL; l++) {
            double pd_fac = kernel[l];
            // mu represents m for the crystal symmetry
            for (mu = 1; mu <= mLinEq[l]; mu++) {
               //nu represents n for the sample symmetry
               for (nu = 1; nu <= nLinEq[l]; nu++) {
                  //size_t l, mu, nu;
                  gshe_complex Tlmn;
                  // Tlmn is a structure that contains both the real and imaginary parts of the coeffs, as well as the derivatives
                  Tlmn = symmetricT(DF, iSymCombo, CS, SS, l, mu, nu, sphi1, cphi1, sPhi, cPhi, sphi2, cphi2);
                  Tlmn = sc_mult(DF, Tlmn, pd_fac);

                  mF[cnt + i * L] = Tlmn.re;
                  cnt++;

                  if (incr > 1) {
                     mF[cnt + i * L] = Tlmn.im;
                     cnt++;
                  }
               }
            }
         }
      }

      free(sphi1); free(cphi1);
      free(sPhi); free(cPhi);
      free(sphi2); free(cphi2);
   }

   //mexPrintf("exit mf\n");

   free(mLinEq);
   free(nLinEq);
}

void calcMF(double* mF, double* phi, double* alpha, double* kernel, size_t CS, size_t SS, size_t maxL, size_t N) {
   size_t DF, iSymCombo, incr;
   size_t L;
   double sum = 0.0, iA;
   //Symmetry Assignments
   size_t* mLinEq, * nLinEq;

   if (maxL > RANK_MAX) maxL = RANK_MAX;

   mLinEq = (size_t*)calloc(maxL + 1, sizeof(size_t));
   nLinEq = (size_t*)calloc(maxL + 1, sizeof(size_t));

   iSymCombo = getSymCombo(CS, SS);

   NumLinEq(iSymCombo, CS, maxL, mLinEq);
   NumLinEq(iSymCombo, SS, maxL, nLinEq);

   if (iSymCombo == 4 || iSymCombo == 3) {
      incr = 1;
   }
   else {
      incr = 2;
   }

   L = 0;

   for (size_t l = 0; l <= maxL; l++) {
      for (size_t mu = 1; mu <= mLinEq[l]; mu++) {
         for (size_t nu = 1; nu <= nLinEq[l]; nu++) {
            L += incr;
         }
      }
   }

   {
      double A;
      kahansum(N, alpha[k], A)
         iA = 1.0 / A;
      //for (size_t i = 0; i < N; i++) alpha[i] *= iA;
   }


   DF = 0;

   {
      for (size_t j = 0; j < L; j++) {
         mF[j] = 0.0;
      }

      double* pmF;
      
      

#pragma omp parallel default(none) firstprivate(DF, iSymCombo, CS, SS, maxL, L, N, incr, mLinEq, nLinEq, iA) shared(mF, kernel, phi, pmF, alpha)
      {
         double* sphi1;
         double* cphi1;
         double* sPhi;
         double* cPhi;
         double* sphi2;
         double* cphi2;
         double* mFc;
         const size_t nthread = omp_get_num_threads();
         const size_t ithread = omp_get_thread_num();
         size_t  j, l, mu, nu;
         int i;

         sphi1 = (double*)calloc(maxL + 1, sizeof(double));
         cphi1 = (double*)calloc(maxL + 1, sizeof(double));
         sPhi = (double*)calloc(maxL + 1, sizeof(double));
         cPhi = (double*)calloc(maxL + 1, sizeof(double));
         sphi2 = (double*)calloc(maxL + 1, sizeof(double));
         cphi2 = (double*)calloc(maxL + 1, sizeof(double));
         mFc = (double*)malloc(L * sizeof(double));

         for (j = 0; j < L; j++) mFc[j] = 0.0;
         
         //printf("%d\n", nthread);

#pragma omp single
         {
            pmF = (double*)malloc(L * nthread * sizeof(double));
         }

         //double*rmF = &pmF[L * ithread];

         for (i = 0; i < L; i++) {
            pmF[i + L * ithread] = 0.0;
         }

#pragma omp for// schedule(static)
         for (i = 0; i < N; i++) {

            double calpha;
            calpha = alpha[i] * iA;
            //pre-calculate all the sines and cosines we will need
            for (j = 0; j <= maxL; j++) {
               sphi1[j] = sin((double)j * phi[3 * i]);
               cphi1[j] = cos((double)j * phi[3 * i]);
               sPhi[j] = sin((double)j * phi[1 + 3 * i]);
               cPhi[j] = cos((double)j * phi[1 + 3 * i]);
               sphi2[j] = sin((double)j * phi[2 + 3 * i]);
               cphi2[j] = cos((double)j * phi[2 + 3 * i]);
            }

            // cnt counts the index of the fourier coeff
            size_t cnt;
            cnt = 0;
            for (l = 0; l <= maxL; l++) {
               double pd_fac = kernel[l];
               // mu represents m for the crystal symmetry
               for (mu = 1; mu <= mLinEq[l]; mu++) {
                  //nu represents n for the sample symmetry
                  for (nu = 1; nu <= nLinEq[l]; nu++) {
                     //size_t l, mu, nu;
                     double kahany, kahant;
                     gshe_complex Tlmn;
                     // Tlmn is a structure that contains both the real and imaginary parts of the coeffs, as well as the derivatives
                     Tlmn = symmetricT(DF, iSymCombo, CS, SS, l, mu, nu, sphi1, cphi1, sPhi, cPhi, sphi2, cphi2);
                     Tlmn = sc_mult(DF, Tlmn, pd_fac);

                     kahany = (calpha * Tlmn.re) - mFc[cnt];
                     kahant = pmF[cnt + L * ithread] + kahany;

                     mFc[cnt] = (kahant - pmF[cnt + L * ithread]) - kahany;
                     pmF[cnt + L * ithread] = kahant;
                     cnt++;

                     if (incr > 1) {
                        kahany = (calpha * Tlmn.im) - mFc[cnt];
                        kahant = pmF[cnt + L * ithread] + kahany;

                        mFc[cnt] = (kahant - pmF[cnt + L * ithread]) - kahany;
                        pmF[cnt + L * ithread] = kahant;
                        cnt++;
                     }
                  }
               }
            }
         }
#pragma omp for
         for (i = 0; i < L; i++) {
            double ssum;
            ssum = 0.0;
            kahansum(nthread, pmF[i + k * L], ssum)
               mF[i] = ssum;
         }

         free(sphi1); free(cphi1);
         free(sPhi); free(cPhi);
         free(sphi2); free(cphi2);
         free(mFc);
      }
      free(pmF);
   }


   mF[0] = 1.0;

   free(mLinEq);
   free(nLinEq);
}

void calcMFK(double* mF, double* phi, double* alpha, double* kernel, size_t CS, size_t SS, size_t maxL, size_t N) {
   size_t DF, iSymCombo, incr;
   size_t L;
   double sum = 0.0, iA;
   //Symmetry Assignments
   size_t* mLinEq, * nLinEq;

   if (maxL > RANK_MAX) maxL = RANK_MAX;

   mLinEq = (size_t*)calloc(maxL + 1, sizeof(size_t));
   nLinEq = (size_t*)calloc(maxL + 1, sizeof(size_t));

   iSymCombo = getSymCombo(CS, SS);

   NumLinEq(iSymCombo, CS, maxL, mLinEq);
   NumLinEq(iSymCombo, SS, maxL, nLinEq);

   if (iSymCombo == 4 || iSymCombo == 3) {
      incr = 1;
   }
   else {
      incr = 2;
   }

   L = 0;

   for (size_t l = 0; l <= maxL; l++) {
      if (kernel[l] != 0.0) {
         for (size_t mu = 1; mu <= mLinEq[l]; mu++) {
            for (size_t nu = 1; nu <= nLinEq[l]; nu++) {
               L += incr;
            }
         }
      }
   }

   {
      double A;
      kahansum(N, alpha[k], A)
         iA = 1.0 / A;
      //for (size_t i = 0; i < N; i++) alpha[i] *= iA;
   }


   DF = 0;

   {
      for (size_t j = 0; j < L; j++) {
         mF[j] = 0.0;
      }

      double* pmF;

#pragma omp parallel default(none) firstprivate(DF, iSymCombo, CS, SS, maxL, L, N, incr, mLinEq, nLinEq, iA) shared(mF, kernel, phi, pmF, alpha)
      {
         double* sphi1;
         double* cphi1;
         double* sPhi;
         double* cPhi;
         double* sphi2;
         double* cphi2;
         double* mFc;
         const size_t nthread = omp_get_num_threads();;
         const size_t ithread = omp_get_thread_num();
         size_t  j, l, mu, nu;
         int i;

         sphi1 = (double*)calloc(maxL + 1, sizeof(double));
         cphi1 = (double*)calloc(maxL + 1, sizeof(double));
         sPhi = (double*)calloc(maxL + 1, sizeof(double));
         cPhi = (double*)calloc(maxL + 1, sizeof(double));
         sphi2 = (double*)calloc(maxL + 1, sizeof(double));
         cphi2 = (double*)calloc(maxL + 1, sizeof(double));
         mFc = (double*)malloc(L * sizeof(double));

         for (j = 0; j < L; j++) mFc[j] = 0.0;

#pragma omp single
         {
            pmF = (double*)malloc(L * nthread * sizeof(double));
         }

         //double*rmF = &pmF[L * ithread];

         for (i = 0; i < L; i++) {
            pmF[i + L * ithread] = 0.0;
         }

#pragma omp for// schedule(static)
         for (i = 0; i < N; i++) {

            double calpha;
            calpha = alpha[i] * iA;
            //pre-calculate all the sines and cosines we will need
            for (j = 0; j <= maxL; j++) {
               sphi1[j] = sin((double)j * phi[3 * i]);
               cphi1[j] = cos((double)j * phi[3 * i]);
               sPhi[j] = sin((double)j * phi[1 + 3 * i]);
               cPhi[j] = cos((double)j * phi[1 + 3 * i]);
               sphi2[j] = sin((double)j * phi[2 + 3 * i]);
               cphi2[j] = cos((double)j * phi[2 + 3 * i]);
            }

            // cnt counts the index of the fourier coeff
            size_t cnt;
            cnt = 0;
            for (l = 0; l <= maxL; l++) {
               double pd_fac = kernel[l];
               if (pd_fac != 0.0) {
                  // mu represents m for the crystal symmetry
                  for (mu = 1; mu <= mLinEq[l]; mu++) {
                     //nu represents n for the sample symmetry
                     for (nu = 1; nu <= nLinEq[l]; nu++) {
                        //size_t l, mu, nu;
                        double kahany, kahant;
                        gshe_complex Tlmn;
                        // Tlmn is a structure that contains both the real and imaginary parts of the coeffs, as well as the derivatives
                        Tlmn = symmetricT(DF, iSymCombo, CS, SS, l, mu, nu, sphi1, cphi1, sPhi, cPhi, sphi2, cphi2);
                        Tlmn = sc_mult(DF, Tlmn, pd_fac);

                        kahany = (calpha * Tlmn.re) - mFc[cnt];
                        kahant = pmF[cnt + L * ithread] + kahany;

                        mFc[cnt] = (kahant - pmF[cnt + L * ithread]) - kahany;
                        pmF[cnt + L * ithread] = kahant;
                        cnt++;

                        if (incr > 1) {
                           kahany = (calpha * Tlmn.im) - mFc[cnt];
                           kahant = pmF[cnt + L * ithread] + kahany;

                           mFc[cnt] = (kahant - pmF[cnt + L * ithread]) - kahany;
                           pmF[cnt + L * ithread] = kahant;
                           cnt++;
                        }
                     }
                  }
               }
            }
         }
#pragma omp for
         for (i = 0; i < L; i++) {
            double ssum;
            ssum = 0.0;
            kahansum(nthread, pmF[i + k * L], ssum)
               mF[i] = ssum;
         }

         free(sphi1); free(cphi1);
         free(sPhi); free(cPhi);
         free(sphi2); free(cphi2);
         free(mFc);
      }
      free(pmF);
   }
   //mexPrintf("exit mf\n");

   free(mLinEq);
   free(nLinEq);
}

void calcMFR(double* mFR, double* phi, double* alpha, size_t CS, size_t SS, size_t maxL, size_t N) {
   size_t DF, iSymCombo, incr;
   size_t L;
   double sum = 0.0, iA;
   //Symmetry Assignments
   size_t* mLinEq, * nLinEq;

   if (maxL > RANK_MAX) maxL = RANK_MAX;

   mLinEq = (size_t*)calloc(maxL + 1, sizeof(size_t));
   nLinEq = (size_t*)calloc(maxL + 1, sizeof(size_t));

   iSymCombo = getSymCombo(CS, SS);

   NumLinEq(iSymCombo, CS, maxL, mLinEq);
   NumLinEq(iSymCombo, SS, maxL, nLinEq);

   if (iSymCombo == 4 || iSymCombo == 3) {
      incr = 1;
   }
   else {
      incr = 2;
   }

   L = 0;

   for (size_t l = 0; l <= maxL; l++) {
      for (size_t mu = 1; mu <= mLinEq[l]; mu++) {
         for (size_t nu = 1; nu <= nLinEq[l]; nu++) {
            L += incr;
         }
      }
   }

   double* mF = (double*)calloc(L, sizeof(double));

   {
      double A;
      kahansum(N, alpha[k], A)
         iA = 1.0 / A;
      //for (size_t i = 0; i < N; i++) alpha[i] *= iA;
   }


   DF = 0;

   {
      for (size_t j = 0; j < L; j++) {
         mF[j] = 0.0;
      }

      double* pmF;

#pragma omp parallel default(none) firstprivate(DF, iSymCombo, CS, SS, maxL, L, N, incr, mLinEq, nLinEq, iA) shared(mF, phi, pmF, alpha)
      {
         double* sphi1;
         double* cphi1;
         double* sPhi;
         double* cPhi;
         double* sphi2;
         double* cphi2;
         double* mFc;
         const size_t nthread = omp_get_num_threads();;
         const size_t ithread = omp_get_thread_num();
         size_t  j, l, mu, nu;
         int i;

         sphi1 = (double*)calloc(maxL + 1, sizeof(double));
         cphi1 = (double*)calloc(maxL + 1, sizeof(double));
         sPhi = (double*)calloc(maxL + 1, sizeof(double));
         cPhi = (double*)calloc(maxL + 1, sizeof(double));
         sphi2 = (double*)calloc(maxL + 1, sizeof(double));
         cphi2 = (double*)calloc(maxL + 1, sizeof(double));
         mFc = (double*)malloc(L * sizeof(double));

         for (j = 0; j < L; j++) mFc[j] = 0.0;

#pragma omp single
         {
            pmF = (double*)malloc(L * nthread * sizeof(double));
         }

         //double*rmF = &pmF[L * ithread];

         for (i = 0; i < L; i++) {
            pmF[i + L * ithread] = 0.0;
         }

#pragma omp for// schedule(static)
         for (i = 0; i < N; i++) {

            double calpha;
            calpha = alpha[i] * iA;
            //pre-calculate all the sines and cosines we will need
            for (j = 0; j <= maxL; j++) {
               sphi1[j] = sin((double)j * phi[3 * i]);
               cphi1[j] = cos((double)j * phi[3 * i]);
               sPhi[j] = sin((double)j * phi[1 + 3 * i]);
               cPhi[j] = cos((double)j * phi[1 + 3 * i]);
               sphi2[j] = sin((double)j * phi[2 + 3 * i]);
               cphi2[j] = cos((double)j * phi[2 + 3 * i]);
            }

            // cnt counts the index of the fourier coeff
            size_t cnt;
            cnt = 0;
            for (l = 0; l <= maxL; l++) {
               // mu represents m for the crystal symmetry
               for (mu = 1; mu <= mLinEq[l]; mu++) {
                  //nu represents n for the sample symmetry
                  for (nu = 1; nu <= nLinEq[l]; nu++) {
                     //size_t l, mu, nu;
                     double kahany, kahant;
                     gshe_complex Tlmn;
                     // Tlmn is a structure that contains both the real and imaginary parts of the coeffs, as well as the derivatives
                     Tlmn = symmetricT(DF, iSymCombo, CS, SS, l, mu, nu, sphi1, cphi1, sPhi, cPhi, sphi2, cphi2);

                     kahany = (calpha * Tlmn.re) - mFc[cnt];
                     kahant = pmF[cnt + L * ithread] + kahany;

                     mFc[cnt] = (kahant - pmF[cnt + L * ithread]) - kahany;
                     pmF[cnt + L * ithread] = kahant;
                     cnt++;

                     if (incr > 1) {
                        kahany = (calpha * Tlmn.im) - mFc[cnt];
                        kahant = pmF[cnt + L * ithread] + kahany;

                        mFc[cnt] = (kahant - pmF[cnt + L * ithread]) - kahany;
                        pmF[cnt + L * ithread] = kahant;
                        cnt++;
                     }
                  }
               }
            }
         }
#pragma omp for
         for (i = 0; i < L; i++) {
            double ssum;
            ssum = 0.0;
            kahansum(nthread, pmF[i + k * L], ssum)
               mF[i] = ssum;
         }

         free(sphi1); free(cphi1);
         free(sPhi); free(cPhi);
         free(sphi2); free(cphi2);
         free(mFc);
      }
      free(pmF);
   }

   L = 0;
   for (size_t l = 0; l <= maxL; l++) {
      size_t LL = 0;
      for (size_t mu = 1; mu <= mLinEq[l]; mu++) {
         for (size_t nu = 1; nu <= nLinEq[l]; nu++) {
            LL += incr;
            L += incr;
         }
      }
      double ssum = 0.0;
      kahansum(LL, mF[k + (L - LL)] * mF[k + (L - LL)], ssum);
      mFR[l] = (ssum);
   }
   //mexPrintf("exit mf\n");

   free(mF);
   free(mLinEq);
   free(nLinEq);
}

void calcPHIS(double* E, double* J, double* H, double* phi, double* alpha, double* mF, double* kernel, size_t CS, size_t SS, size_t maxL, size_t N) {
   size_t DF, iSymCombo, incr;
   size_t nvars, L, LL;
   double sum = 0.0, iE, iA;
   //Symmetry Assignments
   size_t* mLinEq, * nLinEq;
   //Gl
   double* dG, * Gl;
   double* mFR;

   if (maxL > RANK_MAX) maxL = RANK_MAX;

   mLinEq = (size_t*)calloc(maxL + 1, sizeof(size_t));
   nLinEq = (size_t*)calloc(maxL + 1, sizeof(size_t));

   iSymCombo = getSymCombo(CS, SS);

   NumLinEq(iSymCombo, CS, maxL, mLinEq);
   NumLinEq(iSymCombo, SS, maxL, nLinEq);

   if (iSymCombo == 4 || iSymCombo == 3) {
      incr = 1;
   }
   else {
      incr = 2;
   }

   L = 0;
   for (size_t l = 0; l <= maxL; l++) {
      //mexPrintf("K(%i) = %.3f\n",l, kernel[l]);
      if (kernel[l] != 0.0) {
         for (size_t mu = 1; mu <= mLinEq[l]; mu++) {
            for (size_t nu = 1; nu <= nLinEq[l]; nu++) {
               L += incr;
            }
         }
      }
   }

   mFR = (double*)calloc(L, sizeof(double));
   L = 0;
   LL = 0;
   for (size_t l = 0; l <= maxL; l++) {
      for (size_t mu = 1; mu <= mLinEq[l]; mu++) {
         for (size_t nu = 1; nu <= nLinEq[l]; nu++) {

            if (kernel[l] != 0.0) {
               mFR[L] = mF[LL];
               if (incr == 2) {
                  mFR[L + 1] = mF[LL + 1];
               }
               L += incr;
            }
            LL += incr;
         }
      }
   }

   //mexPrintf("LL = %i; L = %i\n", LL, L);
   //for (size_t l = 0; l < L; l++) mexPrintf("LR(%i) = %i;\n", l, Lremap[l]);
   //return;
      //

   nvars = N * 3;

   // Allocate memory for each of the Coeffs
   Gl = (double*)calloc(L, sizeof(double));
   dG = (double*)calloc(nvars * L, sizeof(double));

   {
      double A;
      kahansum(N, alpha[k], A)
         iA = 1.0 / A;
      //for (size_t i = 0; i < N; i++) alpha[i] *= iA;
   }


   DF = 0;

   {
      for (size_t j = 0; j < L; j++) {
         Gl[j] = 0.0;
      }

      double* pGl;

#pragma omp parallel default(none) firstprivate(DF, iSymCombo, CS, SS, maxL, L, N, incr, mLinEq, nLinEq, mFR, iA, kernel) shared(Gl, phi, pGl, alpha)
      {
         double* sphi1 = (double*)calloc(maxL + 1, sizeof(double));
         double* cphi1 = (double*)calloc(maxL + 1, sizeof(double));
         double* sPhi = (double*)calloc(maxL + 1, sizeof(double));
         double* cPhi = (double*)calloc(maxL + 1, sizeof(double));
         double* sphi2 = (double*)calloc(maxL + 1, sizeof(double));
         double* cphi2 = (double*)calloc(maxL + 1, sizeof(double));
         double* Glc = (double*)calloc(L, sizeof(double));
         const size_t nthread = omp_get_num_threads();//1;//
         const size_t ithread = omp_get_thread_num();//0;//
         size_t j, l, mu, nu;
         int i;

         for (j = 0; j < L; j++) Glc[j] = 0.0;

#pragma omp single
         {
            pGl = (double*)malloc(L * nthread * sizeof(double));
         }

         double* rGl = &pGl[L * ithread];

         for (i = 0; i < L; i++) {
            pGl[i + L * ithread] = 0.0;
         }

#pragma omp for schedule(static)
         for (i = 0; i < N; i++) {

            double calpha = alpha[i] * iA;
            //pre-calculate all the sines and cosines we will need
            for (j = 0; j <= maxL; j++) {
               sphi1[j] = sin((double)j * phi[3 * i]);
               cphi1[j] = cos((double)j * phi[3 * i]);
               sPhi[j] = sin((double)j * phi[1 + 3 * i]);
               cPhi[j] = cos((double)j * phi[1 + 3 * i]);
               sphi2[j] = sin((double)j * phi[2 + 3 * i]);
               cphi2[j] = cos((double)j * phi[2 + 3 * i]);
            }

            // cnt counts the index of the fourier coeff
            size_t cnt = 0;
            for (l = 0; l <= maxL; l++) {
               double pd_fac = kernel[l];
               if (pd_fac != 0.0) {
                  // mu represents m for the crystal symmetry
                  for (mu = 1; mu <= mLinEq[l]; mu++) {
                     //nu represents n for the sample symmetry
                     for (nu = 1; nu <= nLinEq[l]; nu++) {
                        //size_t l, mu, nu;
                        double kahany, kahant;
                        gshe_complex Tlmn;
                        // Tlmn is a structure that contains both the real and imaginary parts of the coeffs, as well as the derivatives
                        Tlmn = symmetricT(DF, iSymCombo, CS, SS, l, mu, nu, sphi1, cphi1, sPhi, cPhi, sphi2, cphi2);
                        Tlmn = sc_mult(DF, Tlmn, pd_fac);

                        kahany = (calpha * (Tlmn.re - mFR[cnt])) - Glc[cnt];
                        kahant = rGl[cnt] + kahany;

                        Glc[cnt] = (kahant - rGl[cnt]) - kahany;
                        rGl[cnt] = kahant;
                        cnt++;

                        if (incr > 1) {
                           kahany = (calpha * (Tlmn.im - mFR[cnt])) - Glc[cnt];
                           kahant = pGl[cnt + ithread * L] + kahany;

                           Glc[cnt] = (kahant - rGl[cnt]) - kahany;
                           rGl[cnt] = kahant;
                           cnt++;
                        }
                     }
                  }
               }
            }
         }
#pragma omp for
         for (i = 0; i < L; i++) {
            double ssum;
            ssum = 0.0;
            kahansum(nthread, pGl[i + k * L], ssum)
               Gl[i] = ssum;
         }

         free(sphi1); free(cphi1);
         free(sPhi); free(cPhi);
         free(sphi2); free(cphi2);
         free(Glc);
      }
      free(pGl);
   }

   kahansum(L, Gl[k] * Gl[k], E[0]);

   E[0] = sqrt(E[0]);
   iE = 1.0 / E[0];

   DF = 2;

   for (int i = 0; i < nvars; i++) {
      for (int j = i; j < nvars; j++) {
         H[j + i * nvars] = 0.0;
      }
   }

#pragma omp parallel  default(none)  firstprivate(DF, iSymCombo, CS, SS, maxL, L, N, incr, mLinEq, nLinEq, nvars, iE, iA, kernel) shared(Gl, dG, J, H, phi, alpha)
   {
      double* sphi1 = (double*)calloc(maxL + 1, sizeof(double));
      double* cphi1 = (double*)calloc(maxL + 1, sizeof(double));
      double* sPhi = (double*)calloc(maxL + 1, sizeof(double));
      double* cPhi = (double*)calloc(maxL + 1, sizeof(double));
      double* sphi2 = (double*)calloc(maxL + 1, sizeof(double));
      double* cphi2 = (double*)calloc(maxL + 1, sizeof(double));
#pragma omp for
      for (size_t i = 0; i < N; i++) {
         double sum6[6], sum6c[6], calpha;



         gshe_complex Tlmn;

         size_t j, l, mu, nu, cnt;

         calpha = alpha[i] * iA;

         //pre-calculate all the sines and cosines we will need
         for (j = 0; j <= maxL; j++) {
            sphi1[j] = sin((double)j * phi[3 * i]);
            cphi1[j] = cos((double)j * phi[3 * i]);
            sPhi[j] = sin((double)j * phi[1 + 3 * i]);
            cPhi[j] = cos((double)j * phi[1 + 3 * i]);
            sphi2[j] = sin((double)j * phi[2 + 3 * i]);
            cphi2[j] = cos((double)j * phi[2 + 3 * i]);
         }

         for (j = 0; j < 6; j++) {
            sum6[j] = 0.0;
            sum6c[j] = 0.0;
         }

         // cnt counts the index of the fourier coeff
         cnt = 0;

         for (l = 0; l <= maxL; l++) {
            double pd_fac = kernel[l];
            if (pd_fac != 0.0) {
               // mu represents m for the crystal symmetry
               for (mu = 1; mu <= mLinEq[l]; mu++) {
                  //nu represents n for the sample symmetry
                  for (nu = 1; nu <= nLinEq[l]; nu++) {
                     double kahany, kahant;
                     // Tlmn is a structure that contains both the real and imaginary parts of the coeffs, as well as the derivatives
                     Tlmn = symmetricT(DF, iSymCombo, CS, SS, l, mu, nu, sphi1, cphi1, sPhi, cPhi, sphi2, cphi2);
                     Tlmn = sc_mult(DF, Tlmn, pd_fac);

                     dG[cnt + 3 * i * L] = Tlmn.dre[0] * calpha;
                     dG[1 + cnt + 3 * i * L] = Tlmn.dim[0] * calpha;
                     dG[cnt + (3 * i + 1) * L] = Tlmn.dre[1] * calpha;//
                     dG[1 + cnt + (3 * i + 1) * L] = Tlmn.dim[1] * calpha;
                     dG[cnt + (3 * i + 2) * L] = Tlmn.dre[2] * calpha;//
                     dG[1 + cnt + (3 * i + 2) * L] = Tlmn.dim[2] * calpha;

                     //H[n, m] += ddG/ddphi

                     for (j = 0; j < 6; j++) {
                        kahany = Gl[cnt] * Tlmn.ddre[j] - sum6c[j];
                        kahant = sum6[j] + kahany;

                        sum6c[j] = (kahant - sum6[j]) - kahany;
                        sum6[j] = kahant;
                     }
                     cnt++;

                     for (j = 0; j < 6; j++) {
                        kahany = Gl[cnt] * Tlmn.ddim[j] - sum6c[j];
                        kahant = sum6[j] + kahany;

                        sum6c[j] = (kahant - sum6[j]) - kahany;
                        sum6[j] = kahant;
                     }
                     cnt++;

                  }
               }
            }
         }
         // mm = [0 0 0 1 1 2]
         // nn = [0 1 2 1 2 2]
         //(nn + 3 * i) + (mm + 3 * i) * nvars
         H[(3 * i) + (3 * i) * nvars] = sum6[0] * calpha;
         H[(1 + 3 * i) + (3 * i) * nvars] = sum6[1] * calpha;
         H[(2 + 3 * i) + (3 * i) * nvars] = sum6[2] * calpha;
         H[(1 + 3 * i) + (1 + 3 * i) * nvars] = sum6[3] * calpha;
         H[(2 + 3 * i) + (1 + 3 * i) * nvars] = sum6[4] * calpha;
         H[(2 + 3 * i) + (2 + 3 * i) * nvars] = sum6[5] * calpha;

      }
      free(sphi1); free(cphi1);
      free(sPhi); free(cPhi);
      free(sphi2); free(cphi2);
   }

   {
      int i;
#pragma omp parallel for default(none)  firstprivate(L, nvars, iE) shared(Gl, dG, J)
      for (i = 0; i < nvars; i++) {
         double ssum;
         kahansum(L, Gl[k] * dG[k + i * L], ssum);
         J[i] = ssum * iE;
      }
   }


   {
      //H = H - J*J'
      char* uplo = "L";
      ptrdiff_t blas_nvars, incx;
      double beta;
      blas_nvars = nvars;
      incx = 1;
      beta = -1.0;
      dsyr(uplo, &blas_nvars, &beta, J, &incx, H, &blas_nvars);
   }

   {
      //H = H + dG*dG'
      char* uplo = "L";
      char* trans = "T";
      ptrdiff_t blas_nvars, blas_l;
      double beta = 1.0;
      blas_nvars = nvars;
      blas_l = L;
      dsyrk(uplo, trans, &blas_nvars, &blas_l, &beta, dG, &blas_l, &beta, H, &blas_nvars);
   }

#pragma omp parallel default(none) firstprivate(nvars, iE) shared(H)
   {
#pragma omp for //collapse(2)
      for (size_t i = 0; i < nvars; i++) {
         for (size_t j = i; j < nvars; j++) {
            H[j + i * nvars] *= iE;
            H[i + j * nvars] = H[j + i * nvars];
         }
      }

   }
   /* Put T into the output for matlab and set the array sizes */

   free(Gl);
   free(dG);

   free(mFR);

   free(mLinEq);
   free(nLinEq);
}

//void calcPD(double* pd, double* mF, double *T, size_t CS, size_t SS, size_t maxL, size_t N){
//
//
//}

void calcFULL(double* E, double* J, double* H, double* phi, double* mF, double* kernel, size_t CS, size_t SS, size_t maxL, size_t N) {
   size_t DF, iSymCombo, incr;
   size_t LL;
   size_t nvars, L;
   double sum = 0.0, iE, iA;
   //Symmetry Assignments
   size_t* mLinEq, * nLinEq;
   //Gl
   double* dG, * Gl, * mFR;

   if (maxL > RANK_MAX) maxL = RANK_MAX;

   mLinEq = (size_t*)calloc(maxL + 1, sizeof(size_t));
   nLinEq = (size_t*)calloc(maxL + 1, sizeof(size_t));

   iSymCombo = getSymCombo(CS, SS);

   NumLinEq(iSymCombo, CS, maxL, mLinEq);
   NumLinEq(iSymCombo, SS, maxL, nLinEq);

   if (iSymCombo == 4 || iSymCombo == 3) {
      incr = 1;
   }
   else {
      incr = 2;
   }

   L = 0;
   for (size_t l = 0; l <= maxL; l++) {
      for (size_t mu = 1; mu <= mLinEq[l]; mu++) {
         for (size_t nu = 1; nu <= nLinEq[l]; nu++) {
            L += incr;
         }
      }
   }

   mFR = (double*)calloc(L, sizeof(double));
   L = 0;
   LL = 0;
   for (size_t l = 0; l <= maxL; l++) {
      for (size_t mu = 1; mu <= mLinEq[l]; mu++) {
         for (size_t nu = 1; nu <= nLinEq[l]; nu++) {

            if (kernel[l] != 0.0) {
               mFR[L] = mF[LL];
               if (incr == 2) {
                  mFR[L + 1] = mF[LL + 1];
               }
               L += incr;
            }
            LL += incr;
         }
      }
   }

   nvars = 4 * N;

   // Allocate memory for each of the Coeffs
   Gl = (double*)calloc(L, sizeof(double));
   dG = (double*)calloc(nvars * L, sizeof(double));


   {
      double A;
      kahansum(N, phi[3 + 4 * k], A);
      iA = 1.0 / A;
   }


   DF = 0;

   {

      double* pGl;

#pragma omp parallel default(none) firstprivate(DF, iSymCombo, CS, SS, maxL, L, N, incr, mLinEq, nLinEq, mFR, iA) shared(Gl, phi, pGl, kernel)
      {
         double* sphi1 = (double*)calloc(maxL + 1, sizeof(double));
         double* cphi1 = (double*)calloc(maxL + 1, sizeof(double));
         double* sPhi = (double*)calloc(maxL + 1, sizeof(double));
         double* cPhi = (double*)calloc(maxL + 1, sizeof(double));
         double* sphi2 = (double*)calloc(maxL + 1, sizeof(double));
         double* cphi2 = (double*)calloc(maxL + 1, sizeof(double));
         double* Glc = (double*)malloc(L * sizeof(double));
         const size_t nthread = omp_get_num_threads();//1;//
         const size_t ithread = omp_get_thread_num();//0;//
         size_t j, l, mu, nu;
         int i;

         for (j = 0; j < L; j++) Glc[j] = 0.0;

#pragma omp single
         {
            pGl = (double*)malloc(L * nthread * sizeof(double));
         }

         //double* rGl = &pGl[L * ithread];

         for (i = 0; i < L; i++) {
            pGl[i + L * ithread] = 0.0;
         }

#pragma omp for schedule(static)
         for (i = 0; i < N; i++) {

            double calpha = phi[3 + 4 * i] * iA;
            //pre-calculate all the sines and cosines we will need
            for (j = 0; j <= maxL; j++) {
               sphi1[j] = sin((double)j * phi[4 * i]);
               cphi1[j] = cos((double)j * phi[4 * i]);
               sPhi[j] = sin((double)j * phi[1 + 4 * i]);
               cPhi[j] = cos((double)j * phi[1 + 4 * i]);
               sphi2[j] = sin((double)j * phi[2 + 4 * i]);
               cphi2[j] = cos((double)j * phi[2 + 4 * i]);
            }

            // cnt counts the index of the fourier coeff
            size_t cnt = 0;
            //for  schedule(static) default(none)  firstprivate(L, incr, calpha, DF, iSymCombo, CS, SS, sphi1, cphi1, sPhi, cPhi, sphi2, cphi2) shared(Gl, Glc, mF, lmn)
            for (l = 0; l <= maxL; l++) {
               double pd_fac = kernel[l];
               if (pd_fac != 0.0) {
                  // mu represents m for the crystal symmetry
                  for (mu = 1; mu <= mLinEq[l]; mu++) {
                     //nu represents n for the sample symmetry
                     for (nu = 1; nu <= nLinEq[l]; nu++) {
                        //size_t l, mu, nu;
                        double kahany, kahant;
                        gshe_complex Tlmn;
                        // Tlmn is a structure that contains both the real and imaginary parts of the coeffs, as well as the derivatives
                        Tlmn = symmetricT(DF, iSymCombo, CS, SS, l, mu, nu, sphi1, cphi1, sPhi, cPhi, sphi2, cphi2);
                        Tlmn = sc_mult(DF, Tlmn, pd_fac);

                        kahany = (calpha * (Tlmn.re - mFR[cnt])) - Glc[cnt];
                        kahant = pGl[cnt + L * ithread] + kahany;

                        Glc[cnt] = (kahant - pGl[cnt + L * ithread]) - kahany;
                        pGl[cnt + L * ithread] = kahant;
                        cnt++;

                        if (incr > 1) {
                           kahany = (calpha * (Tlmn.im - mFR[cnt])) - Glc[cnt];
                           kahant = pGl[cnt + ithread * L] + kahany;

                           Glc[cnt] = (kahant - pGl[cnt + L * ithread]) - kahany;
                           pGl[cnt + L * ithread] = kahant;
                           cnt++;
                        }
                     }
                  }
               }
            }
         }
#pragma omp for
         for (i = 0; i < L; i++) {
            double ssum;
            ssum = 0.0;
            kahansum(nthread, pGl[i + k * L], ssum);
            Gl[i] = ssum;
         }

         free(sphi1); free(cphi1);
         free(sPhi); free(cPhi);
         free(sphi2); free(cphi2);
         free(Glc);
      }
      free(pGl);
   }

   kahansum(L, Gl[k] * Gl[k], E[0]);

   E[0] = sqrt(E[0]);
   iE = 1.0 / E[0];

   DF = 2;

   for (int i = 0; i < nvars; i++) {
      for (int j = i; j < nvars; j++) {
         H[j + i * nvars] = 0.0;
      }
   }

#pragma omp parallel  default(none)  firstprivate(DF, iSymCombo, CS, SS, maxL, L, N, incr, mLinEq, nLinEq, nvars, iE, iA, mFR) shared(Gl, dG, kernel, J, H, phi)
   {
      double* sphi1 = (double*)calloc(maxL + 1, sizeof(double));
      double* cphi1 = (double*)calloc(maxL + 1, sizeof(double));
      double* sPhi = (double*)calloc(maxL + 1, sizeof(double));
      double* cPhi = (double*)calloc(maxL + 1, sizeof(double));
      double* sphi2 = (double*)calloc(maxL + 1, sizeof(double));
      double* cphi2 = (double*)calloc(maxL + 1, sizeof(double));
      int i;
      double sum6[6], sum6c[6], tempTA[3], tempTAc[3];

#pragma omp for
      for (i = 0; i < N; i++) {
         gshe_complex Tlmn;

         size_t j, l, mu, nu, cnt;

         //pre-calculate all the sines and cosines we will need
         for (j = 0; j <= maxL; j++) {
            sphi1[j] = sin((double)j * phi[4 * i]);
            cphi1[j] = cos((double)j * phi[4 * i]);
            sPhi[j] = sin((double)j * phi[1 + 4 * i]);
            cPhi[j] = cos((double)j * phi[1 + 4 * i]);
            sphi2[j] = sin((double)j * phi[2 + 4 * i]);
            cphi2[j] = cos((double)j * phi[2 + 4 * i]);
         }

         double calpha = phi[3 + 4 * i] * iA;

         for (j = 0; j < 6; j++) {
            sum6[j] = 0.0;
            sum6c[j] = 0.0;
         }

         for (j = 0; j < 3; j++) {
            tempTA[j] = 0.0;
            tempTAc[j] = 0.0;
         }
         // cnt counts the index of the fourier coeff
         cnt = 0;

         for (l = 0; l <= maxL; l++) {
            double pd_fac = kernel[l]; //sqrt(l*2 + 1);
            if (pd_fac != 0.0) {
               // mu represents m for the crystal symmetry
               for (mu = 1; mu <= mLinEq[l]; mu++) {
                  //nu represents n for the sample symmetry
                  for (nu = 1; nu <= nLinEq[l]; nu++) {
                     double kahany, kahant;
                     // Tlmn is a structure that contains both the real and imaginary parts of the coeffs, as well as the derivatives
                     Tlmn = symmetricT(DF, iSymCombo, CS, SS, l, mu, nu, sphi1, cphi1, sPhi, cPhi, sphi2, cphi2);
                     Tlmn = sc_mult(DF, Tlmn, pd_fac);


                     dG[cnt + 4 * i * L] = Tlmn.dre[0] * calpha;
                     dG[1 + cnt + 4 * i * L] = Tlmn.dim[0] * calpha;
                     dG[cnt + (4 * i + 1) * L] = Tlmn.dre[1] * calpha;//
                     dG[1 + cnt + (4 * i + 1) * L] = Tlmn.dim[1] * calpha;
                     dG[cnt + (4 * i + 2) * L] = Tlmn.dre[2] * calpha;//
                     dG[1 + cnt + (4 * i + 2) * L] = Tlmn.dim[2] * calpha;
                     dG[cnt + (4 * i + 3) * L] = (Tlmn.re - mFR[cnt] - Gl[cnt]) * iA;
                     dG[1 + cnt + (4 * i + 3) * L] = (Tlmn.im - mFR[cnt + 1] - Gl[cnt + 1]) * iA;


                     //H[n, m] += ddG/ddphi
                     for (j = 0; j < 6; j++) {
                        kahany = Gl[cnt] * Tlmn.ddre[j] - sum6c[j];
                        kahant = sum6[j] + kahany;

                        sum6c[j] = (kahant - sum6[j]) - kahany;
                        sum6[j] = kahant;
                     }

                     //H[n, m] += ddG/(da*dphi)
                     for (j = 0; j < 3; j++) {
                        kahany = Gl[cnt] * Tlmn.dre[j] * (1.0 - calpha) - tempTAc[j];//
                        kahant = tempTA[j] + kahany;

                        tempTAc[j] = (kahant - tempTA[j]) - kahany;
                        tempTA[j] = kahant;
                     }
                     cnt++;


                     //H[n, m] += ddG/ddphi
                     for (j = 0; j < 6; j++) {
                        kahany = Gl[cnt] * Tlmn.ddim[j] - sum6c[j];
                        kahant = sum6[j] + kahany;

                        sum6c[j] = (kahant - sum6[j]) - kahany;
                        sum6[j] = kahant;
                     }

                     //H[n, m] += ddG/(da*dphi)
                     for (j = 0; j < 3; j++) {
                        kahany = Gl[cnt] * Tlmn.dim[j] * (1.0 - calpha) - tempTAc[j];//
                        kahant = tempTA[j] + kahany;

                        tempTAc[j] = (kahant - tempTA[j]) - kahany;
                        tempTA[j] = kahant;
                     }
                     cnt++;

                  }
               }
            }
         }
         // mm = [0 0 0 1 1 2]
         // nn = [0 1 2 1 2 2]
         // (nn + 4 * i) + (mm + 4 * i) * nvars
         // H[n, m] += ddG/ddphi
         H[(4 * i) + (4 * i) * nvars] = sum6[0] * calpha;
         H[(1 + 4 * i) + (4 * i) * nvars] = sum6[1] * calpha;
         H[(2 + 4 * i) + (4 * i) * nvars] = sum6[2] * calpha;
         H[(3 + 4 * i) + (4 * i) * nvars] = tempTA[0] * iA;// H[n, m] = ddG/(da*dphi)
         H[(1 + 4 * i) + (1 + 4 * i) * nvars] = sum6[3] * calpha;
         H[(2 + 4 * i) + (1 + 4 * i) * nvars] = sum6[4] * calpha;
         H[(3 + 4 * i) + (1 + 4 * i) * nvars] = tempTA[1] * iA;// H[n, m] = ddG/(da*dphi)
         H[(2 + 4 * i) + (2 + 4 * i) * nvars] = sum6[5] * calpha;
         H[(3 + 4 * i) + (2 + 4 * i) * nvars] = tempTA[2] * iA;// H[n, m] = ddG/(da*dphi)

      }
      free(sphi1); free(cphi1);
      free(sPhi); free(cPhi);
      free(sphi2); free(cphi2);
   }

#pragma omp parallel default(none)  firstprivate(L, nvars) shared(Gl, dG, J)
   {
      int i;
#pragma omp for 
      for (i = 0; i < nvars; i++) {
         double ssum;
         kahansum(L, Gl[k] * dG[k + i * L], ssum)
            J[i] = ssum;
      }
   }

#pragma omp parallel default(none) firstprivate(J, nvars, iA) shared(H)
   {
      int i;
#pragma omp for
      for (i = 0; i < nvars; i++) {
         const double temp = J[i] * iA;
         for (int j = (i + 7 - i % 4); j < nvars; j += 4) {
            H[j + i * nvars] -= temp;
         }
         if ((i + 1) % 4 == 0) {
            H[i + i * nvars] -= 2.0 * (J[i]) * iA;
            for (int j = i + 1; j < nvars; j++) {
               H[j + i * nvars] -= J[j] * iA;
            }
         }
      }
      //#pragma omp for
      //for (i = 3; i < nvars; i+=4) {
      //    H[i + i*nvars] -= (J[i] + J[i])*iA;
      //    for (int j = i+1; j < nvars; j++){
      //        H[j + i*nvars] -= J[j]*iA;
      //    }
      //}
   }

#pragma omp parallel default(none) firstprivate(nvars, iE) shared(J)
   {
      int i;
#pragma omp for
      for (i = 0; i < nvars; i++) {
         J[i] *= iE;
      }
   }


   {
      //H = H - J*J'
      char* uplo = "L";
      ptrdiff_t blas_nvars, incx;
      double beta;
      blas_nvars = nvars;
      incx = 1;
      beta = -1.0;
      dsyr(uplo, &blas_nvars, &beta, J, &incx, H, &blas_nvars);
   }

   {
      //H = H + dG*dG'
      char* uplo = "L";
      char* trans = "T";
      ptrdiff_t blas_nvars, blas_l;
      double beta = 1.0;
      blas_nvars = nvars;
      blas_l = L;
      dsyrk(uplo, trans, &blas_nvars, &blas_l, &beta, dG, &blas_l, &beta, H, &blas_nvars);
   }



#pragma omp parallel default(none) firstprivate(nvars, iE) shared(H)
   {
      //H = H/E
      size_t i;
#pragma omp for //collapse(2)
      for (i = 0; i < nvars; i++) {
         for (size_t j = i; j < nvars; j++) {
            H[j + i * nvars] *= iE;
            H[i + j * nvars] = H[j + i * nvars];
         }
      }

   }
   /* Put T into the output for matlab and set the array sizes */

   free(Gl);
   free(dG);

   free(mFR);

   free(mLinEq);
   free(nLinEq);

}

void calcSTD(double* E, double* J, double* H, double* phi, size_t N, double K) {
   double A, A2, iA, iA2, iA3, iA4, sA2, aE, iN, iE, JJ1, JJ2, HH1, HH2, HH3;
   //double *alpha;
   size_t nvars = 4 * N;

   iN = 1.0 / ((double)N);
   kahansum(N, phi[3 + 4 * k], A)
      iA = 1.0 / A;


   kahansum(N, (phi[3 + 4 * k] * iA - iN) * (phi[3 + 4 * k] * iA - iN), aE)

      aE = sqrt(aE);
   iE = 1.0 / aE;
   E[0] += K * aE;

   if (aE < 1e-14) iE = 1.0;

   A2 = A * A;
   iA2 = 1.0 / A2;
   iA3 = pow(A, -3.0);
   iA4 = pow(A, -4.0);

   kahansum(N, phi[3 + 4 * k] * phi[3 + 4 * k], sA2)


      JJ2 = -iE * iA3 * sA2;
   JJ1 = iE * iA2;
   kahansum(N, ((1.5 * phi[3 + 4 * k] * iA - iN) * iA * phi[3 + 4 * k]), HH1)
      HH1 += iN;
   HH1 *= 2.0 * iA2;

   HH1 *= iE;
   HH2 = -iE * 2.0 * iA3;

   HH3 = -iE * JJ1 * JJ1;

   HH1 -= iE * JJ2 * JJ2;
   HH2 -= iE * JJ1 * JJ2;

   //mexPrintf("std omp %.4f %.4f %.4f %.4f\n", K, HH1, HH2, HH3);




#pragma omp parallel default(none) firstprivate(N, phi, nvars, K, iE, iA2, JJ1, JJ2, HH1, HH2, HH3) shared(H, J)
   {

      size_t j;
#pragma omp for
      for (size_t i = 3; i < nvars; i += 4) {
         J[i] += K * (JJ1 * phi[i] + JJ2);
         H[i + i * nvars] += K * iE * iA2;
      }

#pragma omp for collapse(2)
      for (size_t i = 3; i < nvars; i += 4) {
         for (j = 3; j < nvars; j += 4) {
            H[j + i * nvars] += K * (HH1 + (phi[i] + phi[j]) * HH2 + phi[i] * phi[j] * HH3);
         }
      }
   }
   //mexPrintf("out std omp\n");

   //free(Sk);
   //free(dSk);
   //free(ddSk);
}

void calcE(double* E, double* phi, double* mF, size_t CS, size_t SS, size_t maxL, size_t N) {
   size_t DF, iSymCombo, incr;
   size_t nvars, L;
   double sum = 0.0, iA;
   //Symmetry Assignments
   size_t* mLinEq, * nLinEq;
   //Gl
   double* Gl;

   if (maxL > RANK_MAX) maxL = RANK_MAX;

   mLinEq = (size_t*)calloc(maxL + 1, sizeof(size_t));
   nLinEq = (size_t*)calloc(maxL + 1, sizeof(size_t));

   iSymCombo = getSymCombo(CS, SS);

   NumLinEq(iSymCombo, CS, maxL, mLinEq);
   NumLinEq(iSymCombo, SS, maxL, nLinEq);

   if (iSymCombo == 4 || iSymCombo == 3) {
      incr = 1;
   }
   else {
      incr = 2;
   }

   L = 0;
   for (size_t l = 0; l <= maxL; l++) {
      for (size_t mu = 1; mu <= mLinEq[l]; mu++) {
         for (size_t nu = 1; nu <= nLinEq[l]; nu++) {
            L += incr;
         }
      }
   }

   nvars = 4 * N;

   // Allocate memory for each of the Coeffs
   Gl = (double*)calloc(L, sizeof(double));


   {
      double A;
      kahansum(N, phi[3 + 4 * k], A)
         iA = 1.0 / A;
   }


   DF = 0;

   {

      double* pGl;

#pragma omp parallel default(none) firstprivate(DF, iSymCombo, CS, SS, maxL, L, N, incr, mLinEq, nLinEq, mF, iA) shared(Gl, phi, pGl)
      {
         double* sphi1 = (double*)calloc(maxL + 1, sizeof(double));
         double* cphi1 = (double*)calloc(maxL + 1, sizeof(double));
         double* sPhi = (double*)calloc(maxL + 1, sizeof(double));
         double* cPhi = (double*)calloc(maxL + 1, sizeof(double));
         double* sphi2 = (double*)calloc(maxL + 1, sizeof(double));
         double* cphi2 = (double*)calloc(maxL + 1, sizeof(double));
         double* Glc = (double*)malloc(L * sizeof(double));
         const size_t nthread = omp_get_num_threads();//1;//
         const size_t ithread = omp_get_thread_num();//0;//
         size_t j, l, mu, nu;
         int i;

         for (j = 0; j < L; j++) Glc[j] = 0.0;

#pragma omp single
         {
            pGl = (double*)malloc(L * nthread * sizeof(double));
         }

         double* rGl = &pGl[L * ithread];

         for (i = 0; i < L; i++) {
            rGl[i] = 0.0;
         }

#pragma omp for schedule(static)
         for (i = 0; i < N; i++) {

            double calpha = phi[3 + 4 * i] * iA;
            //pre-calculate all the sines and cosines we will need
            for (j = 0; j <= maxL; j++) {
               sphi1[j] = sin((double)j * phi[4 * i]);
               cphi1[j] = cos((double)j * phi[4 * i]);
               sPhi[j] = sin((double)j * phi[1 + 4 * i]);
               cPhi[j] = cos((double)j * phi[1 + 4 * i]);
               sphi2[j] = sin((double)j * phi[2 + 4 * i]);
               cphi2[j] = cos((double)j * phi[2 + 4 * i]);
            }

            // cnt counts the index of the fourier coeff
            size_t cnt = 0;
            //for  schedule(static) default(none)  firstprivate(L, incr, calpha, DF, iSymCombo, CS, SS, sphi1, cphi1, sPhi, cPhi, sphi2, cphi2) shared(Gl, Glc, mF, lmn)
            for (l = 0; l <= maxL; l++) {
               // mu represents m for the crystal symmetry
               for (mu = 1; mu <= mLinEq[l]; mu++) {
                  //nu represents n for the sample symmetry
                  for (nu = 1; nu <= nLinEq[l]; nu++) {
                     //size_t l, mu, nu;
                     double kahany, kahant;
                     gshe_complex Tlmn;
                     // Tlmn is a structure that contains both the real and imaginary parts of the coeffs, as well as the derivatives
                     Tlmn = symmetricT(DF, iSymCombo, CS, SS, l, mu, nu, sphi1, cphi1, sPhi, cPhi, sphi2, cphi2);

                     kahany = (calpha * (Tlmn.re - mF[cnt])) - Glc[cnt];
                     kahant = rGl[cnt] + kahany;

                     Glc[cnt] = (kahant - rGl[cnt]) - kahany;
                     rGl[cnt] = kahant;
                     cnt++;

                     if (incr > 1) {
                        kahany = (calpha * (Tlmn.im - mF[cnt])) - Glc[cnt];
                        kahant = pGl[cnt + ithread * L] + kahany;

                        Glc[cnt] = (kahant - rGl[cnt]) - kahany;
                        rGl[cnt] = kahant;
                        cnt++;
                     }
                  }
               }
            }
         }
#pragma omp for
         for (i = 0; i < L; i++) {
            double ssum;
            ssum = 0.0;
            kahansum(nthread, pGl[i + k * L], ssum)
               Gl[i] = ssum;
         }

         free(sphi1); free(cphi1);
         free(sPhi); free(cPhi);
         free(sphi2); free(cphi2);
         free(Glc);
      }
      free(pGl);
   }

   kahansum(L, Gl[k] * Gl[k], E[0])

      E[0] = sqrt(E[0]);

   free(Gl);

   free(mLinEq);
   free(nLinEq);

   return;
}

void calcESTD(double* E, double* phi, size_t N, double K) {
   double A, iA, aE, iN;
   //double *alpha;
   size_t nvars = 4 * N;

   iN = 1.0 / ((double)N);
   kahansum(N, phi[3 + 4 * k], A)
      iA = 1.0 / A;
   kahansum(N, (phi[3 + 4 * k] * iA - iN) * (phi[3 + 4 * k] * iA - iN), aE)
      aE = sqrt(aE);
   E[0] += K * aE;
   //mexPrintf("out std omp\n");
}

void calcDIFFPHI(double* E, double* J, double* H, double* phi, double* alpha, double* mF, double* kernel, const size_t CS, const size_t SS, size_t const maxL, size_t const N, const size_t nvars, const size_t Z, double r) {
   size_t DF, iSymCombo, incr;
   size_t  L;
   double sum = 0.0, iA;
   //Symmetry Assignments
   size_t* mLinEq, * nLinEq;



   mLinEq = (size_t*)calloc(maxL + 1, sizeof(size_t));
   nLinEq = (size_t*)calloc(maxL + 1, sizeof(size_t));

   iSymCombo = getSymCombo(CS, SS);

   NumLinEq(iSymCombo, CS, maxL, mLinEq);
   NumLinEq(iSymCombo, SS, maxL, nLinEq);

   if (iSymCombo == 4 || iSymCombo == 3) {
      incr = 1;
   }
   else {
      incr = 2;
   }

   L = 0;
   for (size_t l = 0; l <= maxL; l++) {
      if (kernel[l] > 0.0) {
         for (size_t mu = 1; mu <= mLinEq[l]; mu++) {
            for (size_t nu = 1; nu <= nLinEq[l]; nu++) {
               L += incr;
            }
         }
      }
   }

   double* EE = (double*)calloc(Z, sizeof(double));
   double* mT = (double*)calloc(L, sizeof(double));
   double* dG = (double*)calloc(3 * N * L, sizeof(double));
   double* ddG = (double*)calloc(6 * N * L, sizeof(double));

   {
      double A;
      kahansum(N, alpha[k], A)
         iA = 1.0 / A;
   }

   DF = 2;



   for (size_t j = 0; j < L; j++) {
      mT[j] = 0.0;
   }

   double* pmT;

#pragma omp parallel default(none) firstprivate(DF, iSymCombo, CS, SS, maxL, L, N, incr, mLinEq, nLinEq, iA) shared(phi, alpha, kernel, mT, dG, ddG, pmT)
   {
      const size_t nthread = omp_get_num_threads();
      const size_t ithread = omp_get_thread_num();
      size_t  i, j, l, mu, nu;
      
      
      
      double* sphi1 = (double*)calloc(maxL + 1, sizeof(double));
      double* cphi1 = (double*)calloc(maxL + 1, sizeof(double));
      double* sPhi = (double*)calloc(maxL + 1, sizeof(double));
      double* cPhi = (double*)calloc(maxL + 1, sizeof(double));
      double* sphi2 = (double*)calloc(maxL + 1, sizeof(double));
      double* cphi2 = (double*)calloc(maxL + 1, sizeof(double));
      double* mTc = (double*)malloc(L * sizeof(double));

      for (j = 0; j < L; j++) mTc[j] = 0.0;

#pragma omp single
      {
         pmT = (double*)malloc(L * nthread * sizeof(double));
      }

      for (i = 0; i < L; i++) pmT[i + L * ithread] = 0.0;

#pragma omp for// schedule(static)
      for (i = 0; i < N; i++) {

         double calpha;
         calpha = alpha[i] * iA;
         //pre-calculate all the sines and cosines we will need
         for (j = 0; j <= maxL; j++) {
            sphi1[j] = sin((double)j * phi[3 * i]);
            cphi1[j] = cos((double)j * phi[3 * i]);
            sPhi[j] = sin((double)j * phi[1 + 3 * i]);
            cPhi[j] = cos((double)j * phi[1 + 3 * i]);
            sphi2[j] = sin((double)j * phi[2 + 3 * i]);
            cphi2[j] = cos((double)j * phi[2 + 3 * i]);
         }

         // cnt counts the index of the fourier coeff
         size_t cnt;
         cnt = 0;
         for (l = 0; l <= maxL; l++) {
            double pd_fac = kernel[l];
            if (pd_fac != 0.0) {
               // mu represents m for the crystal symmetry
               for (mu = 1; mu <= mLinEq[l]; mu++) {
                  //nu represents n for the sample symmetry
                  for (nu = 1; nu <= nLinEq[l]; nu++) {
                     //size_t l, mu, nu;
                     double kahany, kahant;
                     gshe_complex Tlmn;
                     // Tlmn is a structure that contains both the real and imaginary parts of the coeffs, as well as the derivatives
                     Tlmn = symmetricT(DF, iSymCombo, CS, SS, l, mu, nu, sphi1, cphi1, sPhi, cPhi, sphi2, cphi2);
                     Tlmn = sc_mult(DF, Tlmn, pd_fac);
                     kahany = (calpha * Tlmn.re) - mTc[cnt];
                     kahant = pmT[cnt + L * ithread] + kahany;

                     mTc[cnt] = (kahant - pmT[cnt + L * ithread]) - kahany;
                     pmT[cnt + L * ithread] = kahant;

                     //dG
                     dG[cnt + 3 * i * L] = Tlmn.dre[0] * calpha;
                     dG[1 + cnt + 3 * i * L] = Tlmn.dim[0] * calpha;
                     dG[cnt + (3 * i + 1) * L] = Tlmn.dre[1] * calpha;//
                     dG[1 + cnt + (3 * i + 1) * L] = Tlmn.dim[1] * calpha;
                     dG[cnt + (3 * i + 2) * L] = Tlmn.dre[2] * calpha;//
                     dG[1 + cnt + (3 * i + 2) * L] = Tlmn.dim[2] * calpha;

                     //ddG
                     ddG[cnt + 6 * i * L] = Tlmn.ddre[0] * calpha;
                     ddG[1 + cnt + 6 * i * L] = Tlmn.ddim[0] * calpha;
                     ddG[cnt + (6 * i + 1) * L] = Tlmn.ddre[1] * calpha;//
                     ddG[1 + cnt + (6 * i + 1) * L] = Tlmn.ddim[1] * calpha;
                     ddG[cnt + (6 * i + 2) * L] = Tlmn.ddre[2] * calpha;//
                     ddG[1 + cnt + (6 * i + 2) * L] = Tlmn.ddim[2] * calpha;
                     ddG[cnt + (6 * i + 3) * L] = Tlmn.ddre[3] * calpha;
                     ddG[1 + cnt + (6 * i + 3) * L] = Tlmn.ddim[3] * calpha;
                     ddG[cnt + (6 * i + 4) * L] = Tlmn.ddre[4] * calpha;//
                     ddG[1 + cnt + (6 * i + 4) * L] = Tlmn.ddim[4] * calpha;
                     ddG[cnt + (6 * i + 5) * L] = Tlmn.ddre[5] * calpha;//
                     ddG[1 + cnt + (6 * i + 5) * L] = Tlmn.ddim[5] * calpha;

                     cnt++;

                     if (incr > 1) {
                        kahany = (calpha * Tlmn.im) - mTc[cnt];
                        kahant = pmT[cnt + L * ithread] + kahany;

                        mTc[cnt] = (kahant - pmT[cnt + L * ithread]) - kahany;
                        pmT[cnt + L * ithread] = kahant;
                        cnt++;
                     }
                  }
               }
            }
         }
      }
#pragma omp for
      for (i = 0; i < L; i++) {
         double ssum;
         ssum = 0.0;
         kahansum(nthread, pmT[i + k * L], ssum);
         mT[i] = ssum;
      }

      free(sphi1); free(cphi1);
      free(sPhi); free(cPhi);
      free(sphi2); free(cphi2);
      free(mTc);
   }
   free(pmT);


   double* Gz = (double*)calloc(Z * L, sizeof(double));
   double* Jz = (double*)calloc(Z * nvars, sizeof(double));
   double* Er1 = (double*)calloc(Z, sizeof(double));
   double* ddGz = (double*)calloc(6 * N * Z, sizeof(double));
   double* ddG2z = (double*)calloc(6 * N, sizeof(double));



#pragma omp parallel for collapse(2)
   for (size_t i = 0; i < Z; i++) {
      for (size_t j = 0; j < L; j++) {
         Gz[j + i * L] = mT[j] - mF[j + i * L];
      }
   }

#pragma omp parallel for
   for (size_t z = 0; z < Z; z++) {
      double ssum = 0.0;
      kahansum(L, Gz[k + z * L] * Gz[k + z * L], ssum);
      EE[z] = sqrt(ssum);
      Er1[z] = pow(ssum, (r - 1.0) * 0.5);
   }

   {

      char* transa;
      char* transb;
      char* uplo;
      ptrdiff_t blas_nvars, blas_l, blas_z, blas_one, blas_dvars;
      double zero, one, ssum1, HH2, E1r1, HH1;
      one = 1.0;
      zero = 0.0;
      blas_nvars = nvars;
      blas_dvars = 6 * N;
      blas_l = L;
      blas_z = Z;
      blas_one = 1;
      transa = "t";
      transb = "n";
      uplo = "L";

      kahansum(Z, pow(EE[k], r), ssum1);
      kahansum(Z, pow(EE[k], r - 2.0), HH2);

      E1r1 = pow(ssum1, (1.0 / r) - 1.0);

      E[0] = pow(ssum1, (1.0 / r));


#pragma omp parallel for
      for (size_t z = 0; z < Z; z++) {
         double iE = 1.0 / EE[z];

         for (size_t j = 0; j < L; j++) Gz[j + z * L] *= iE;
      }

      //Jz = dG*Gz
      dgemm(transa, transb, &blas_nvars, &blas_z, &blas_l, &one, dG, &blas_l, Gz, &blas_l, &zero, Jz, &blas_nvars);

      //J = sum(Ez^r)^(1/r-1)*Jz*Ez^(r-1)
      dgemv(transb, &blas_nvars, &blas_z, &E1r1, Jz, &blas_nvars, Er1, &blas_one, &zero, J, &blas_one);

      //ddGz = ddG*Gz
      dgemm(transa, transb, &blas_dvars, &blas_z, &blas_l, &one, ddG, &blas_l, Gz, &blas_l, &zero, ddGz, &blas_dvars);

      //ddGZ2 = sum(Ez^r)^(1/r-1)*ddGz*Ez^(r-1)
      dgemv(transb, &blas_dvars, &blas_z, &E1r1, ddGz, &blas_dvars, Er1, &blas_one, &zero, ddG2z, &blas_one);

      HH1 = E1r1 * HH2;

      // Hz = dGz'*dGz
      dsyrk(uplo, transa, &blas_nvars, &blas_l, &HH1, dG, &blas_l, &zero, H, &blas_nvars);


      // Jz'Jz from creating ddE and dD'*dD can be collapsed into the same matrix multiply
      for (size_t z = 0; z < Z; z++) {
         double iEsq = pow(EE[z], (r - 2.0) * 0.5);

         for (size_t j = 0; j < nvars; j++) Jz[j + z * nvars] *= iEsq;
      }

      HH1 = (r - 2.0) * E1r1;

      // E^(r-2)*Jz'*Jz == (E^(r/2-1)*Jz')*(E^(r/2-1)*Jz)
      // H += D^(r-2)*J'*J
      // H += D^(r-1)*(J'*J)/E
      dsyrk(uplo, transb, &blas_nvars, &blas_z, &HH1, Jz, &blas_nvars, &one, H, &blas_nvars);


      // H -= dR'*dR
      HH1 = -(r - 1.0) / E[0];
      dsyr(uplo, &blas_nvars, &HH1, J, &blas_one, H, &blas_nvars);



#pragma omp parallel for 
      for (size_t i = 0; i < N; i++) {
         H[(3 * i) + (3 * i) * nvars] += ddG2z[6 * i];
         H[(1 + 3 * i) + (3 * i) * nvars] += ddG2z[1 + 6 * i];
         H[(2 + 3 * i) + (3 * i) * nvars] += ddG2z[2 + 6 * i];
         H[(1 + 3 * i) + (1 + 3 * i) * nvars] += ddG2z[3 + 6 * i];
         H[(2 + 3 * i) + (1 + 3 * i) * nvars] += ddG2z[4 + 6 * i];
         H[(2 + 3 * i) + (2 + 3 * i) * nvars] += ddG2z[5 + 6 * i];
      }
   }


#pragma omp parallel default(none) firstprivate(nvars) shared(H)
   {
#pragma omp for //collapse(2)
      for (size_t i = 0; i < nvars; i++) {
         for (size_t j = i; j < nvars; j++) {
            H[i + j * nvars] = H[j + i * nvars];
         }
      }

   }

   free(mLinEq);
   free(nLinEq);
   free(Gz);
   free(Jz);
   free(ddGz);
   free(ddG2z);
   free(Er1);
   free(EE);
   free(mT);
   free(dG);
   free(ddG);
}

void calcEDIFFPHI(double* E, double* phi, double* alpha, double* mF, const size_t CS, const size_t SS, size_t const maxL, size_t const N, const size_t nvars, const size_t Z, double r) {
   size_t DF, iSymCombo, incr;
   size_t  L;
   double sum = 0.0, iA;
   //Symmetry Assignments
   size_t* mLinEq, * nLinEq;



   mLinEq = (size_t*)calloc(maxL + 1, sizeof(size_t));
   nLinEq = (size_t*)calloc(maxL + 1, sizeof(size_t));

   iSymCombo = getSymCombo(CS, SS);

   NumLinEq(iSymCombo, CS, maxL, mLinEq);
   NumLinEq(iSymCombo, SS, maxL, nLinEq);

   if (iSymCombo == 4 || iSymCombo == 3) {
      incr = 1;
   }
   else {
      incr = 2;
   }

   L = 0;
   for (size_t l = 0; l <= maxL; l++) {
      for (size_t mu = 1; mu <= mLinEq[l]; mu++) {
         for (size_t nu = 1; nu <= nLinEq[l]; nu++) {
            L += incr;
         }
      }
   }

   double* EE = (double*)calloc(Z, sizeof(double));
   double* mT = (double*)calloc(L, sizeof(double));

   {
      double A;
      kahansum(N, alpha[k], A)
         iA = 1.0 / A;
   }

   DF = 2;



   for (size_t j = 0; j < L; j++) {
      mT[j] = 0.0;
   }

   double* pmT;

#pragma omp parallel default(none) firstprivate(DF, iSymCombo, CS, SS, maxL, L, N, incr, mLinEq, nLinEq, iA) shared(phi, alpha, mT, pmT)
   {
      const size_t nthread = omp_get_num_threads();;
      const size_t ithread = omp_get_thread_num();
      size_t  i, j, l, mu, nu;

      double* sphi1 = (double*)calloc(maxL + 1, sizeof(double));
      double* cphi1 = (double*)calloc(maxL + 1, sizeof(double));
      double* sPhi = (double*)calloc(maxL + 1, sizeof(double));
      double* cPhi = (double*)calloc(maxL + 1, sizeof(double));
      double* sphi2 = (double*)calloc(maxL + 1, sizeof(double));
      double* cphi2 = (double*)calloc(maxL + 1, sizeof(double));
      double* mTc = (double*)malloc(L * sizeof(double));

      for (j = 0; j < L; j++) mTc[j] = 0.0;

#pragma omp single
      {
         pmT = (double*)malloc(L * nthread * sizeof(double));
      }

      for (i = 0; i < L; i++) pmT[i + L * ithread] = 0.0;

#pragma omp for// schedule(static)
      for (i = 0; i < N; i++) {

         double calpha;
         calpha = alpha[i] * iA;
         //pre-calculate all the sines and cosines we will need
         for (j = 0; j <= maxL; j++) {
            sphi1[j] = sin((double)j * phi[3 * i]);
            cphi1[j] = cos((double)j * phi[3 * i]);
            sPhi[j] = sin((double)j * phi[1 + 3 * i]);
            cPhi[j] = cos((double)j * phi[1 + 3 * i]);
            sphi2[j] = sin((double)j * phi[2 + 3 * i]);
            cphi2[j] = cos((double)j * phi[2 + 3 * i]);
         }

         // cnt counts the index of the fourier coeff
         size_t cnt;
         cnt = 0;
         for (l = 0; l <= maxL; l++) {
            // mu represents m for the crystal symmetry
            for (mu = 1; mu <= mLinEq[l]; mu++) {
               //nu represents n for the sample symmetry
               for (nu = 1; nu <= nLinEq[l]; nu++) {
                  //size_t l, mu, nu;
                  double kahany, kahant;
                  gshe_complex Tlmn;
                  // Tlmn is a structure that contains both the real and imaginary parts of the coeffs, as well as the derivatives
                  Tlmn = symmetricT(DF, iSymCombo, CS, SS, l, mu, nu, sphi1, cphi1, sPhi, cPhi, sphi2, cphi2);

                  kahany = (calpha * Tlmn.re) - mTc[cnt];
                  kahant = pmT[cnt + L * ithread] + kahany;

                  mTc[cnt] = (kahant - pmT[cnt + L * ithread]) - kahany;
                  pmT[cnt + L * ithread] = kahant;

                  cnt++;

                  if (incr > 1) {
                     kahany = (calpha * Tlmn.im) - mTc[cnt];
                     kahant = pmT[cnt + L * ithread] + kahany;

                     mTc[cnt] = (kahant - pmT[cnt + L * ithread]) - kahany;
                     pmT[cnt + L * ithread] = kahant;
                     cnt++;
                  }
               }
            }
         }
      }
#pragma omp for
      for (i = 0; i < L; i++) {
         double ssum;
         ssum = 0.0;
         kahansum(nthread, pmT[i + k * L], ssum)
            mT[i] = ssum;
      }

      free(sphi1); free(cphi1);
      free(sPhi); free(cPhi);
      free(sphi2); free(cphi2);
      free(mTc);
   }
   free(pmT);

#ifdef RDEBUG
   mexPrintf("mT = [");
   for (int jj = 0; jj < L; jj++) mexPrintf(" %.4f", mT[jj]);
   mexPrintf("];\n");
#endif

   //#pragma omp parallel default(none) firstprivate(L, N, nvars, Z, iA) shared(EE, J, mF, mT, dG, ddG, pmJ)
   {
      const size_t nthread = omp_get_num_threads();
      const size_t ithread = omp_get_thread_num();
      double* Gl = (double*)calloc(L, sizeof(double));
      size_t l;

      //#pragma omp for// schedule(static)
      for (size_t z = 0; z < Z; z++) {


         {
            //double kahanc, kahany, kahant;
            //kahanc = 0.0;
            double ssum = 0.0;
            for (l = 0; l < L; l++) {
               Gl[l] = mT[l] - mF[l + z * L];

               //kahany = Gl[l] * Gl[l] - kahanc;
               //kahant = ssum + kahany;
               //kahanc = (kahant - ssum) - kahany;
               //ssum = kahant;
            }

            kahansum(L, Gl[k] * Gl[k], ssum);

            EE[z] = sqrt(ssum);
         }
      }



      free(Gl);
   }

   {
      double ssum = 0.0;
      kahansum(Z, pow(EE[k], r), ssum);
      E[0] = pow(ssum, (1.0 / r));
   }






   //mexPrintf("exit mf\n");


   free(mLinEq);
   free(nLinEq);
   free(EE);
   free(mT);
}

void calcMAXPHI(double* E, double* J, double* H, double* phi, double* alpha, double* F, double* kernel, const size_t CS, const size_t SS, size_t const maxL, size_t const N, const size_t nvars, const size_t Z, double r) {
   size_t DF, iSymCombo, incr;
   size_t  L;
   double sum = 0.0, iA;
   //Symmetry Assignments
   size_t* mLinEq, * nLinEq;

   mLinEq = (size_t*)calloc(maxL + 1, sizeof(size_t));
   nLinEq = (size_t*)calloc(maxL + 1, sizeof(size_t));

   iSymCombo = getSymCombo(CS, SS);

   NumLinEq(iSymCombo, CS, maxL, mLinEq);
   NumLinEq(iSymCombo, SS, maxL, nLinEq);

   double K = 6.0;

   if (iSymCombo == 4 || iSymCombo == 3) {
      incr = 1;
   }
   else {
      incr = 2;
   }

   L = 0;
   for (size_t l = 0; l <= maxL; l++) {
      for (size_t mu = 1; mu <= mLinEq[l]; mu++) {
         for (size_t nu = 1; nu <= nLinEq[l]; nu++) {
            L += incr;
         }
      }
   }
   double* mT = (double*)calloc(L, sizeof(double));
   double* dG = (double*)calloc(3 * N * L, sizeof(double));
   double* ddG = (double*)calloc(6 * N * L, sizeof(double));

   {
      double A;
      kahansum(N, alpha[k], A)
         iA = 1.0 / A;
   }

   DF = 2;



   for (size_t j = 0; j < L; j++) {
      mT[j] = 0.0;
   }

   double* pmT;

#pragma omp parallel default(none) firstprivate(DF, iSymCombo, CS, SS, maxL, L, N, incr, mLinEq, nLinEq, iA) shared(phi, alpha, kernel, dG, ddG, mT, pmT)
   {
      const size_t nthread = omp_get_num_threads();
      const size_t ithread = omp_get_thread_num();
      size_t  i, j, l, mu, nu;

      double* sphi1 = (double*)calloc(maxL + 1, sizeof(double));
      double* cphi1 = (double*)calloc(maxL + 1, sizeof(double));
      double* sPhi = (double*)calloc(maxL + 1, sizeof(double));
      double* cPhi = (double*)calloc(maxL + 1, sizeof(double));
      double* sphi2 = (double*)calloc(maxL + 1, sizeof(double));
      double* cphi2 = (double*)calloc(maxL + 1, sizeof(double));
      double* mTc = (double*)malloc(L * sizeof(double));

      for (j = 0; j < L; j++) mTc[j] = 0.0;

#pragma omp single
      {
         pmT = (double*)malloc(L * nthread * sizeof(double));
      }

      for (i = 0; i < L; i++) pmT[i + L * ithread] = 0.0;

#pragma omp for// schedule(static)
      for (i = 0; i < N; i++) {

         double calpha;
         calpha = alpha[i] * iA;
         //pre-calculate all the sines and cosines we will need
         for (j = 0; j <= maxL; j++) {
            sphi1[j] = sin((double)j * phi[3 * i]);
            cphi1[j] = cos((double)j * phi[3 * i]);
            sPhi[j] = sin((double)j * phi[1 + 3 * i]);
            cPhi[j] = cos((double)j * phi[1 + 3 * i]);
            sphi2[j] = sin((double)j * phi[2 + 3 * i]);
            cphi2[j] = cos((double)j * phi[2 + 3 * i]);
         }

         // cnt counts the index of the fourier coeff
         size_t cnt;
         cnt = 0;
         for (l = 0; l <= maxL; l++) {
            double pd_fac = kernel[l];
            // mu represents m for the crystal symmetry
            for (mu = 1; mu <= mLinEq[l]; mu++) {
               //nu represents n for the sample symmetry
               for (nu = 1; nu <= nLinEq[l]; nu++) {
                  //size_t l, mu, nu;
                  double kahany, kahant;
                  gshe_complex Tlmn;
                  // Tlmn is a structure that contains both the real and imaginary parts of the coeffs, as well as the derivatives
                  Tlmn = symmetricT(DF, iSymCombo, CS, SS, l, mu, nu, sphi1, cphi1, sPhi, cPhi, sphi2, cphi2);
                  Tlmn = sc_mult(DF, Tlmn, pd_fac);

                  kahany = (calpha * Tlmn.re) - mTc[cnt];
                  kahant = pmT[cnt + L * ithread] + kahany;

                  mTc[cnt] = (kahant - pmT[cnt + L * ithread]) - kahany;
                  pmT[cnt + L * ithread] = kahant;

                  //dG
                  dG[cnt + 3 * i * L] = Tlmn.dre[0] * calpha;
                  dG[1 + cnt + 3 * i * L] = Tlmn.dim[0] * calpha;
                  dG[cnt + (3 * i + 1) * L] = Tlmn.dre[1] * calpha;//
                  dG[1 + cnt + (3 * i + 1) * L] = Tlmn.dim[1] * calpha;
                  dG[cnt + (3 * i + 2) * L] = Tlmn.dre[2] * calpha;//
                  dG[1 + cnt + (3 * i + 2) * L] = Tlmn.dim[2] * calpha;

                  //ddG
                  ddG[cnt + 6 * i * L] = Tlmn.ddre[0] * calpha;
                  ddG[1 + cnt + 6 * i * L] = Tlmn.ddim[0] * calpha;
                  ddG[cnt + (6 * i + 1) * L] = Tlmn.ddre[1] * calpha;//
                  ddG[1 + cnt + (6 * i + 1) * L] = Tlmn.ddim[1] * calpha;
                  ddG[cnt + (6 * i + 2) * L] = Tlmn.ddre[2] * calpha;//
                  ddG[1 + cnt + (6 * i + 2) * L] = Tlmn.ddim[2] * calpha;
                  ddG[cnt + (6 * i + 3) * L] = Tlmn.ddre[3] * calpha;
                  ddG[1 + cnt + (6 * i + 3) * L] = Tlmn.ddim[3] * calpha;
                  ddG[cnt + (6 * i + 4) * L] = Tlmn.ddre[4] * calpha;//
                  ddG[1 + cnt + (6 * i + 4) * L] = Tlmn.ddim[4] * calpha;
                  ddG[cnt + (6 * i + 5) * L] = Tlmn.ddre[5] * calpha;//
                  ddG[1 + cnt + (6 * i + 5) * L] = Tlmn.ddim[5] * calpha;


                  cnt++;

                  if (incr > 1) {
                     kahany = (calpha * Tlmn.im) - mTc[cnt];
                     kahant = pmT[cnt + L * ithread] + kahany;

                     mTc[cnt] = (kahant - pmT[cnt + L * ithread]) - kahany;
                     pmT[cnt + L * ithread] = kahant;
                     cnt++;
                  }
               }
            }
         }
      }
#pragma omp for
      for (i = 0; i < L; i++) {
         double ssum;
         ssum = 0.0;
         kahansum(nthread, pmT[i + k * L], ssum);
         mT[i] = ssum;
      }

      free(sphi1); free(cphi1);
      free(sPhi); free(cPhi);
      free(sphi2); free(cphi2);
      free(mTc);
   }
   free(pmT);

   double* p = (double*)malloc(Z * sizeof(double));
   double* dP = (double*)calloc(3 * N * Z, sizeof(double));
   double* ddP = (double*)calloc(6 * N * Z, sizeof(double));
   double* Er1 = (double*)malloc(Z * sizeof(double));

   double* dR = (double*)calloc(3 * N, sizeof(double));
   double* ddR = (double*)calloc(6 * N, sizeof(double));

   {

      char* transa;
      char* transb;
      char* uplo;
      ptrdiff_t blas_nvars, blas_l, blas_z, blas_one, blas_dvars;
      double zero, one, R, HH2, Er, E1r1, HH1;
      one = 1.0;
      zero = 0.0;
      blas_nvars = nvars;
      blas_dvars = 6 * N;
      blas_l = L;
      blas_z = Z;
      blas_one = 1;
      transa = "t";
      transb = "n";
      uplo = "L";

      dgemv(transa, &blas_l, &blas_z, &one, F, &blas_l, mT, &blas_one, &zero, p, &blas_one);

      kahansum(Z, pow(p[k], r), Er);
      E1r1 = pow(Er, (1.0 / r - 1.0));
      R = pow(Er, 1.0 / r);

      E[0] = (R - K) * (R - K);

#pragma omp parallel for
      for (size_t z = 0; z < Z; z++) {
         Er1[z] = pow(p[z], (r - 1.0));
      }

      dgemm(transa, transb, &blas_nvars, &blas_z, &blas_l, &one, dG, &blas_l, F, &blas_l, &zero, dP, &blas_nvars);

      dgemv(transb, &blas_nvars, &blas_z, &E1r1, dP, &blas_nvars, Er1, &blas_one, &zero, dR, &blas_one);

      dgemm(transa, transb, &blas_dvars, &blas_z, &blas_l, &one, ddG, &blas_l, F, &blas_l, &zero, ddP, &blas_dvars);

      dgemv(transb, &blas_dvars, &blas_z, &E1r1, ddP, &blas_dvars, Er1, &blas_one, &zero, ddR, &blas_one);

#pragma omp parallel for
      for (size_t i = 0; i < nvars; i++) {
         J[i] = 2.0 * (R - K) * dR[i];
      }

#pragma omp parallel for
      for (size_t z = 0; z < Z; z++) {
         double iE = pow(p[z], (r - 2.0) * 0.5);
         for (size_t i = 0; i < nvars; i++) {
            dP[i + z * nvars] *= iE;
         }
      }

      HH1 = (r - 1.0) * E1r1;

      dsyrk(uplo, transb, &blas_nvars, &blas_z, &HH1, dP, &blas_nvars, &zero, H, &blas_nvars);

      HH1 = 1.0 / (R - K) - (r - 1.0) / R;

      dsyr(uplo, &blas_nvars, &HH1, dR, &blas_one, H, &blas_nvars);


#pragma omp parallel for 
      for (size_t i = 0; i < N; i++) {
         H[(3 * i) + (3 * i) * nvars] += ddR[6 * i];
         H[(1 + 3 * i) + (3 * i) * nvars] += ddR[1 + 6 * i];
         H[(2 + 3 * i) + (3 * i) * nvars] += ddR[2 + 6 * i];
         H[(1 + 3 * i) + (1 + 3 * i) * nvars] += ddR[3 + 6 * i];
         H[(2 + 3 * i) + (1 + 3 * i) * nvars] += ddR[4 + 6 * i];
         H[(2 + 3 * i) + (2 + 3 * i) * nvars] += ddR[5 + 6 * i];
      }

#pragma omp parallel for 
      for (size_t i = 0; i < nvars; i++) {
         for (size_t j = i; j < nvars; j++) {
            H[j + i * nvars] *= 2.0 * (R - K);
            H[i + j * nvars] = H[j + i * nvars];
         }
      }



   }




   //mexPrintf("exit mf\n");


   free(mLinEq);
   free(nLinEq);
   free(mT);
   free(p);
   free(dP);
   free(ddP);
   free(dG);
   free(ddG);
   free(dR);
}

void calcEMAXPHI(double* E, double* phi, double* alpha, double* F, double* kernel, const size_t CS, const size_t SS, size_t const maxL, size_t const N, const size_t nvars, const size_t Z, double r) {
   size_t DF, iSymCombo, incr;
   size_t  L;
   double sum = 0.0, iA;
   //Symmetry Assignments
   size_t* mLinEq, * nLinEq;

   mLinEq = (size_t*)calloc(maxL + 1, sizeof(size_t));
   nLinEq = (size_t*)calloc(maxL + 1, sizeof(size_t));

   iSymCombo = getSymCombo(CS, SS);

   NumLinEq(iSymCombo, CS, maxL, mLinEq);
   NumLinEq(iSymCombo, SS, maxL, nLinEq);

   double K = 6;

   if (iSymCombo == 4 || iSymCombo == 3) {
      incr = 1;
   }
   else {
      incr = 2;
   }

   L = 0;
   for (size_t l = 0; l <= maxL; l++) {
      for (size_t mu = 1; mu <= mLinEq[l]; mu++) {
         for (size_t nu = 1; nu <= nLinEq[l]; nu++) {
            L += incr;
         }
      }
   }
   double* mT = (double*)calloc(L, sizeof(double));


   {
      double A;
      kahansum(N, alpha[k], A)
         iA = 1.0 / A;
   }

   DF = 2;



   for (size_t j = 0; j < L; j++) {
      mT[j] = 0.0;
   }

   double* pmT;

#pragma omp parallel default(none) firstprivate(DF, iSymCombo, CS, SS, maxL, L, N, incr, mLinEq, nLinEq, iA) shared(phi, alpha, kernel, mT, pmT)
   {
      const size_t nthread = omp_get_num_threads();
      const size_t ithread = omp_get_thread_num();
      size_t  i, j, l, mu, nu;

      double* sphi1 = (double*)calloc(maxL + 1, sizeof(double));
      double* cphi1 = (double*)calloc(maxL + 1, sizeof(double));
      double* sPhi = (double*)calloc(maxL + 1, sizeof(double));
      double* cPhi = (double*)calloc(maxL + 1, sizeof(double));
      double* sphi2 = (double*)calloc(maxL + 1, sizeof(double));
      double* cphi2 = (double*)calloc(maxL + 1, sizeof(double));
      double* mTc = (double*)malloc(L * sizeof(double));

      for (j = 0; j < L; j++) mTc[j] = 0.0;

#pragma omp single
      {
         pmT = (double*)malloc(L * nthread * sizeof(double));
      }

      for (i = 0; i < L; i++) pmT[i + L * ithread] = 0.0;

#pragma omp for// schedule(static)
      for (i = 0; i < N; i++) {

         double calpha;
         calpha = alpha[i] * iA;
         //pre-calculate all the sines and cosines we will need
         for (j = 0; j <= maxL; j++) {
            sphi1[j] = sin((double)j * phi[3 * i]);
            cphi1[j] = cos((double)j * phi[3 * i]);
            sPhi[j] = sin((double)j * phi[1 + 3 * i]);
            cPhi[j] = cos((double)j * phi[1 + 3 * i]);
            sphi2[j] = sin((double)j * phi[2 + 3 * i]);
            cphi2[j] = cos((double)j * phi[2 + 3 * i]);
         }

         // cnt counts the index of the fourier coeff
         size_t cnt;
         cnt = 0;
         for (l = 0; l <= maxL; l++) {
            double pd_fac = kernel[l];
            // mu represents m for the crystal symmetry
            for (mu = 1; mu <= mLinEq[l]; mu++) {
               //nu represents n for the sample symmetry
               for (nu = 1; nu <= nLinEq[l]; nu++) {
                  //size_t l, mu, nu;
                  double kahany, kahant;
                  gshe_complex Tlmn;
                  // Tlmn is a structure that contains both the real and imaginary parts of the coeffs, as well as the derivatives
                  Tlmn = symmetricT(DF, iSymCombo, CS, SS, l, mu, nu, sphi1, cphi1, sPhi, cPhi, sphi2, cphi2);
                  Tlmn = sc_mult(DF, Tlmn, pd_fac);
                  kahany = (calpha * Tlmn.re) - mTc[cnt];
                  kahant = pmT[cnt + L * ithread] + kahany;

                  mTc[cnt] = (kahant - pmT[cnt + L * ithread]) - kahany;
                  pmT[cnt + L * ithread] = kahant;

                  cnt++;

                  if (incr > 1) {
                     kahany = (calpha * Tlmn.im) - mTc[cnt];
                     kahant = pmT[cnt + L * ithread] + kahany;

                     mTc[cnt] = (kahant - pmT[cnt + L * ithread]) - kahany;
                     pmT[cnt + L * ithread] = kahant;
                     cnt++;
                  }
               }
            }
         }
      }
#pragma omp for
      for (i = 0; i < L; i++) {
         double ssum;
         ssum = 0.0;
         kahansum(nthread, pmT[i + k * L], ssum);
         mT[i] = ssum;
      }

      free(sphi1); free(cphi1);
      free(sPhi); free(cPhi);
      free(sphi2); free(cphi2);
      free(mTc);
   }
   free(pmT);

   double* p = (double*)malloc(L * sizeof(double));

   {

      char* transa;
      char* transb;
      char* uplo;
      ptrdiff_t blas_nvars, blas_l, blas_z, blas_one, blas_dvars;
      double zero, one, R, HH2, Er, HH1;
      one = 1.0;
      zero = 0.0;
      blas_nvars = nvars;
      blas_dvars = 6 * N;
      blas_l = L;
      blas_z = Z;
      blas_one = 1;
      transa = "t";
      transb = "n";
      uplo = "L";

      dgemv(transa, &blas_l, &blas_z, &one, F, &blas_l, mT, &blas_one, &zero, p, &blas_one);

      kahansum(Z, pow(p[k], r), Er);
      R = pow(Er, 1.0 / r);

      E[0] = (R - K) * (R - K);
   }




   //mexPrintf("exit mf\n");


   free(mLinEq);
   free(nLinEq);
   free(mT);
   free(p);
}


