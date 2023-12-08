
clear all
mex COMPFLAGS="$COMPFLAGS /openmp" calcGSH.cpp GSHFun.cpp GSHLibrarySmol.cpp -lmwblas

%% 
maxL = 4;
CS = 16;
SS = 2;

mF = calcGSH(1, rand(3, 10), rand(1, 10), CS, SS, maxL);
