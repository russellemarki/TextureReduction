
clear all
mex COMPFLAGS="$COMPFLAGS /openmp" calcGSH.cpp GSHFun.cpp GSHLibrarySmol.cpp -lmwblas

%% 
maxL = 4;
CS = 16;
SS = 2;

mF = calcGSH(1, rand(3, 10), rand(1, 10), CS, SS, maxL);

% err = @(x) (calcGSH(5, x, mF, CS, SS, maxL)).^2;

xo = rand(3, 10);

alpha = rand(1, 10);

xo = xo(:);

err = @(x) errfun(x, alpha, mF, CS, SS, maxL);

J = gradest(err, xo)';

H = hessian(err, xo);

[E0, J0, H0] = calcGSH(14, xo, alpha, mF, CS, SS, maxL);

test = H0./H;



function[out] = errfun(x, alpha, mF, CS, SS, maxL)
x = reshape(x, [3, length(alpha)]);
x = [x; alpha(:)'];

out = (calcGSH(5, x, mF, CS, SS, maxL)).^2;
end