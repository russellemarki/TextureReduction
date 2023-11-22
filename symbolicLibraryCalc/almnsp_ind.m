function[out] = almnsp_ind(l, m, n, s)
%m:l
%m - l/6 + n + s + l*m*4 + 2*l*n + l^2*m*4 - l^2/2 + (2*l^3)/3 + l^4
l = l/2;
m = m/2;
out = (6*l*l - 2*l - 1)*(l+1)*l/6 + (2*l+1)*((2*l+1)*m + n) + s + 1;