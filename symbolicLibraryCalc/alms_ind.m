function[out] = alms_ind(l, m, n)
l = l/2;
out = (2*l-1)*(2*l+1)*l/3 + m*(2*l+1) + n + 1;