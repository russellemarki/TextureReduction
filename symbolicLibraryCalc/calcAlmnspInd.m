clearvars
close all
clc


maxL = 4;
pow = 5;

% sym = @(z) z;

cnt = 0;
for l = 0:maxL
    for m = 0:2:2*l
        for n = 0:2*l
            for s = 0:2*l
                cnt = cnt + 1;
            end
        end
    end
end

A = sym(zeros(217, cnt));
b = sym(zeros(cnt, 1));

cnt = 0;

for l = 0:maxL
    for m = 0:2:2*l
        for n = 0:2*l
            for s = 0:2*l
                cnt = cnt + 1;

                ncnt = 0;
                for i = 0:pow
                    for j = 0:pow
                        for k = 0:pow
                            ncnt = ncnt + 1;
                            A(ncnt, cnt) = sym(l)^i*sym(m)^j*sym(n)^k;
                        end
                    end
                end
                A(ncnt+1, cnt) = sym(s);

                b(cnt) = sym(cnt-1);
            end
        end
    end
end

%% 

%m:2l
%m + n + s + 4*l*m + 2*l*n + 4*l^2*m - l^2 + 2*l^4

%m:2:2l
%m/2 - l/6 + n + s + 2*l*m + 2*l*n + 2*l^2*m - l^2/2 + (2*l^3)/3 + l^4

%m:l
%m - l/6 + n + s + l*m*4 + 2*l*n + l^2*m*4 - l^2/2 + (2*l^3)/3 + l^4

x = (transpose(A))\b;
% 
%  out = abs((A')*dx - b);
% 
%  figure
%  histogram(out)

% x(abs(x) < 1e-4) = 0;
% 
% x(2) = 49/25;
% x(21) = 73/400;
% x(26) = -1/12;
% x(51) = 1/8;
% x(52) = 7/20;
% x(66) = 1/2;

%% 

l = sym('l', 'real');
m = sym('m', 'real');
n = sym('n', 'real');
s = sym('s', 'real');

AA = sym(zeros(ncnt, 1));

ncnt = 0;
for i = 0:pow
    for j = 0:pow
        for k = 0:pow
            ncnt = ncnt + 1;
            AA(ncnt, 1) = l^i*m^j*n^k;
        end
    end
end
AA(ncnt+1, 1) = s;
%% 
cnt = 0;

test = zeros(20, 1);

for l = 0:2:maxL
    for m = 0:2:l
        for n = 0:l
            for s = 0:l
                cnt = cnt + 1;
                test(cnt) = almnsp_ind(l, m, n, s);
            end
        end
    end
end
 all(diff(test)==1)