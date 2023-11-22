clearvars
close all
clc


maxL = 4;
pow = 5;

% sym = @(z) z;

cnt = 0;
for l = 0:maxL
    for m = 0:2*l
        for n = 0:2*l
            cnt = cnt + 1;
        end
    end
end

A = sym(zeros(20, cnt));
b = sym(zeros(cnt, 1));

cnt = 0;

for l = 0:maxL
    for m = 0:2*l
        for n = 0:2*l
            cnt = cnt + 1;

            ncnt = 0;
            for i = 0:pow
                for j = 0:pow
                    ncnt = ncnt + 1;
                    A(ncnt, cnt) = sym(l)^i*sym(m)^j;
                end
            end
            A(ncnt+1, cnt) = sym(n);

            b(cnt) = sym(cnt-1);
        end
    end
end

%%

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

AA = sym(zeros(1, ncnt));

ncnt = 0;
for i = 0:pow
    for j = 0:pow
        ncnt = ncnt + 1;
        AA(ncnt) = l^i*m^j;
    end
end
AA(ncnt+1) = n;

AA*x
%%
cnt = 0;

test = zeros(20, 1);

for l = 0:2:maxL
    for m = 0:l
        for n = 0:l
            cnt = cnt + 1;
            test(cnt) = alms_ind(l, m, n);
        end
    end
end
all(diff(test)==1)