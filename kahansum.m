function[sum] = kahansum(x, dir)
% Kahan summation algorithm in the 2nd direction if x
% is more than 1-D

if nargin < 2
    dir = 0;
end

[sz1, sz2] = size(x);

 if ((sz1 > 1) && (sz2 > 1))
    sum = zeros(sz1, 1);
    c = zeros(sz1, 1);
    
elseif (sz1 > 1 && sz2 == 1) && dir ~= 1
    
    x = x(:)';
    sum = 0;
    c = 0;
    sz2 = sz1;
    
else
    sum = 0;
    c = 0;
    
end

for i = 1:sz2
    y = x(:, i) - c;
    t = sum + y;
    
    c = (t - sum) - y;
    
    sum = t;
end
end