function[] = writeVPSC(euler, wgt, fname)
  
if nargin == 2
  fname = wgt;
  wgt = euler(4, :);
  euler = euler(1:3, :);
end
  
n = size(euler, 2);

if all(wgt==wgt(1))
    wgt = ones(1, n);
else
    wgt = n*wgt./sum(wgt);
end

fid = fopen(fname, 'w');

fprintf(fid, 'TEXTURE AT STRAIN\n\n\nB  %.0f\n', n);

for i = 1:n
    fprintf(fid, '%24.15g', euler(:, i)*180.0/pi);
    fprintf(fid, '%24.15g\n', wgt(i));
end

fclose(fid);
end