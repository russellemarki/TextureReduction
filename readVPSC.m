function[varargout] = readVPSC(fname)

[euler] = importdata(fname, ' ', 4);
if iscell(euler)
    [euler] = importdata(fname, '\t', 4);
end

if isstruct(euler)
    euler = euler.data;
end

phis = euler(:, 1:3)';

if any(phis(:)>30)
    phis = phis*pi/180;
end
wgt = euler(:, 4);

swgt = kahansum(wgt);

if abs(swgt - 1) < 1e-5
    wgt = wgt./swgt;
end

if nargout == 2
  varargout = {phis, wgt};
  elseif nargout == 1
  phis = [phis; wgt(:)'];
  varargout = {phis};
end