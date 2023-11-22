function[phi] = calcPhi(ori)
%wrap mtex format in case I forget to add transpose
phi = ori.Euler();
phi = phi';
end