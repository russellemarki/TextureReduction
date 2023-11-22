function[ori] = calcOri(phi, CS, SS)
%wrap mtex format in case I forget to add transpose
ori = orientation('Euler', phi', CS, SS);
end