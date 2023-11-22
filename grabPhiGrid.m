function[phi_grid] = grabPhiGrid(CS)
if CS.id == 16
    [phi_grid, ~] = readVPSC('Tex/grid15k_orthtric.txt');
elseif CS.id == 45
    [phi_grid, ~] = readVPSC('Tex/random_cubtric.txt');
elseif CS.id == 40
    [phi_grid, ~] = readVPSC('Tex/random_hextric.txt');
else
    [phi_grid, ~] = readVPSC('Tex/random_allcs.txt');
end
end