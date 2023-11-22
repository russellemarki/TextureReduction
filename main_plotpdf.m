clearvars
close all
clc

%% Inputs

% filenames, texture should be in VPSC format
fname_1 = 'Tex/OrthTric/rand_24L_kernel_001.txt';
fname_2 = 'Tex/Ti_Wu_0.TEX';

CS = crystalSymmetry('orthorhombic'); %Crystal Symmetry
SS = specimenSymmetry('triclinic'); %Sample Symmetry

% bandwidth of fourier coefficients, higher bandwidth will require more crystals
% calcDT has a maximum bandwidth of 34
maxL = 16;

% for plotting purposes
psi = deLaValleePoussinKernel('halfwidth', 1*degree, 'bandwidth', maxL);


[phis, alpha] = readVPSC(fname_1);

% create orientations in mtex format
ori = calcOri(phis, CS, SS);

%load odf of input
odf_1 = calcDensity(ori, 'weights', alpha, 'kernel', psi);


[phis, alpha] = readVPSC(fname_2);

% create orientations in mtex format
ori = calcOri(phis, CS, SS);

%load odf of output
odf_2 = calcDensity(ori, 'weights', alpha, 'kernel', psi);

figure('Units', 'Normalized', 'Position', [0.1850 0.2000 0.6000 0.6157])

subplot(2,1,1)
%rplotPDF can be used with subplots
cbound = rplotPDF(odf_1, 'textVert', 'SR 10000', 'axisLabel', {'TD', 'RD'});
subplot(2,1,2)
rplotPDF(odf_2, 'textVert', 'Compacted Tex', 'colorBound', cbound);%


exportgraphics(gcf,'odf_plotted.png')