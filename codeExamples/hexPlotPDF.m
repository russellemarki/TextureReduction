clearvars
close all
clc

addpath('../')

%% Inputs

% filenames, texture should be in VPSC format
fname_1 = '../Tex/Ti_625_0.TEX';
fname_2 = '../Tex/Ti_Wu_0.TEX';

name_1 = '625 C';
name_2 = 'Wu';

CS = crystalSymmetry('hexagonal'); %Crystal Symmetry
SS = specimenSymmetry('triclinic'); %Sample Symmetry

% bandwidth of fourier coefficients, higher bandwidth will require more crystals
% calcDT has a maximum bandwidth of 34
maxL = 16;

% for plotting purposes
psi = deLaValleePoussinKernel('halfwidth', 1*degree, 'bandwidth', maxL);

%% 

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

%% test

figure('Units', 'Normalized', 'Position', [0.1850 0.2000 0.6000 0.6157])
subplot(2,1,1)
%rplotPDF can be used with subplots
cbound = rplotPDF(odf_1, 'textVert', name_1, 'axisLabel', {'TD', 'RD'});
subplot(2,1,2)
rplotPDF(odf_2, 'textVert', name_2, 'colorBound', cbound);%

exportgraphics(gcf,'test.png')

%% test_sidebyside

c = 0.03;
figure('Units', 'Normalized', 'Position', [0.0505 0.5481 1.0495 0.3046])

subplot('Position',[(0.0+c) 0.0 (0.5-c) 1.0])
cbound1 = rplotPDF(odf_1, 'textVert', name_1, 'iColorBar', false);%

subplot('Position',[(0.5) 0.0 (0.5) 1.0])
rplotPDF(odf_2, 'textVert', name_2,'colorBound',cbound1);

exportgraphics(gcf,'test_sidebyside.png')