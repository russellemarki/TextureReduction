clearvars
close all
clc

addpath('../')

%% Inputs

% filenames, texture should be in VPSC format
fname_in = '../Tex/DP590_initial.txt';
fname_out = 'test.txt';

CS = crystalSymmetry('cubic'); %Crystal Symmetry
SS = specimenSymmetry('triclinic'); %Sample Symmetry

% bandwidth of fourier coefficients, higher bandwidth will require more crystals
maxL = 16;

% number of crystals to be output
N = 95;

% maxL for plotting purposes
maxL_plot = maxL;

miscinputs = {'fixAlpha', true};

%% Running the Code

[phis_target, alpha_target] = readVPSC(fname_in);

% create orientations in mtex format
ori_target = calcOri(phis_target, CS, SS);

mF_target = calcGSH(1, phis_target, alpha_target, CS.id, SS.id, maxL);

% Run the compaction code
tic
[ori, alpha] = fourier2ori(mF_target, CS, SS, maxL, N, miscinputs{:});
toc

phis = calcPhi(ori);

%write the compacted texture to a VPSC formatted text file
writeVPSC(phis, alpha, fname_out)

%% Plotting PDFs

%for plotting purposes
psi = DirichletKernel(maxL_plot);

%load odf of input
odf_in = calcDensity(ori_target, 'weights', alpha_target, 'kernel', psi);

%load odf of output
odf_out = calcDensity(ori, 'weights', alpha, 'kernel', psi);

figure('Units', 'Normalized', 'Position', [0.1850 0.2000 0.6000 0.6157])
subplot(2,1,1)
%rplotPDF can be used with subplots
cbound = rplotPDF(odf_in, 'textVert', 'Input Tex', 'axisLabel', {'TD', 'RD'});
subplot(2,1,2)
rplotPDF(odf_out, 'textVert', 'Compacted Tex', 'colorBound', cbound);%