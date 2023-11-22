clearvars
close all
clc


%% Inputs

CS = crystalSymmetry('cubic'); %Crystal Symmetry
SS = specimenSymmetry('triclinic'); %Sample Symmetry

% bandwidth of fourier coefficients, higher bandwidth worse for fewer
% crystals
% calcDT has a maximum bandwidth of 64
maxL = 32;

psik = 1.0./(1.0*(0:maxL) + 1.0);


%% Running the Code

mF_target = calcGSH(1, phis_target, alpha_target, CS.id, SS.id, maxL, psik);

% create orientations in mtex format
ori_target = calcOri(phis_target, CS, SS);

NN = [50 100 200 400 800 1600];%

for j = 1:length(NN)
    % calculate the target mean Fourier Coeffs
    mF = calcGSH(1, [1; 1; 1], 1, CS.id, SS.id, maxL);

    mF(2:end) = 0.0;
    mF(1) = 1.0;

    N = NN(j);
    for i = 1:200
        miscinputs = {'MaxIter', 2000, 'kernel', psik};

        fname_out = sprintf('Tex/CubTric/unidif1l1_%iL_%iN_%03i.txt',maxL, N, i);
        fname_in = sprintf('Tex/OrthTric/diff_16L_%i_%03i.txt', N, i);

        % Run the compaction code
        tic
        [ori, alpha] = fourier2ori(mF, CS, SS, maxL, iGuess, miscinputs{:});
        toc

        phis = calcPhi(ori);

        %write the compacted texture to a VPSC formatted text file
        writeVPSC(phis, alpha, fname_out)
    end
end
