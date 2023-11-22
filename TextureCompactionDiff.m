clearvars
close all
clc

set(0, 'defaultAxesTickLabelInterpreter','latex'); set(0, 'defaultLegendInterpreter','latex'); set(0, 'defaultTextInterpreter','latex'); set(0,'DefaultLineLineWidth',1.5); set(0,'DefaultAxesFontSize',14); set(0,'DefaultLegendFontSize',12,'DefaultLegendFontSizeMode','manual')
 
figure('Units', 'Normalized', 'Position', [0.3448 0.2037 1.0672 0.6157])

NN = [50 100 200 400 800 1600];%

for ii = 1:length(NN)
    maxL = 16;
    miscinputs = {'MaxIter', 250};

    CS = crystalSymmetry('orthorhombic'); %Crystal Symmetry
    SS = specimenSymmetry('triclinic'); %Sample Symmetry

    N = NN(ii);

    kk = 200;

    phi_grid = grabPhiGrid(CS);

    psi_scl = 1.0./(1:maxL + 1);
    psi_scl(9:end) = 0;
    %psi_scl(maxL_scl + 1) = 1;
    psi_inv = 0* psi_scl;
    psi_inv(9:end) = 1;

    mF = calcGSH(1, ones(3, 1), 1, CS.id, SS.id, maxL);
    mF(2:end) = 0.0;

    %Initialize Czl
    pF = calcGSH(11, ones(3, 1), 1, CS.id, SS.id, maxL, psi_scl);
    pF = [pF*0.0, zeros(size(pF, 1), kk)];

    dpow = 2;

    for i = 1:kk

        fname_out = sprintf('Tex/OrthTric/diff_%iL_%i_%03i.txt',maxL, N, i);
        pname_out = sprintf('Tex/OrthTric/diff_%iL_%i_%03i.png',maxL, N, i);

        indr = randperm(size(phi_grid, 2));
        iGuess = [phi_grid(:, indr(1:N)); ones(1, N)./N];

        tic
        % Run the compaction code
        [phis, alpha] = fourier2ori_max(mF, pF(:, 1:i), CS.id, SS.id, maxL, psi_scl, psi_inv, dpow, iGuess);
        toc

        pF_new = calcGSH(11, phis, alpha, CS.id, SS.id, maxL, psi_scl);

        %add Cl to Czl
        pF(:, i+1) = pF_new;

        %write the compacted texture to a VPSC formatted text file
        writeVPSC(phis, alpha, fname_out)

        if false
            subplot('Position', [0.02 0 0.7 1])
            
            psi = DirichletKernel(maxL);

            ori = calcOri(phis, CS, SS);

            %load odf of input
            odf = calcDensity(ori, 'weights', alpha, 'kernel', psi);
            %rplotPDF can be used with subplots
            cbound = rplotPDF(odf, 'textVert', sprintf('%iLs; %iL; %iN', maxL, maxL, N), 'axisLabel', {'RD', 'TD'});%, 'textVert', sprintf('fit=%i plot=%i N=%i', maxL, maxL_plot, N));

            lt = calcGSH(12, phis, alpha, CS.id, SS.id, maxL);

            subplot('Position', [0.75 0.1 0.22 0.8])

            semilogy(1:maxL, lt(2:end), '.','MarkerSize',60)
            xlabel('max L')
            title('L $$\left|\Sigma\right|$$')

            drawnow
%             exportgraphics(gcf,pname_out)
        end
    end
end