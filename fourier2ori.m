function[ori, alpha] = fourier2ori(mF, CS, SS, maxL, N, varargin)
%FOURIER2ORI will convert the input mF Fourier Coeffs to a set of
%orientations
%
%    mF is the mean Fourier coefficients
%
%    calcRF is a function handle that calculates fourier coefficients from
%    orientations
%
%    if N is a scalar it is the number of crystals used in the compaction,
%    if N is a vector it is the initial guess to be used for the
%    compaction.
% -------------------------------------------------------------------------
%    author:      Russell Marki
%    affiliation: Marko Knezevic, University of New Hampshire
%    email:       rem1022@wildcats.unh.edu

%fmincon - Use a single error measure for matlab's optimization toolbox,
%inputs a Hessian matrix to expedite convergence

p = inputParser;
addOptional(p, 'fixAlpha', false)
addOptional(p, 'fixPhi', false)
addOptional(p, 'maxIter', 2000)
addOptional(p, 'iGuess', 'sample')
addOptional(p, 'kernel', ones(1, maxL+1))
addOptional(p, 'Display', true)

parse(p,varargin{:})

fixAlpha = p.Results.fixAlpha;
fixPhi = p.Results.fixPhi;
maxIter = p.Results.maxIter;
initialGuess = p.Results.iGuess;
kernel = p.Results.kernel;
dispFlag = p.Results.Display;

% Scale the TDI max possible distance in an ODF
TDI_scl = grabTDIScale(CS, maxL);

fcnt = 0;


if numel(N) > 1
    %if N is actually xo
    Xo = N;
    N = numel(N)/4;
else
    %initialize phis for unconstrained solver
    if dispFlag
        fprintf('calculate initial guess\n')
    end
    Xo = MF2PHI(initialGuess, mF, CS, SS, maxL, N);
end

Xo = reshape(Xo, [4, N]);
phis = Xo(1:3, :);
alpha = Xo(4, :);

%define an initial scaling factor for alpha, mostly for matlab's steps
alpha_scl = N;

%rescale alpha
alpha = alpha./kahansum(alpha)*alpha_scl;
Xo(4, :) = alpha;

if maxIter > 0

    if fixAlpha%~fixPhi%
        if fixAlpha
            maxIterPhi = maxIter;
        else
            maxIterPhi = 100;
        end
        %first optimize Phi
        options = optimoptions('fminunc','MaxIterations',maxIterPhi,'MaxFunctionEvaluations',1000000,...
            'Display','off','SpecifyObjectiveGradient',true,'Algorithm','trust-region',...
            'HessianFcn', 'objective');

        %we optimize phis for the unconstrained solver
        x = phis(:);

        %add alpha to calcEUncon
        if dispFlag
            fprintf('phis trust-region\n')
            fprintf(' Iteration No          TDI\n')
        end
        x = fminunc(@(x) calcEJH_wrapper(x, mF, alpha, true, false, TDI_scl, kernel, CS.id, SS.id, maxL, dispFlag), x, options);
        x = reshape(x, [3, N]);

        x = [x; alpha(:)'];
        Xo = x(:);
    end



    if ~fixPhi && ~fixAlpha
        %set lower and upper bounds on alpha
        min_alpha = alpha_scl/(N*100);
        max_alpha = alpha_scl*100/(N);

        lb = zeros(4*N, 1) - Inf;
        lb(4:4:end) = 0;

        ub = zeros(4*N, 1) + Inf;
        %     ub(4:4:end) = max_alpha;

        %change stopping criterion here and with TDI_stop
        options = optimoptions('fmincon','MaxIterations',maxIter,'MaxFunctionEvaluations',1000000,...
            'Display','off','SpecifyObjectiveGradient',true,'Algorithm','trust-region-reflective',...%interior-point
            'HessianFcn','objective',...
            'FunctionTolerance', 1e-10, 'StepTolerance', 1e-10);

        if dispFlag
            fprintf('phis & alpha trust-region-reflective\n')
            fprintf(' Iteration No          TDI\n')
        end

        x = fmincon(@(x) calcEJH_wrapper(x, mF, [], false, false, TDI_scl, kernel, CS.id, SS.id, maxL, dispFlag), Xo(:), [], [], [], [], lb, ub, [], options);
    end

else

    x = Xo;
end


%reshaping x back into regular format
x = reshape(x, [4 N]);

%converts phis [x(:, 1:3)] to the fundamental zone
ori = calcOri(x(1:3, :), CS, SS);

ori = ori2eulerFZ(ori);

%rescales alpha to 1
alpha = x(4, :)';
alpha = alpha/kahansum(alpha);

%--------------------------------------------------------------------------

    function[E, J, H] = calcEJH_wrapper(x, mF, alpha, fixAlpha, fixPhi, TDI_scl, kernel, CS, SS, maxL, dispFlag)

        if fixPhi
            %             [E, J, H, TDI, STD] = calcEJH_alpha(x, in2, in3, in4);

        elseif fixAlpha
            [E, J, H] = calcGSH(3, x, alpha, mF, CS, SS, maxL, kernel);
            TDI = E;
        else
            [E, J, H] = calcGSH(2, x, mF, CS, SS, maxL, kernel);
            TDI = E;
        end

        TDI = TDI*TDI_scl;
        fcnt = fcnt + 1;

        if dispFlag
            fprintf('%8.0f  %20.4e\n', fcnt, TDI)
        end

    end

end