function[phis, alpha] = fourier2ori_max(mF, pF, CS, SS, maxL, psi_scl, psi_inv, dpow, X0)
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

maxIter = 250;

N = numel(X0)/4;
X0 = reshape(X0, [4, N]);
phis = X0(1:3, :);
alpha = X0(4, :);

%Only output every step if code is not executing in parallel
isOnWorker = ~isempty(getCurrentTask());


%initialize function variable
fcnt = 0;
%first use unconstrained solver
options = optimoptions('fminunc','MaxIterations',maxIter,'MaxFunctionEvaluations',1000000,...
    'Display','off','SpecifyObjectiveGradient',true,'Algorithm','trust-region',...
    'HessianFcn', 'objective');

%we optimize phis for the unconstrained solver
x = phis(:);

%add alpha to calcEUncon
if ~isOnWorker
    fprintf('phis trust-region\n')
    fprintf(' Iteration No          DIFF             SCL\n')
end
% options = optimoptions('fminunc','MaxIterations',10,'MaxFunctionEvaluations',1000000,...
%     'Display','off','SpecifyObjectiveGradient',true,'Algorithm','trust-region',...
%     'HessianFcn', 'objective');
% 
% [x] = fminunc(@(x) calcEJH_wrapper(x, alpha, mF, CS, SS, maxL, psi_scl, psi_inv, dpow, pF, true, isOnWorker), x, options);

options = optimoptions('fminunc','MaxIterations',maxIter,'MaxFunctionEvaluations',1000000,...
    'Display','off','SpecifyObjectiveGradient',true,'Algorithm','trust-region',...
    'HessianFcn', 'objective');

[x] = fminunc(@(x) calcEJH_wrapper(x, alpha, mF, CS, SS, maxL, psi_scl, psi_inv, dpow, pF, false, isOnWorker), x, options);


%reshaping x back into regular format
x = reshape(x, [3, numel(x)/3]);

phis = x(1:3, :);

%rescales alpha to 1
alpha = alpha/kahansum(alpha);

%--------------------------------------------------------------------------

    function[E, J, H] = calcEJH_wrapper(x, alpha, mF, CS, SS, maxL, psi_scl, psi_inv, dpow, pF, iFirst, isOnWorker)

        r = -20;
        K = 1e-6;
        if false%iFirst
            % Set the psi_inv to zero
            [E, J, H] = calcGSH(3, x, alpha, mF, CS, SS, maxL, psi_inv);

            E1 = 0;
            E2 = E;

        else
            % Set the difference to the max
            [E1, J1, H1] = calcGSH(7, x, alpha, pF, CS, SS, maxL, r, psi_scl);
            % Set the psi_inv to zero
            [E2, J2, H2] = calcGSH(3, x, alpha, mF, CS, SS, maxL, psi_inv);

            n1 = -dpow;
            n2 = 1;
            %E = K1*E1 + K2*E2;
            %J = K1*J1 + K2*J2;
            %H = K1*H1 + K2*H2;
            E = (E1)^n1 * (K + E2)^n2;
            J = J1.*(n1*E1^(n1 - 1)*(K + E2)^n2) + J2.*(n2*E1^n1*(K + E2)^(n2 - 1));
            H = (n1*E1^(n1 - 1)*(K + E2)^n2).*H1 + (n2*E1^n1*(K + E2)^(n2 - 1)).*H2 + n1*E1^(n1 - 2)*(n1 - 1)*(K + E2)^n2*(J1'.*J1) + (n1*n2*E1^(n1 - 1)*(K + E2)^(n2 - 1))*(J2'.*J1 + J1'.*J2) + (J2'.*J2)*(n2*E1^n1*(n2 - 1)*(K + E2)^(n2 - 2));
%             %H = (n1*E1^(n1 - 1)*(K + E2)^n2).*H1 + (n2*E1^n1*(K + E2)^(n2 - 1)).*H2 + n1*E1^(n1 - 2)*(n1 - 1)*(K + E2)^n2*(J1'.*J1) + n1*n2*E1^(n1 - 1)*(K + E2)^(n2 - 1)*(J2'.*J1) + (n1*n2*E1^(n1 - 1)*(K + E2)^(n2 - 1))*(J1'.*J2) + (J2'.*J2)*(n2*E1^n1*(n2 - 1)*(K + E2)^(n2 - 2));
        end
        fcnt = fcnt + 1;

        if ~isOnWorker%mod(fcnt, IterOut)==0
            fprintf('%8.0f %20.4e %20.4e\n', fcnt, E1, E2)
        end

    end

end

%--------------------------------------------------------------------------