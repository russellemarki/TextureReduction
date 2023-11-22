function[iGuess] = MF2PHI(method, mF, CS, SS, maxL, N, phi_grid)

switch method
    case 'grid'
        [iGuess] = MF2PHI_grid(mF, CS, SS, maxL, N);
    case 'sample'
        [iGuess] = MF2PHI_sample(mF, CS, SS, maxL, N);
    case 'input_sample'
        if nargin < 7
            phi_grid = grabPhiGrid(CS);
        end
        iGuess = MF2PHI_input_sample(mF, CS, SS, maxL, N, phi_grid);
end

end

function[iGuess] = MF2PHI_grid(mF, CS, SS, maxL, N)
%get the L values, the reduction indices, and the negative indices

NN = N*1.1;

% the global grid
S3G_global = equispacedSO3Grid(CS, SS, 'points', NN);

S3G_global = S3G_global(:);

first = true;

while(numel(S3G_global.phi1) <= N || first)
    if first
        first = false;
    else
        while(numel(S3G_global.phi1) < N*1.1)
            NN = NN + nnz(~ind);
            
            S3G_global = equispacedSO3Grid(CS, SS, 'points', NN);
            
            S3G_global = S3G_global(:);
        end
    end
    
    phi_grid = S3G_global.Euler;
    phi_grid = phi_grid';
    
    %calc F for entire grid
    F = calcGSH(9, phi_grid, CS.id, SS.id, maxL, 2*(0:maxL) + 1);
    
    
    %calculate probability density
    pd = sum(F.*mF);
    
    ind = pd > 0;
    
    pd = pd(ind);
    
    S3G_global = S3G_global(ind);
end

[~, ind] = sort(pd);

ind = ind(numel(pd) - N + 1 : end);

pd = pd(ind);

ori = S3G_global(ind);

ori = ori2eulerFZ(ori);

%convert to Euler angles
phi_out = calcPhi(ori);

alpha_out = reshape(pd./sum(pd), [1, numel(pd)]);

iGuess = [phi_out; alpha_out];
end

function[iGuess] = MF2PHI_sample(mF, CS, SS, maxL, N)
%get the L values, the reduction indices, and the negative indices

% the global grid
S3G_global = equispacedSO3Grid(CS, SS, 'points', 10*N);%'resolution', 5*degree);%

if numel(S3G_global.phi1) < N
    S3G_global = equispacedSO3Grid(CS, SS, 'resolution', 2*degree);
end

S3G_global = S3G_global(:);

S3G_global = ori2eulerFZ(S3G_global);

phi_grid = calcPhi(S3G_global);

%calc F for entire grid
F = calcGSH(9, phi_grid, CS.id, SS.id, maxL, 2*(0:maxL) + 1);

%calculate probability density
pd = sum(F.*mF, 1)';

%only take real part
pd = real(pd);

ind = pd > 0;

pd = pd(ind);

S3G_global = S3G_global(ind);

if false
    ind = randperm(length(pd));
else
    [~, ind] = sort(pd);
    ind = [flipud(ind(1:2:end));ind(2:2:end)];
end

pd = pd(ind);

S3G_global = S3G_global(ind);

%discretely sample the probability density
[ind, w] = discretesample(pd, N);

%take those probabilities
ori = S3G_global(ind);

phi_out = calcPhi(ori);
alpha_out = reshape(w./sum(w), [1 numel(w)]);

iGuess = [phi_out; alpha_out];
end

function[iGuess] = MF2PHI_input_sample(mF, CS, SS, maxL, N, phi_grid)

%calc F for entire grid
F = calcGSH(9, phi_grid, CS.id, SS.id, maxL, 2*(0:maxL) + 1);

%calculate probability density
pd = sum(F.*mF, 1)';

%only take real part
pd = real(pd);

ind = pd > 0;

pd = pd(ind);

phi_grid = phi_grid(:, ind);

ind = randperm(length(pd));

pd = pd(ind);

phi_grid = phi_grid(:, ind);

%discretely sample the probability density
[ind, w] = discretesample(pd, N);

phi_out = phi_grid(:, ind);

alpha_out = reshape(w./sum(w), [1, numel(w)]);

iGuess = [phi_out; alpha_out];
end


function [x, w] = discretesample(p, N,varargin)
% Samples from a discrete distribution
%
%   x = discretesample(p, n)
%       independently draws n samples (with replacement) from the
%       distribution specified by p, where p is a probability array
%       whose elements sum to 1.
%
%       Suppose the sample space comprises K distinct objects, then
%       p should be an array with K elements. In the output, x(i) = k
%       means that the k-th object is drawn at the i-th trial.
%
%       I have changed this to be more specific to material science
%
%   Remarks
%   -------
%       - This function is mainly for efficient sampling in non-uniform
%         distribution, which can be either parametric or non-parametric.
%
%       - The function is implemented based on histc, which has been
%         highly optimized by mathworks. The basic idea is to divide
%         the range [0, 1] into K bins, with the length of each bin
%         proportional to the probability mass. And then, n values are
%         drawn from a uniform distribution in [0, 1], and the bins that
%         these values fall into are picked as results.
%
%       - This function can also be employed for continuous distribution
%         in 1D/2D dimensional space, where the distribution can be
%         effectively discretized.
%
%       - This function can also be useful for sampling from distributions
%         which can be considered as weighted sum of "modes".
%         In this type of applications, you can first randomly choose
%         a mode, and then sample from that mode. The process of choosing
%         a mode according to the weights can be accomplished with this
%         function.
%
%   Examples
%   --------
%       % sample from a uniform distribution for K objects.
%       p = ones(1, K) / K;
%       x = discretesample(p, n);
%
%       % sample from a non-uniform distribution given by user
%       x = discretesample([0.6 0.3 0.1], n);
%
%       % sample from a parametric discrete distribution with
%       % probability mass function given by f.
%       p = f(1:K);
%       x = discretesample(p, n);
%

%   Created by Dahua Lin, On Oct 27, 2008
%   Poorly modified by Russell Marki, On October 19, 2020; check mtex for the original
%

% parse and verify input arguments

assert(isfloat(p), 'discretesample:invalidarg', ...
    'p should be an array with floating-point value type.');

assert(isnumeric(N) && isscalar(N) && N >= 0 && N == fix(N), ...
    'discretesample:invalidarg', ...
    'n should be a nonnegative integer scalar.');

% process p if necessary

K = numel(p);
if ~isequal(size(p), [1, K])
    p = reshape(p, [1, K]);
end

% construct the bins

edges = [0, cumsum(p)];
s = edges(end);
if abs(s - 1) > eps
    edges = edges * (1 / s);
end

%10-19-2020
%create a function that would find the number of times necessary to sample
%the bins for the requisite #of crystals

rv = linspace(0 + 1/(2*N), 1 - 1/(2*N), N);
% rv = rand(1, N);
% % % rv = rv(2:end-1);

c = histc(rv, edges);
ce = c(end);
c = c(1:end-1);
c(end) = c(end) + ce;

xv = find(c);

%the weights of each crystal
w = c(xv);

nout = numel(xv) - N;

%split the necessary crystals in two
[~, ind] = sort(w);

ind = ind(end + nout + 1 : end);

wcut = floor(w(ind)/2);

%subtract the weights from the crystals
w(ind) = w(ind) - wcut;

%add the crystals as separate indices
x = [xv, xv(ind)];
w = [w, wcut];
end
