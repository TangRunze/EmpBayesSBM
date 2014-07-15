function [] = ebsbm(nBlock, gStart, gEnd, scaleCovarianceStart, ...
    scaleCovarianceEnd, hasLabel, nCore, nBurnIn, nConverge, ...
    dimLatentPosition, isHomophily, isIdentifiable, theta)

% EmpBayesSBM implements an empirical Bayes methodology for estimation of
% block memberships of vertices in a random graph drawn from the stochastic
% blockmodel.
% 
% There are two versions of EmpBayesSBM: simulation & real data.
% And this is the real data version.
% 
% EBSBM(nBlock, gStart, gEnd, scaleCovarianceStart, scaleCovarianceEnd, ...
%     hasLabel)
%
% nBlock selects the number of blocks in the stochastic blockmodel.
%
% gStart/gEnd selects the starting/ending graph in the simulation.
%
% scaleCovarianceStart/scaleCovarianceEnd selects the starting/ending scale
% of the covariance matrices.
% scaleCovarianceStart/scaleCovarianceEnd should be integers between 1 and 
% 5.
%       scaleCovariance = 1:  0
%       scaleCovariance = 2:  n
%       scaleCovariance = 3:  n*n_k
%       scaleCovariance = 4:  n^2
%       scaleCovariance = 5:  Infinity
% 
% hasLabel = 1 means true labels are known. Otherwise truelabels are
% unknown.
% Default hasLabel = 0.
%
% EXAMPLE: 
%       EBSBM(3, 1, 10, 4, 5, 1)
%
% EBSBM(nBlock, gStart, gEnd, scaleCovarianceStart, scaleCovarianceEnd, ...
%     hasLabel, nCore)
% 
% nCore selects the number of cores in parallel programming.
% nCore = 1 do NOT run the code parallel.
% nCore should be integers between 1 and 12.
% Default nCore = 1.
% 
% EXAMPLE:
%       EBSBM(3, 1, 10, 4, 5, 1, 1)     NOT parallel
%       EBSBM(3, 1, 10, 4, 5, 1, 12)    use 12 cores
%
% EBSBM(nBlock, gStart, gEnd, scaleCovarianceStart, scaleCovarianceEnd, ...
%     hasLabel, nCore, nBurnIn, nConverge)
%
% nBurnIn selects the number of iterations in burn-in part.
% Default nBurnIn = 19000.
%
% nConverge selects the number of iterations for analysis after the burn-in
% part.
% Default nConverge = 500.
%
% EXAMPLE:
%       EBSBM(3, 1, 10, 4, 5, 1, 1, 10000, 1000)
%
% EBSBM(nBlock, gStart, gEnd, scaleCovarianceStart, scaleCovarianceEnd, ...
%     hasLabel, nCore, nBurnIn, nConverge, dimLatentPosition)
%
% dimLatentPosition selects the dimension of latent positions.
% dimLatentPosition should be integers between 1 and nBlock.
% Default dimLatentPosition = nBlock.
% 
% EXAMPLE:
%       EBSBM(3, 1, 10, 4, 5, 1, 1, 10000, 1000, 2)
%
% EBSBM(nBlock, gStart, gEnd, scaleCovarianceStart, scaleCovarianceEnd, ...
%     hasLabel, nCore, nBurnIn, nConverge, dimLatentPosition, ...
%     isHomophily, isIdentifiable)
%
% isHomophily controls the constraints. 
% If isHomophily = 1, then it enforces homophily:
%       <nu_i,nu_j> <= <nu_i,nu_i>  for any i, j
% Default isHomophily = 1.
%
% isIdentifiable controls the constraints.
% If isIdentifiable = 1, then it enforces identifiability:
%       <nu_i,nu_i> >= <nu_j,nu_j>  for any i > j
% Default isIdentifiable = 0.
%
% EXAMPLE:
%       EBSBM(3, 1, 10, 4, 5, 1, 1, 10000, 1000, 2, 1, 1)
% 
% EBSBM(nBlock, gStart, gEnd, scaleCovarianceStart, scaleCovarianceEnd, ...
%     hasLabel, nCore, nBurnIn, nConverge, dimLatentPosition, ...
%     isHomophily, isIdentifiable, theta)
% 
% theta selects the hyperparameters for the prior distribution for rho.
% Default theta = [1, 1, ..., 1] of length nBlock.
%
% EXAMPLE:
%       EBSBM(3, 1, 10, 4, 5, 1, 1, 10000, 1000, 2, 1, 1, [1, 3, 2])

% Inference for the stochastic blockmodel is an interesting direction in
% the statistical community, as well as in various application domains such
% as social networks, citation networks, brain connectivity networks, etc.
% Recent theoretical developments have shown that a random dot product
% latent position graph formulation of the stochastic blockmodel informs a
% mixture of normal distributions for the adjacency spectral embedding.
% We employ this new theory to provide an empirical Bayes methodology for
% estimation of block memberships of vertices in a random graph drawn from
% the stochastic blockmodel. The posterior inference is conducted using a
% Metropolis-within-Gibbs algorithm.

% Author: Runze Tang
% Johns Hopkins University
% Email: tangrunze@gmail.com
% Website: https://github.com/wonderforyou/EmpBayesSBM
% Oct 2013; Last revision: 15-July-2014

%% --- Default Parameter Setting ---
if (nargin < 13)
    % hyperparameters for the prior distribution for rho (1-by-nBlock)
    theta = ones(1, nBlock);
end

if (nargin < 12)
    % isHomophily controls the constraints. 
    % If isHomophily = 1, then it enforces homophily:
    %       <nu_i,nu_j> <= <nu_i,nu_i>  for any i, j
    isHomophily = 1;
    % isIdentifiable controls the constraints.
    % If isIdentifiable = 1, then it enforces identifiability:
    %       <nu_i,nu_i> >= <nu_j,nu_j>  for any i > j
    isIdentifiable = 0;
end

if (nargin < 10)
    % dimLatentPosition selects the dimension of latent positions.
    dimLatentPosition = nBlock;
end

if (nargin < 9)
    % number of iterations in burn-in part
    nBurnIn = 19000;
    % Take the last NConverge iterations as the results after the burn-in part
    nConverge = 500;
end

if (nargin < 7)
    % nCore selects the number of cores in parallel programming.
    % nCore = 1 do NOT run the code parallel.
    nCore = 1;
end

if ~ismember(nargin, [6, 7, 9, 10, 12, 13])
    error('Invalid number of input arguments!')
end

if ((ceil(nBlock) ~= floor(nBlock)) || (nBlock <= 0))
    error('Number of blocks should be a positive integer!')
end

if ((ceil(gStart) ~= floor(gStart)) || (ceil(gEnd) ~= floor(gEnd)) || ...
        (gStart <= 0) || (gEnd <= 0))
    error('gStart/gEnd should be positive integers!')
end

if (gStart > gEnd)
    error('gStart should be less or equal to gEnd!')
end

if ((ceil(scaleCovarianceStart) ~= floor(scaleCovarianceStart)) || ...
        (ceil(scaleCovarianceEnd) ~= floor(scaleCovarianceEnd)) || ...
        (scaleCovarianceStart <= 0) || (scaleCovarianceStart > 5) || ...
        (scaleCovarianceEnd <= 0) || (scaleCovarianceEnd > 5))
    error(['scaleCovarianceStart/scaleCovarianceEnd should be' ...
        'positive integers between 1 and 5!'])
end

if (scaleCovarianceStart > scaleCovarianceEnd)
    error(['scaleCovarianceStart should be less or equal to' ...
        'scaleCovarianceEnd!'])
end

if ((ceil(nCore) ~= floor(nCore)) || (nCore <= 0) || (nCore >12))
    error('Number of cores should be a positive integer between 1 and 12!')
end

if ((ceil(nBurnIn) ~= floor(nBurnIn)) || (nBurnIn <= 0))
    error(['Number of iterations in burn-in part should be' ...
        'a positive integer!'])
end

if ((ceil(nConverge) ~= floor(nConverge)) || (nConverge <= 0))
    error(['Number of iterations for analysis after burn-in part' ...
        'should be a positive integer!'])
end

if ((ceil(dimLatentPosition) ~= floor(dimLatentPosition)) || ...
        (dimLatentPosition <= 0) || (dimLatentPosition > nBlock))
    error(['Dimension of latent positions should be positive integers' ...
        'and less or equal to number of blocks!'])
end

if ((isHomophily ~= 0) && (isHomophily ~= 1))
    error('isHomophily should be either 0 or 1!')
end

if ((isIdentifiable ~= 0) && (isIdentifiable ~= 1))
    error('isIdentifiable should be either 0 or 1!')
end

if (any(theta <= 0))
    error('Hyperparameters should be positive!')
end

%% --- Parallel Computing ---
if isempty(gcp('nocreate'))
    parpool(nCore);
end

%% --- Parameters Setting ---
% modelType selects the type of model.
% In practise, modelType could only be 2.
%       modelType = 1:    Gold
%       modelType = 2:    ASGE
%       modelType = 3:    ASGE1
%       modelType = 4:    Flat
modelType = 2;

%% Monte Carlo Simulation

for iGraph = gStart:gEnd
    % Generate data if there does not exist one, otherwise read the
    % existing data.
    [nVertex, adjMatrix, muHat, sigmaHat, tauHat, pTauHat, tauStar] = ...
        datareader(nBlock, dimLatentPosition, iGraph, hasLabel);
    % nVertex selects the number of vertices in the graph
    % Adjacency matrix A (n-by-n) symmetric
    % cluster means mu_hat (K-by-d)
    % cluster covariances Sigma_hat (d-by-d-by-K)
    % classification tau_hat (1-by-n)
    % posterior p_tau_hat (K-by-n)
    
    % Reorder the blocks when they do not satisfy the identifiability
    % constraints.
    [muHat, sigmaHat, tauHat, pTauHat] = changeorder(muHat, sigmaHat, ...
        tauHat, pTauHat, nBlock);
    
    % Project the estimated mu_hat onto the neareast point in the feasible
    % region.
    if (~checkconstraints(muHat, isHomophily, isIdentifiable))
        tmpMuHat = muHat';
        tmpMuHat = tmpMuHat(:);
        tmpMuHat = fmincon(@(x) projectionobjectivefun(x,tmpMuHat), ...
            tmpMuHat, [], [], [], [], ...
            - ones(dimLatentPosition*nBlock, 1), ...
            ones(dimLatentPosition*nBlock, 1), @(x) ...
            projectionconditionfun(x, nBlock, dimLatentPosition, ...
            isHomophily, isIdentifiable));
        muHat = reshape(tmpMuHat, [dimLatentPosition, nBlock])';
        clear tmpMuHat;
    end
    
    % Run the algorithm to estimate the block membership of vertices
    for scaleCovariance = scaleCovarianceStart:scaleCovarianceEnd
        savefile = ['./results/results-SBM-model' ...
            num2str(modelType) '-scale' num2str(scaleCovariance) ...
            '-graph' num2str(iGraph) '.mat'];
        if exist(savefile, 'file') == 0
            if (modelType == 1) || (scaleCovariance == 5)
                [errorRate, tau, tauResult] = mcmc1chain(nVertex, ...
                    nBlock, dimLatentPosition, adjMatrix, muHat, ...
                    sigmaHat, tauHat, 0, tauStar, 0, 0, theta, nBurnIn, ...
                    nConverge*2, scaleCovariance, modelType, ...
                    isHomophily, isIdentifiable);
            else
                [errorRate, tau, tauResult] = mcmc2chains(nVertex, ...
                    nBlock, dimLatentPosition, adjMatrix, muHat, ...
                    sigmaHat, tauHat, pTauHat, 0, tauStar, 0, 0, theta, ...
                    nBurnIn, nConverge, scaleCovariance, modelType, ...
                    isHomophily, isIdentifiable);
            end
            parsave(savefile, errorRate, tau, tauResult);
        end
    end
end

%% --- Close Parallel Computing ---
delete(gcp('nocreate'))
