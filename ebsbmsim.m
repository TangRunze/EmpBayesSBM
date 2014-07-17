function [] = ebsbmsim(nVertex, nBlock, epsilonInB, gStart, gEnd, ...
    modelTypeStart, modelTypeEnd, scaleCovarianceStart, ...
    scaleCovarianceEnd, nCore, nBurnIn, nConverge, rho, ...
    dimLatentPosition, isHomophily, isIdentifiable, theta)

% EmpBayesSBM implements an empirical Bayes methodology for estimation of
% block memberships of vertices in a random graph drawn from the stochastic
% blockmodel.
% 
% There are two versions of EmpBayesSBM: simulation & real data.
% And this is the simulation version.
% 
% EBSBMSIM(nVertex, nBlock, epsilonInB, gStart, gEnd, ...
%     modelTypeStart, modelTypeEnd, scaleCovarianceStart,
%     scaleCovarianceEnd)
%
% nVertex selects the number of vertices in the graph.
%
% nBlock selects the number of blocks in the stochastic blockmodel.
%
% epsilonInB controls the true model. The probability matrix
%       B = (0.5 - epsilonInB)*J + 2*epsilonInB*I
% epsilonInB should be inside [0, 0.5].
%
% gStart/gEnd selects the starting/ending graph in the simulation.
%
% modelTypeStart/modelTypeEnd selects the starting/ending type of model.
% modelTypeStart/modelTypeEnd should be integers between 1 and 4.
%       modelType = 1:    Gold
%       modelType = 2:    ASGE
%       modelType = 3:    ASGE1
%       modelType = 4:    Flat
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
% EXAMPLE: 
%       EBSBMSIM(150, 3, 0.1, 1, 10, 1, 2, 4, 5)
%
% EBSBMSIM(nVertex, nBlock, epsilonInB, gStart, gEnd, ...
%     modelTypeStart, modelTypeEnd, scaleCovarianceStart, ...
%     scaleCovarianceEnd, nCore)
% 
% nCore selects the number of cores in parallel programming.
% nCore = 1 do NOT run the code parallel.
% nCore should be integers between 1 and 12.
% Default nCore = 1.
% 
% EXAMPLE:
%       EBSBMSIM(150, 3, 0.1, 1, 10, 1, 2, 4, 5, 1)     NOT parallel
%       EBSBMSIM(150, 3, 0.1, 1, 10, 1, 2, 4, 5, 12)    use 12 cores
%
% EBSBMSIM(nVertex, nBlock, epsilonInB, gStart, gEnd, ...
%     modelTypeStart, modelTypeEnd, scaleCovarianceStart, ...
%     scaleCovarianceEnd, nCore, nBurnIn, nConverge)
%
% nBurnIn selects the number of iterations in burn-in part.
% Default nBurnIn = 19000.
%
% nConverge selects the number of iterations for analysis after the burn-in
% part.
% Default nConverge = 500.
%
% EXAMPLE:
%       EBSBMSIM(150, 3, 0.1, 1, 10, 1, 2, 4, 5, 1, 10000, 1000)
%
% EBSBMSIM(nVertex, nBlock, epsilonInB, gStart, gEnd, ...
%     modelTypeStart, modelTypeEnd, scaleCovarianceStart, ...
%     scaleCovarianceEnd, nCore, nBurnIn, nConverge, rho)
%
% rho selects the block proportions.
% rho should sum up to 1.
% Default rho = [1/nBlock, 1/nBlock, ..., 1/nBlock] of length nBlock.
% 
% EXAMPLE:
%       EBSBMSIM(150, 3, 0.1, 1, 10, 1, 2, 4, 5, 1, 10000, 1000, ...
%           [1/2, 1/3, 1/6])
%
% EBSBMSIM(nVertex, nBlock, epsilonInB, gStart, gEnd, ...
%     modelTypeStart, modelTypeEnd, scaleCovarianceStart, ...
%     scaleCovarianceEnd, nCore, nBurnIn, nConverge, rho, ...
%     dimLatentPosition)
%
% dimLatentPosition selects the dimension of latent positions.
% dimLatentPosition should be integers between 1 and nBlock.
% Default dimLatentPosition = nBlock.
% 
% EXAMPLE:
%       EBSBMSIM(150, 3, 0.1, 1, 10, 1, 2, 4, 5, 1, 10000, 1000, ...
%           [1/2, 1/3, 1/6], 2)
%
% EBSBMSIM(nVertex, nBlock, epsilonInB, gStart, gEnd, ...
%     modelTypeStart, modelTypeEnd, scaleCovarianceStart, ...
%     scaleCovarianceEnd, nCore, nBurnIn, nConverge, rho, ...
%     dimLatentPosition, isHomophily, isIdentifiable)
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
%       EBSBMSIM(150, 3, 0.1, 1, 10, 1, 2, 4, 5, 1, 10000, 1000, ...
%           [1/2, 1/3, 1/6], 2, 1, 1)
% 
% EBSBMSIM(nVertex, nBlock, epsilonInB, gStart, gEnd, ...
%     modelTypeStart, modelTypeEnd, scaleCovarianceStart, ...
%     scaleCovarianceEnd, nCore, nBurnIn, nConverge, rho, ...
%     dimLatentPosition, isHomophily, isIdentifiable, theta)
% 
% theta selects the hyperparameters for the prior distribution for rho.
% Default theta = [1, 1, ..., 1] of length nBlock.
%
% EXAMPLE:
%       EBSBMSIM(150, 3, 0.1, 1, 10, 1, 2, 4, 5, 1, 10000, 1000, ...
%           [1/2, 1/3, 1/6], 2, 1, 1, [1, 3, 2])

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
% Oct 2013; Last revision: 17-July-2014

%% --- Default Parameter Setting ---
if (nargin < 17)
    % hyperparameters for the prior distribution for rho (1-by-nBlock)
    theta = ones(1, nBlock);
end

if (nargin < 16)
    % isHomophily controls the constraints. 
    % If isHomophily = 1, then it enforces homophily:
    %       <nu_i,nu_j> <= <nu_i,nu_i>  for any i, j
    isHomophily = 1;
    % isIdentifiable controls the constraints.
    % If isIdentifiable = 1, then it enforces identifiability:
    %       <nu_i,nu_i> >= <nu_j,nu_j>  for any i > j
    isIdentifiable = 0;
end

if (nargin < 14)
    % dimLatentPosition selects the dimension of latent positions.
    dimLatentPosition = nBlock;
end

if (nargin < 13)
    % true block proportion
    rho = repmat(1/nBlock, 1, nBlock);
end

if (nargin < 12)
    % number of iterations in burn-in part
    nBurnIn = 19000;
    % Take the last NConverge iterations as the results after the burn-in part
    nConverge = 500;
end

if (nargin < 10)
    % nCore selects the number of cores in parallel programming.
    % nCore = 1 do NOT run the code parallel.
    nCore = 1;
end

if ~ismember(nargin, [9, 10, 12, 13, 14, 16, 17])
    error('Invalid number of input arguments!')
end

if ((ceil(nVertex) ~= floor(nVertex)) || (nVertex <= 0))
    error('Number of vertices should be a positive integer!')
end

if ((ceil(nBlock) ~= floor(nBlock)) || (nBlock <= 0))
    error('Number of blocks should be a positive integer!')
end

if ((epsilonInB < 0) || (epsilonInB > 0.5))
    error('Epsilon should be inside [0, 0.5]!')
end

if ((ceil(gStart) ~= floor(gStart)) || (ceil(gEnd) ~= floor(gEnd)) || ...
        (gStart <= 0) || (gEnd <= 0))
    error('gStart/gEnd should be positive integers!')
end

if (gStart > gEnd)
    error('gStart should be less or equal to gEnd!')
end

if ((ceil(modelTypeStart) ~= floor(modelTypeStart)) || ...
        (ceil(modelTypeEnd) ~= floor(modelTypeEnd)) || ...
        (modelTypeStart <= 0) || (modelTypeStart > 4) || ...
        (modelTypeEnd <= 0) || (modelTypeEnd > 4))
    error(['modelTypeStart/modelTypeEnd should be positive integers' ...
        'between 1 and 4!'])
end

if (modelTypeStart > modelTypeEnd)
    error('modelTypeStart should be less or equal to modelTypeEnd!')
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

if ((abs(sum(rho) - 1) > 1e-6) || any(rho < 0))
    error('Block proportions should be nonnegative numbers sum up to 1!')
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
% if isempty(gcp('nocreate'))
%     parpool(nCore);
% end
delete(gcp('nocreate'))
parpool(nCore);

%% --- Parameters Setting ---

% block probability matrix
B = (0.5 - epsilonInB)*ones(nBlock, nBlock) + 2*epsilonInB*eye(nBlock);

% true tau_star (1-by-nVertex)
tauStar = [];
nVectorStar = nVertex*rho;
for i = 1:nBlock
    tauStar = [tauStar, i*ones(1, nVectorStar(i))];
end

% true nu (nBlock-by-dimLatentPosition)
% nuStar = chol(B)';
[U, S, ~] = svds(B, dimLatentPosition);
nuStar = U*sqrt(S);
if (~checkconstraints(nuStar, isHomophily, isIdentifiable))
    error(['The true latent positions nu derived from the probability ' ...
        'matrix B does not satisfy the constraints!']);
end

% true spectral graph embedding Xhat (nBlock-by-dimLatentPosition)
xHat = asge(B, dimLatentPosition);

% true covariance matrix Sigma_star
% (dimLatentPosition-by-dimLatentPosition-by-nBlock)
sigmaStar = covariancecalculator(xHat, rho);

%% Monte Carlo Simulation

parfor iGraph = gStart:gEnd
    % Generate data if there does not exist one, otherwise read the
    % existing data.
    [adjMatrix, muHat, sigmaHat, tauHat, pTauHat] = ...
        datagenerator(nVertex, nBlock, dimLatentPosition, B, rho, ...
        epsilonInB, iGraph);
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
    end
    
    % Run the algorithm to estimate the block membership of vertices
    for modelType = modelTypeStart:modelTypeEnd
        for scaleCovariance = scaleCovarianceStart:scaleCovarianceEnd
            saveFile = ['./results/results-SBM-model' ...
                num2str(modelType) '-scale' num2str(scaleCovariance) ...
                '-sim-n' num2str(nVertex) '-eps' num2str(epsilonInB) ...
                '-graph' num2str(iGraph) '.mat'];
            if exist(saveFile, 'file') == 0
                if (modelType == 1) || (scaleCovariance == 5)
                    [errorRate, tau, tauConfidence, tauResult] = ...
                        mcmc1chain(nVertex, nBlock, dimLatentPosition, ...
                        adjMatrix, muHat, sigmaHat, tauHat, rho, ...
                        tauStar, nuStar, sigmaStar, theta, nBurnIn, ...
                        nConverge*2, scaleCovariance, modelType, ...
                        isHomophily, isIdentifiable);
                else
                    [errorRate, tau, tauConfidence, tauResult] = ...
                        mcmc2chains(nVertex, nBlock, dimLatentPosition, ...
                        adjMatrix, muHat, sigmaHat, tauHat, pTauHat, ...
                        rho, tauStar, nuStar, sigmaStar, theta, nBurnIn,...
                        nConverge, scaleCovariance, modelType, ...
                        isHomophily, isIdentifiable);
                end
                parsave(saveFile, errorRate, tau, tauConfidence, ...
                    tauResult);
            end
        end
    end
end

%% --- Close Parallel Computing ---
delete(gcp('nocreate'))
