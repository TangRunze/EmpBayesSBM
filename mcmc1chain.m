function [errorRateMap, tauMap, tauConfidence, tauResult] = ...
    mcmc1chain(nVertex, nBlock, dimLatentPosition, adjMatrix, muHat, ...
    sigmaHat, tauHat, rho, tauStar, nuStar, sigmaStar, theta, nBurnIn, ...
    nConverge, scaleCovariance, modelType, isHomophily, isIdentifiable)
% Inference using 1 MCMC chain.

rng shuffle;

nMetropolistHastings = 10;

%% --- Initialization ---
% Take the last NConverge iterations as the results after the burn-in part
tauResult = zeros(nConverge, nVertex);

% calculation of T0 (1-by-K), true number of vertices in each block
nVectorHat0 = 0;
if (tauStar ~= 0)
    nVectorHat0 = zeros(1, nBlock);
    for iBlock = 1:nBlock
        nVectorHat0(iBlock) = sum(tauStar == iBlock);
    end
end

% initialization of tau, use the estimated tauHat from asge directly
tau = tauHat;

% calculation of T (1-by-K), number of vertices in each block
nVectorHat = zeros(1, nBlock);
for iBlock = 1:nBlock
    nVectorHat(iBlock) = sum(tau == iBlock);
end

% Generate valid initialial nu (K-by-d)
nu = nugenerator(nVertex, nBlock, dimLatentPosition, nuStar, sigmaStar, ...
    muHat, sigmaHat, scaleCovariance, modelType, isHomophily, ...
    isIdentifiable, 1, nVectorHat0);

% pre-calculation of f = <nu_i,nu_j>^p*(1-<nu_i,nu_j>)^(1-p) (KxKx2)
f = fcalculator(nu);

%% Burn-in Part
iter = 0;
while (iter < nBurnIn)
    iter = iter + 1;
    % --- Gibbs Sampling of tau ---
    [tau, nVectorHat] = updatetau(tau, adjMatrix, nVectorHat, theta, ...
        rho, f, nVertex, nBlock, modelType);
    % --- Metropolis-Hasting of nu ---
    if (scaleCovariance ~= 5)
        for iMetropolisHastings = 1:nMetropolistHastings
            [nu, f] = updatenu(nVertex, tau, nu, muHat, sigmaHat, ...
                nuStar, sigmaStar, adjMatrix, f, nBlock, ...
                dimLatentPosition, scaleCovariance, modelType, ...
                isHomophily, isIdentifiable, nVectorHat0);
        end
    end
end

%% Ergodic Average
iter = 0;
while (iter < nConverge)
    iter = iter + 1;
    % --- Gibbs Sampling of tau ---
    [tau, nVectorHat] = updatetau(tau, adjMatrix, nVectorHat, theta, ...
        rho, f, nVertex, nBlock, modelType);
    % --- Metropolis-Hasting of nu ---
    if (scaleCovariance ~= 5)
        for iMetropolisHastings = 1:nMetropolistHastings
            [nu, f] = updatenu(nVertex, tau, nu, muHat, sigmaHat, ...
                nuStar, sigmaStar, adjMatrix, f, nBlock, ...
                dimLatentPosition, scaleCovariance, modelType, ...
                isHomophily, isIdentifiable, nVectorHat0);
        end
    end
    tauResult(iter, :) = tau;
end

tauMap = mode(tauResult);

% Calculate confidence of estimation.
tauConfidence = sum(tauResult == repmat(tauMap, nConverge, 1))/nConverge;

errorRateMap = nVertex;

if (tauStar ~= 0)
    permutation = perms(1:nBlock);
    for iFactorial = 1:factorial(nBlock)
        position = permutation(iFactorial, :);
        tmpTau = tauMap;
        for jBlock = 1:nBlock
            nv = (tauMap == position(jBlock));
            tmpTau(nv) = jBlock;
        end
        if sum(tauStar ~= tmpTau) < errorRateMap
            errorRateMap = sum(tauStar ~= tmpTau);
        end
    end
end

errorRateMap = errorRateMap/nVertex;


