function [errorRateMap, tauMap, tauResult] = MC(nVertex, nBlock, ...
    dimLatentPosition, adjMatrix, muHat, sigmaHat, tauHat, rho, ...
    tauStar, nuStar, sigmaStar, theta, nBurnIn, nConverge, ...
    scaleCovariance, modelType, isHomophily, isIdentifiable)
% Inference using 1 MCMC chain.

rng shuffle;

nMetropolistHastings = 10;

%% --- Initialization ---
% Take the last NConverge iterations as the results after the burn-in part
tauResult = zeros(nConverge, nVertex);

% initialization of tau, use the estimated tauHat from asge directly
tau = tauHat;

% calculation of T (1-by-K), number of vertices in each block
nVectorHat = zeros(1, nBlock);
for iBlock = 1:nBlock
    nVectorHat(iBlock) = sum(tau == iBlock);
end

% Generate valid initialial nu (K-by-d)
nu = nuGenerator(nVertex, nBlock, dimLatentPosition, nuStar, sigmaStar, ...
    muHat, sigmaHat, scaleCovariance, modelType, isHomophily, ...
    isIdentifiable, 1);

% pre-calculation of f = <nu_i,nu_j>^p*(1-<nu_i,nu_j>)^(1-p) (KxKx2)
f = fCalculator(nu);

%% Burn-in Part
iter = 0;
while (iter < nBurnIn)
    iter = iter + 1;
    % --- Gibbs Sampling of tau ---
    [tau, nVectorHat] = tauUpdate(tau, adjMatrix, nVectorHat, theta, ...
        rho, f, nVertex, nBlock, modelType);
    % --- Metropolis-Hasting of nu ---
    if (scaleCovariance ~= 5)
        for iMetropolisHastings = 1:nMetropolistHastings
            [nu, f] = nuUpdate(nVertex, tau, nu, muHat, sigmaHat, ...
                nuStar, sigmaStar, adjMatrix, f, nBlock, ...
                dimLatentPosition, scaleCovariance, modelType, ...
                isHomophily, isIdentifiable);
        end
    end
end

%% Ergodic Average
iter = 0;
while (iter < nConverge)
    iter = iter + 1;
    % --- Gibbs Sampling of tau ---
    [tau, nVectorHat] = tauUpdate(tau, adjMatrix, nVectorHat, theta, ...
        rho, f, nVertex, nBlock, modelType);
    % --- Metropolis-Hasting of nu ---
    if (scaleCovariance ~= 5)
        for iMetropolisHastings = 1:nMetropolistHastings
            [nu, f] = nuUpdate(nVertex, tau, nu, muHat, sigmaHat, ...
                nuStar, sigmaStar, adjMatrix, f, nBlock, ...
                dimLatentPosition, scaleCovariance, modelType, ...
                isHomophily, isIdentifiable);
        end
    end
    tauResult(iter, :) = tau;
end

tauMap = mode(tauResult);

errorRateMap = nVertex;
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

errorRateMap = errorRateMap/nVertex;


