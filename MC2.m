function [errorRateMap, tauMap, tauResult] = MC2(nVertex, nBlock, ...
    dimLatentPosition, adjMatrix, muHat, sigmaHat, tauHat, ...
    pTauHat, rho, tauStar, nuStar, sigmaStar, theta, nBurnIn, ...
    nConverge, scaleCovariance, modelType, isHomophily, isIdentifiable)
% Inference using 2 MCMC chains.

rng shuffle;

nMetropolistHastings = 10;

%% Metropolis-Hasting within Gibbs Sampling

isConverge = 0;
tmpTau1 = zeros(nConverge, nVertex);
tmpTau2 = zeros(nConverge, nVertex);
Q = zeros(nConverge, 2);
nboot = nConverge;

maxIteration = nBurnIn + nConverge;
iterOut = 0;
isRestart = 0;

while (~isConverge) || (iterOut <= nConverge)
    % If it is the beginning of the chain, initialize the parameters. Chain
    % 1 is using a fixed initialization, so when restart=1, it is better to
    % use what we have as initial.
    if (iterOut == 0)
        % Initialization of chain 2 based on the EM posterior prob
        tmp = mnrnd(1, pTauHat');
        tau2 = nBlock + 1 - sum(cumsum(tmp'));
        
        % calculation of T2 (1-by-K), number of vertices in each block for
        % the 2nd chain
        nVectorHat2 = zeros(1, nBlock);
        for iBlock = 1:nBlock
            nVectorHat2(iBlock) = sum(tau2 == iBlock);
        end
        
        % Generate valid initialial nu2 (K-by-d) for chain 2
        nu2 = nuGenerator(nVertex, nBlock, dimLatentPosition, nuStar, ...
            sigmaStar, muHat, sigmaHat, scaleCovariance, modelType, ...
            isHomophily, isIdentifiable, 1);
        
        % pre-calculation of
        % f2 = <nu2_i,nu2_j>^p*(1-<nu2_i,nu2_j>)^(1-p) (KxKx2)
        % for chain 2
        f2 = fCalculator(nu2);
        
        if (isRestart == 0)
            % Initialization of chain 1 based on asge estimator
            tau1 = tauHat;
            
            % calculation of T1 (1-by-K), number of vertices in each block
            % for the 1st chain
            nVectorHat1 = zeros(1,nBlock);
            for iBlock = 1:nBlock
                nVectorHat1(iBlock) = sum(tau1 == iBlock);
            end
            
            % Generate valid initialial nu1 (K-by-d) for chain 1
            nu1 = nuGenerator(nVertex, nBlock, dimLatentPosition, ...
                nuStar, sigmaStar, muHat, sigmaHat, scaleCovariance, ...
                modelType, isHomophily, isIdentifiable, 1);
            
            % pre-calculation of
            % f1 = <nu1_i,nu1_j>^p*(1-<nu1_i,nu1_j>)^(1-p) (KxKx2)
            % for chain 1
            f1 = fCalculator(nu1);
        end
    end

    iterIn = 0;
    while (iterIn < nConverge)
        iterIn = iterIn + 1;
        % --- Gibbs Sampling of tau ---
        [tau1, nVectorHat1] = tauUpdate(tau1, adjMatrix, nVectorHat1, ...
            theta, rho, f1, nVertex, nBlock, modelType);
        [tau2, nVectorHat2] = tauUpdate(tau2, adjMatrix, nVectorHat2, ...
            theta, rho, f2, nVertex, nBlock, modelType);
        % --- Metropolis-Hasting of nu ---
        if scaleCovariance~=5
            for iMetropolistHastings = 1:nMetropolistHastings
                [nu1, f1] = nuUpdate(nVertex, tau1, nu1, muHat, ...
                    sigmaHat, nuStar, sigmaStar, adjMatrix, f1, nBlock, ...
                    dimLatentPosition, scaleCovariance, modelType, ...
                    isHomophily, isIdentifiable);
                [nu2, f2] = nuUpdate(nVertex, tau2, nu2, muHat, ...
                    sigmaHat, nuStar, sigmaStar, adjMatrix, f2, nBlock, ...
                    dimLatentPosition, scaleCovariance, modelType, ...
                    isHomophily, isIdentifiable);
            end
        end
        tmpTau1(iterIn, :) = tau1;
        tmpTau2(iterIn, :) = tau2;
    end
    iterOut = iterOut + nConverge;
    
    % --- GR test ---
    Q(:, 1) = mean(repmat(tauStar, nConverge, 1) ~= tmpTau1, 2);
    Q(:, 2) = mean(repmat(tauStar, nConverge, 1) ~= tmpTau2, 2);
    % Potential Scale Reduction Factor
    % Rhat = psrf(Q(:,1),Q(:,2));
    bootRhat = bootstrap(Q(:,1),Q(:,2), nboot);
    bootRhat = sort(bootRhat);
    if bootRhat(0.95*nboot) < 1.1
        isConverge = 1;
    elseif (bootRhat(0.95*nboot) >= 1.1) && (iterOut >= maxIteration)
        iterOut = 0;
        isRestart = 1;
    end
end

%% Calculate Result

tauResult = tmpTau1;

tauMap = mode(tauResult);

errorRateMap = nVertex;
permutation = perms(1:nBlock);
for iFactorial = 1:factorial(nBlock)
    pos = permutation(iFactorial,:);
    tau_tmp = tauMap;
    for jBlock = 1:nBlock
        nv = (tauMap == pos(jBlock));
        tau_tmp(nv) = jBlock;
    end
    if sum(tauStar ~= tau_tmp) < errorRateMap
        errorRateMap = sum(tauStar ~= tau_tmp);
    end
end

errorRateMap = errorRateMap/nVertex;


