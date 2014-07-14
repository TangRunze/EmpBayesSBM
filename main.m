
% EmpBayesSBM implements an empirical Bayes methodology for estimation of
% block memberships of vertices in a random graph drawn from the stochastic
% blockmodel.

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
% Oct 2013; Last revision: 14-July-2014

clear all;
close all;

%% --- Parallel Computing ---
% if isempty(gcp('nocreate'))
%     parpool(12);
% end

% s = matlabpool('size');
% if s == 0
%     matlabpool(12);
% end

%% --- Parameters Setting ---

% head and tail of the Monte Carlo replicates
gStart = 1;
gEnd = 1;

% number of vertices
nVertex = 150;

% number of blocks
nBlock = 3;

% dimension of latent position nu
dimLatentPosition = nBlock;

% parameters in constraints checking
% If isHomophily = 1, then it enforces homophily:
% <nu_i,nu_j> <= <nu_i,nu_i>
% If isIdentifiable = 1, then it enforces identifiability:
% <nu_i,nu_i> >= <nu_j,nu_j> for i > j
isHomophily = 1;
isIdentifiable = 0;

% parameters in B = (0.5 - eps)*J + 2*eps*I
epsilonInB = 0.1;

% block probability matrix
B = (0.5 - epsilonInB)*ones(nBlock, nBlock) + 2*epsilonInB*eye(nBlock);

% true block proportion
rho = [1/3, 1/3, 1/3];

% hyperparameters for the prior distribution for rho (1-by-K)
theta = ones(1, nBlock);

% true tau_star (1-by-n)
tauStar = [];
nVectorStar = nVertex*rho;
for i = 1:nBlock
    tauStar = [tauStar, i*ones(1, nVectorStar(i))];
end

% true nu (K-by-d)
nuStar = chol(B)';
if (~checkconstraints(nuStar, isHomophily,0))
    error(['The true latent positions nu derived from the probability' ...
        'matrix B does not satisfy the constraints in S!']);
end

% true spectral graph embedding Xhat (K-by-d)
xHat = asge(B, dimLatentPosition);

% true covariance matrix Sigma_star (d-by-d-by-K)
sigmaStar = covariancecalculator(xHat, rho);

%% Monte Carlo Simulation

% type of scaleCovariance, which controls the size of the covariance matrix
% in the prior.
% scaleCovariance = 1:  0
% scaleCovariance = 2:  n
% scaleCovariance = 3:  n*n_k
% scaleCovariance = 4:  n^2
% scaleCovariance = 5:  Infinity
MAX_SCALE_COVARIANCE = 5;

% type of model
% modelType = 1:    Gold
% modelType = 2:    ASGE
% modelType = 3:    ASGE1
% modelType = 4:    Flat
MAX_MODEL_TYPE = 3;

% number of iterations in burn-in part
nBurnIn = 19000;
% nBurnIn = 100;

% Take the last NConverge iterations as the results after the burn-in part
nConverge = 1000;
% nConverge = 100;

for iGraph = gStart:gEnd
    % Generate data if there does not exist one, otherwise read the
    % existing data.
    [adjMatrix, muHat, sigmaHat, tauHat, pTauHat] = ...
        datagenerator(nVertex, nBlock, dimLatentPosition, B, rho, iGraph);
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
    for modelType = 2
        for scaleCovariance = 2
            savefile = ['./results/results-SBM-model' ...
                num2str(modelType) '-scale' num2str(scaleCovariance) ...
                '-graph' num2str(iGraph) '.mat'];
            if exist(savefile, 'file') == 0
                if (modelType == 1) || (scaleCovariance == 5)
                    nConverge = 1000;
                    [errorRate, tau, tauResult] = mcmc1chain(nVertex, ...
                        nBlock, dimLatentPosition, adjMatrix, muHat, ...
                        sigmaHat, tauHat, rho, tauStar, nuStar, ...
                        sigmaStar, theta, nBurnIn, nConverge, ...
                        scaleCovariance, modelType, isHomophily, ...
                        isIdentifiable);
                else
                    nConverge = 500;
                    [errorRate, tau, tauResult] = mcmc2chains(nVertex, ...
                        nBlock, dimLatentPosition, adjMatrix, muHat, ...
                        sigmaHat, tauHat, pTauHat, rho, tauStar, nuStar,...
                        sigmaStar, theta, nBurnIn, nConverge, ...
                        scaleCovariance, modelType, isHomophily, ...
                        isIdentifiable);
                end
                parsave(savefile, errorRate, tau, tauResult);
            end
        end
    end
end

% delete(gcp('nocreate'))
% matlabpool close
