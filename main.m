
% EmpBayesSBM implements an empirical Bayes methodology for estimation of block
% memberships of vertices in a random graph drawn from the stochastic
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
gstart = 1;
gend = 10;

% number of points
n = 150;

% number of blocks
K = 3;

% dimension of nu
d = K;

% number of monte carlo replicates
G = 100;

% constraints parameters
homophily = 1;
identifiability = 0;

% epsilon parameters in B = (0.5 - eps)*J + 2*eps*I
eps = 0.1;

% block probability matrix
B = (0.5 - eps)*ones(K,K) + 2*eps*eye(K);

% true block proportion
rho = [1/3, 1/3, 1/3];

% hyperparameters for the prior distribution for rho (1-by-K)
theta = ones(1,K);

% true tau_star (1-by-n)
tau_star = [];
ni = n*rho;
for i = 1:K
    tau_star = [tau_star, i*ones(1,ni(i))];
end

% true nu (K-by-d)
nu_star = chol(B)';
if ~CheckS(nu_star,homophily,0)
    error('The true nu does not satisfy the constraints in S');
end

% true Sigma_star (d-by-d-by-K)
Xhat = asge(B,K);
Sigma_star = CovarianceCalculator(Xhat,rho);

%% Monte Carlo Simulation

% type of c
% c = 1:  0
% c = 2:  n
% c = 3:  n*n_k
% c = 4:  n^2
% c = 5:  Infinity
MaxC = 5;

% type of model
% model = 1:    Gold
% model = 2:    ASGE
% model = 3:    ASGE1
% model = 4:    Flat
MaxModel = 3;

% number of Burn-In part
NBurnIn = 19000;
% NBurnIn = 10;

% Take the last NConverge iterations as the results after the burn-in part
NConverge = 1000;
% NGS = 10;

for g = gstart:gend

    [A,mu_hat,Sigma_hat,tau_hat,p_tau_hat] = DataGenerator(n,K,d,B,rho,g);
    % Adjacency matrix A (n-by-n) symmetric
    % cluster means mu_hat (K-by-d)
    % cluster covariances Sigma_hat (d-by-d-by-K)
    % classification tau_hat (1-by-n)
    % posterior p_tau_hat (K-by-n)
    
    % Reorder the blocks when they do not satisfy the identifiability
    % constraints.
    [mu_hat,Sigma_hat,tau_hat,p_tau_hat] = ...
        ChangeOrder(mu_hat,Sigma_hat,tau_hat,p_tau_hat,K);
    
    % Project the estimated mu_hat onto the neareast point in the feasible
    % region.
    if (~CheckS(mu_hat,homophily,identifiability))
        tmp = mu_hat';
        tmp = tmp(:);
        tmp = fmincon(@(x) ProjectObj(x,tmp),tmp,[],[],[],[],...
            -ones(d*K,1),ones(d*K,1),...
            @(x) ProjectCon(x,K,d,homophily,identifiability));
        mu_hat = reshape(tmp,[d,K])';
    end
    
    for model = 2
        for c = 1:2
            savefile = ['./results/results-SBM-model' num2str(model) ...
                '-c' num2str(c) '-graph' num2str(g) '.mat'];
            if exist(savefile,'file') == 0
                if (model == 1) || (c == 5)
                    NConverge = 1000;
                    [error_rate,tau,tau_result] = MC(n,K,d,A,mu_hat,...
                        Sigma_hat,tau_hat,rho,tau_star,nu_star,...
                        Sigma_star,theta,NBurnIn,NConverge,c,model,...
                        homophily,identifiability);
                else
                    NConverge = 500;
                    [error_rate,tau,tau_result] = MC2(n,K,d,A,mu_hat,Sigma_hat,...
                        tau_hat,p_tau_hat,rho,tau_star,nu_star,...
                        Sigma_star,theta,NBurnIn,NConverge,c,model,...
                        homophily,identifiability);
                end
                parsave(savefile,error_rate,tau,tau_result);
            end
        end
    end
end

% delete(gcp('nocreate'))
% matlabpool close
