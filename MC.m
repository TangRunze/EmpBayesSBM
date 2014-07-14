function [error_rate_map,tau_map,tau_result] = MC(n,K,d,A,mu_hat,...
    Sigma_hat,tau_hat,rho,tau_star,nu_star,Sigma_star,theta,NBurnIn,...
    NConverge,c,model,homophily,identifiability)
% Inference using 1 MCMC chain.

rng shuffle;

%% --- Initialization ---
% Take the last NConverge iterations as the results after the burn-in part
tau_result = zeros(NConverge,n);

% initialization of tau, use the estimated tau_hat from asge directly
tau = tau_hat;

% calculation of T (1-by-K), number of vertices in each block
T = zeros(1,K);
for i = 1:K
    T(i) = sum(tau == i);
end

% calculation of T (1-by-K), number of vertices in each block for the true
% labels
T0 = zeros(1,K);
for i = 1:K
    T0(i) = sum(tau_star == i);
end

% Generate valid initialial nu (K-by-d)
nu = nuGenerator(K,d,nu_star,Sigma_star,mu_hat,Sigma_hat,...
    model,c,homophily,identifiability,1);

% pre-calculation of f = <nu_i,nu_j>^p*(1-<nu_i,nu_j>)^(1-p) (KxKx2)
f = fCalculator(nu);

%% Burn-in Part
iter = 0;
while (iter < NBurnIn)
    iter = iter + 1;
    % --- Gibbs Sampling of tau ---
    [tau,T] = tauUpdate(tau,A,T,theta,rho,f,n,K,model);
    % --- Metropolis-Hasting of nu ---
    if c~=5
        for NMH = 1:10
            [nu,f] = nuUpdate(tau,nu,mu_hat,Sigma_hat,nu_star,...
                Sigma_star,A,f,K,d,c,model,homophily,identifiability);
        end
    end
end

%% Ergodic Average
iter = 0;
while (iter < NConverge)
    iter = iter + 1;
    % --- Gibbs Sampling of tau ---
    [tau,T] = tauUpdate(tau,A,T,theta,rho,f,n,K,model);
    % --- Metropolis-Hasting of nu ---
    if c~=5
        for NMH = 1:10
            [nu,f] = nuUpdate(tau,nu,mu_hat,Sigma_hat,nu_star,...
                Sigma_star,A,f,K,d,c,model,homophily,identifiability);
        end
    end
    tau_result(iter,:) = tau;
end

tau_map = mode(tau_result);

error_rate_map = n;
permutation = perms(1:K);
for i = 1:factorial(K)
    pos = permutation(i,:);
    tau_tmp = tau_map;
    for j = 1:K
        nv = (tau_map == pos(j));
        tau_tmp(nv) = j;
    end
    if sum(tau_star~=tau_tmp) < error_rate_map
        error_rate_map = sum(tau_star~=tau_tmp);
    end
end

error_rate_map = error_rate_map/n;


