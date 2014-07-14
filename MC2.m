function [error_rate_map,tau_map,tau_result] = ...
    MC2(n,K,d,A,mu_hat,Sigma_hat,tau_hat,p_tau_hat,rho,tau_star,nu_star,...
    Sigma_star,theta,NBurnIn,NConverge,c,model,homophily,identifiability)
% Inference using 2 MCMC chains.

rng shuffle;

%% Metropolis-Hasting within Gibbs Sampling

Converge = 0;
tau_tmp1 = zeros(NConverge,n);
tau_tmp2 = zeros(NConverge,n);
Q = zeros(NConverge,2);
nboot = NConverge;

iter_max = NBurnIn + NConverge;
iter_out = 0;
restart = 0;

while (~Converge) || (iter_out <= NConverge)
    % If it is the beginning of the chain, initialize the parameters. Chain
    % 1 is using a fixed initialization, so when restart=1, it is better to
    % use what we have as initial.
    if (iter_out == 0)
        % Initialization of chain 2 based on the EM posterior prob
        tmp = mnrnd(1,p_tau_hat');
        tau2 = K + 1 - sum(cumsum(tmp'));
        
        % calculation of T2 (1-by-K), number of vertices in each block for
        % the 2nd chain
        T2 = zeros(1,K);
        for i = 1:K
            T2(i) = sum(tau2 == i);
        end
        
        % Generate valid initialial nu2 (K-by-d) for chain 2
        nu2 = nuGenerator(K,d,nu_star,Sigma_star,mu_hat,Sigma_hat,...
            model,c,homophily,identifiability,1);
        
        % pre-calculation of
        % f2 = <nu2_i,nu2_j>^p*(1-<nu2_i,nu2_j>)^(1-p) (KxKx2)
        % for chain 2
        f2 = fCalculator(nu2);
        
        if restart == 0
            % Initialization of chain 1 based on asge estimator
            tau1 = tau_hat;
            
            % calculation of T1 (1-by-K), number of vertices in each block
            % for the 1st chain
            T1 = zeros(1,K);
            for i = 1:K
                T1(i) = sum(tau1 == i);
            end
            
            % Generate valid initialial nu1 (K-by-d) for chain 1
            nu1 = nuGenerator(K,d,nu_star,Sigma_star,mu_hat,Sigma_hat,...
                model,c,homophily,identifiability,1);
            
            % pre-calculation of
            % f1 = <nu1_i,nu1_j>^p*(1-<nu1_i,nu1_j>)^(1-p) (KxKx2)
            % for chain 1
            f1 = fCalculator(nu1);
        end
    end

    iter_in = 0;
    while (iter_in < NConverge)
        iter_in = iter_in + 1;
        % --- Gibbs Sampling of tau ---
        [tau1,T1] = tauUpdate(tau1,A,T1,theta,rho,f1,n,K,model);
        [tau2,T2] = tauUpdate(tau2,A,T2,theta,rho,f2,n,K,model);
        % --- Metropolis-Hasting of nu ---
        if c~=5
            for NMH = 1:10
                [nu1,f1] = nuUpdate(tau1,nu1,mu_hat,Sigma_hat,nu_star,...
                    Sigma_star,A,f1,K,d,c,model,homophily,identifiability);
                [nu2,f2] = nuUpdate(tau2,nu2,mu_hat,Sigma_hat,nu_star,...
                    Sigma_star,A,f2,K,d,c,model,homophily,identifiability);
            end
        end
        tau_tmp1(iter_in,:) = tau1;
        tau_tmp2(iter_in,:) = tau2;
    end
    iter_out = iter_out + NConverge;
    
    % --- GR test ---
    Q(:,1) = mean(repmat(tau_star,NConverge,1) ~= tau_tmp1,2);
    Q(:,2) = mean(repmat(tau_star,NConverge,1) ~= tau_tmp2,2);
    % Potential Scale Reduction Factor
    % Rhat = psrf(Q(:,1),Q(:,2));
    bootRhat = bootstrap(Q(:,1),Q(:,2),nboot);
    bootRhat = sort(bootRhat);
    if bootRhat(0.95*nboot) < 1.1
        Converge = 1;
    elseif (bootRhat(0.95*nboot) >= 1.1) && (iter_out >= iter_max)
        iter_out = 0;
        restart = 1;
    end
end

%% Calculate Result

tau_result = tau_tmp1;

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


