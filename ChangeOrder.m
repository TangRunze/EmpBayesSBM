function [mu_hat,Sigma_hat,tau_hat,p_tau_hat] = ...
    ChangeOrder(mu_hat,Sigma_hat,tau_hat,p_tau_hat,K)
% Reorder the blocks when they do not satisfy the identifiability
% constraints.

tmp = diag(mu_hat*mu_hat');
[~, pos] = sort(tmp);
mu_hat_tmp = mu_hat(pos,:);
Sigma_hat_tmp = Sigma_hat(:,:,pos);
p_tau_hat_tmp = p_tau_hat(pos,:);
tau_hat_tmp = tau_hat;
for j = 1:K
    nv = (tau_hat == pos(j));
    tau_hat_tmp(nv) = j;
end
mu_hat = mu_hat_tmp;
Sigma_hat = Sigma_hat_tmp;
tau_hat = tau_hat_tmp;
p_tau_hat = p_tau_hat_tmp;