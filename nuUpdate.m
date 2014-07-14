function [nu_tilde,f_num] = nuUpdate(tau,nu,mu_hat,Sigma_hat,nu_star,...
    Sigma_star,A,f_den,K,d,c,model,homophily,identifiability)

%% --- Generate valid nu_tilde ---
nu_tilde = nuGenerator(K,d,nu_star,Sigma_star,mu_hat,Sigma_hat,...
    model,c,homophily,identifiability,0);

%% -- pre-calculation of f = <nu_i,nu_j>^p*(1-<nu_i,nu_j>)^(1-p) (KxKx2) --
f_num = fCalculator(nu_tilde);

%% --- Calculate the probability of accepting ---
tmp_num = (1 - A).*log(f_num(tau,tau,1)) + A.*log(f_num(tau,tau,2));
tmp_den = (1 - A).*log(f_den(tau,tau,1)) + A.*log(f_den(tau,tau,2));
tmp = tmp_num - tmp_den;
tmp = triu(tmp,1);
logratio = sum(sum(tmp));

if model == 4
    for i = 1:K
        logratio = logratio + ( - 1/2*(nu(i,:) - mu_hat(i,:))*...
            ((Sigma_hat(:,:,i))\(nu(i,:) - mu_hat(i,:))'));
        logratio = logratio - ( - 1/2*(nu_tilde(i,:) - mu_hat(i,:))*...
            ((Sigma_hat(:,:,i))\(nu_tilde(i,:) - mu_hat(i,:))'));
    end
end

% if c == 1
%     if model == 1
%         for i = 1:K
%             logratio = logratio + ( - 1/2*(nu(i,:) - nu_hat(i,:))*((Sigma_hat(:,:,i))\(nu(i,:) - nu_hat(i,:))'));
%             logratio = logratio - ( - 1/2*(nu_tilde(i,:) - nu_hat(i,:))*((Sigma_hat(:,:,i))\(nu_tilde(i,:) - nu_hat(i,:))'));
%         end
%     else
%         for i = 1:K
%             logratio = logratio + ( - 1/2*(nu(i,:) - nu_hat(i,:))*((Sigma_hat(:,:,i))\(nu(i,:) - nu_hat(i,:))'));
%             logratio = logratio - ( - 1/2*(nu_tilde(i,:) - nu_hat(i,:))*((Sigma_hat(:,:,i))\(nu_tilde(i,:) - nu_hat(i,:))'));
%         end
%     end
% end

% It should be minus nu_tilde, plus nu.

% if c == 2
%     if model == 2
%         for i = 1:K
%             logratio = logratio + ( - 1/2*(nu_tilde(i,:) - nu_hat(i,:))*((Sigma_hat(:,:,i))\(nu_tilde(i,:) - nu_hat(i,:))'));
%             logratio = logratio - ( - 1/2*(nu(i,:) - nu_hat(i,:))*((Sigma_hat(:,:,i))\(nu(i,:) - nu_hat(i,:))'));
%         end
%     elseif model == 3
%         for i = 1:K
%             logratio = logratio + ( - 1/2*(nu_tilde(i,:) - nu_hat(i,:))*((Sigma_star(:,:,i)/n)\(nu_tilde(i,:) - nu_hat(i,:))'));
%             logratio = logratio - ( - 1/2*(nu(i,:) - nu_hat(i,:))*((Sigma_star(:,:,i)/n)\(nu(i,:) - nu_hat(i,:))'));
%         end
%     elseif model == 1
%         for i = 1:K
%             logratio = logratio + ( - 1/2*(nu_tilde(i,:) - nu_star(i,:))*((Sigma_star(:,:,i)/n)\(nu_tilde(i,:) - nu_star(i,:))'));
%             logratio = logratio - ( -1/2*(nu(i,:) - nu_star(i,:))*((Sigma_star(:,:,i)/n)\(nu(i,:) - nu_star(i,:))'));
%         end
%     end
% elseif c == 3
%     tmp = reshape(repmat(T,[d*d,1]),[d,d,K]);
%     if model == 2
%         tmp = Sigma_hat./tmp;
%         for i = 1:K
%             logratio = logratio + ( - 1/2*(nu_tilde(i,:) - nu_hat(i,:))*((tmp(:,:,i))\(nu_tilde(i,:) - nu_hat(i,:))'));
%             logratio = logratio - ( - 1/2*(nu(i,:) - nu_hat(i,:))*((tmp(:,:,i))\(nu(i,:) - nu_hat(i,:))'));
%         end
%     elseif model == 3
%         tmp = Sigma_star./tmp/n;
%         for i = 1:K
%             logratio = logratio + ( - 1/2*(nu_tilde(i,:) - nu_hat(i,:))*((tmp(:,:,i))\(nu_tilde(i,:) - nu_hat(i,:))'));
%             logratio = logratio - ( - 1/2*(nu(i,:) - nu_hat(i,:))*((tmp(:,:,i))\(nu(i,:) - nu_hat(i,:))'));
%         end        
%     elseif model == 1
%         tmp = Sigma_star./tmp/n;
%         for i = 1:K
%             logratio = logratio + ( - 1/2*(nu_tilde(i,:) - nu_star(i,:))*((tmp(:,:,i))\(nu_tilde(i,:) - nu_star(i,:))'));
%             logratio = logratio - ( - 1/2*(nu(i,:) - nu_star(i,:))*((tmp(:,:,i))\(nu(i,:) - nu_star(i,:))'));
%         end        
%     end
% elseif c == 4
%     if model == 2
%         for i = 1:K
%             logratio = logratio + ( - 1/2*(nu_tilde(i,:) - nu_hat(i,:))*((Sigma_hat(:,:,i)/n)\(nu_tilde(i,:) - nu_hat(i,:))'));
%             logratio = logratio - ( - 1/2*(nu(i,:) - nu_hat(i,:))*((Sigma_hat(:,:,i)/n)\(nu(i,:) - nu_hat(i,:))'));
%         end
%     elseif model == 3
%         for i = 1:K
%             logratio = logratio + ( - 1/2*(nu_tilde(i,:) - nu_hat(i,:))*((Sigma_star(:,:,i)/(n^2))\(nu_tilde(i,:) - nu_hat(i,:))'));
%             logratio = logratio - ( - 1/2*(nu(i,:) - nu_hat(i,:))*((Sigma_star(:,:,i)/(n^2))\(nu(i,:) - nu_hat(i,:))'));
%         end
%     elseif model == 1
%         for i = 1:K
%             logratio = logratio + ( - 1/2*(nu_tilde(i,:) - nu_star(i,:))*((Sigma_star(:,:,i)/(n^2))\(nu_tilde(i,:) - nu_star(i,:))'));
%             logratio = logratio - ( - 1/2*(nu(i,:) - nu_star(i,:))*((Sigma_star(:,:,i)/(n^2))\(nu(i,:) - nu_star(i,:))'));
%         end
%     end
% end

ratio = exp(logratio);

%% --- Accept or reject the proposed status according to the ratio ---
flag = rand;

if flag>ratio
    nu_tilde = nu;
    f_num = f_den;
end


