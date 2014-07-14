function [tau,T] = tauUpdate(tau,A,T,theta,rho,f,n,K,model)
% Update block memberships tau.

% version1
% prob = ones(1,K);
% 
% tau_new = zeros(1,n);
% for i = 1:n
%     prob = ones(1,K);
%     for taui = 1:K
%         for j = 1:n
%             if (i~=j)
%                 prob(taui) = prob(taui) * f(taui,tau(j),A(i,j)+1);
%             end
%         end
%         prob = prob/sum(prob);
%         prob(prob<0) = 0;
%     end
%     tau_new(i) = randsample(K,1,true,prob)';
%     T(tau_new(i)) = T(tau_new(i)) + 1;
%     T(tau(i)) = T(tau(i)) - 1;
% end

% version2
% tau_new = zeros(1,n);
% 
% prob1 = zeros(K,n);
% 
% for i = 1:n
%     prob = zeros(K,1);
%     
%     for taui = 1:K
%         tmp = (1 - A(i,:)) .* f(taui,tau(1:n),1) + A(i,:) .* f(taui,tau(1:n),2);
%         tmp(i) = [];
%         prob(taui) = prod(tmp) * Pi_star(taui);
%     end
%     
%     % make sure all greater than zero and normalize
%     prob(prob<0) = 0;
%     prob = prob/sum(prob);
% 
%     prob1(:,i) = prob;
%     
%     tau_new(i) = randsample(K,1,true,prob)';
%     T(tau_new(i)) = T(tau_new(i)) + 1;
%     T(tau(i)) = T(tau(i)) - 1;
% end

% version3

% tau_new = zeros(1,n);
% prob2 = zeros(K,n);

if model == 1
    for i = 1:n
        tmp = (1 - A(i*ones(1,K),:)) .* f(:,tau(1:n),1) + A(i*ones(1,K),:) .* f(:,tau(1:n),2);
        tmp(:,i) = [];
        prob = prod(tmp,2) .* rho(1:K)';
        % make sure all greater than zero and normalize
        prob(prob<0) = 0;
        prob = prob/sum(prob);
        
        %     prob2(:,i) = prob;
        %
        %     tau_new(i) = randsample(K,1,true,prob);
        %     T(tau_new(i)) = T(tau_new(i)) + 1;
        %     T(tau(i)) = T(tau(i)) - 1;
        
        T(tau(i)) = T(tau(i)) - 1;
        
        ind_tmp = 1;
        r = rand;
        while r > prob(ind_tmp)
            ind_tmp = ind_tmp + 1;
            prob(ind_tmp) = prob(ind_tmp) + prob(ind_tmp-1);
        end
        tau(i) = ind_tmp;
%         if rand<prob(1)
%             tau(i) = 1;
%         else
%             tau(i) = 2;
%         end
        
        %     tau(i) = randsample(K,1,true,prob);
        T(tau(i)) = T(tau(i)) + 1;
    end
else
    for i = 1:n
        T(tau(i)) = T(tau(i)) - 1;
        
        f_gamma = T + theta;
        
%         tmp1 = repmat(theta+T,[2,1]) + eye(2);
%         tmp2 = prod(gamma(tmp1),2);
%         tmp2 = tmp2/sum(tmp2);
        
        
        tmp = (1 - A(i*ones(1,K),:)) .* f(:,tau(1:n),1) + A(i*ones(1,K),:) .* f(:,tau(1:n),2);
        tmp(:,i) = [];
        prob = prod(tmp,2) .* f_gamma';
        % make sure all greater than zero and normalize
        prob(prob<0) = 0;
        prob = prob/sum(prob);
        
        %     prob2(:,i) = prob;
        %
        %     tau_new(i) = randsample(K,1,true,prob);
        %     T(tau_new(i)) = T(tau_new(i)) + 1;
        %     T(tau(i)) = T(tau(i)) - 1;
        
        ind_tmp = 1;
        r = rand;
        while r > prob(ind_tmp)
            ind_tmp = ind_tmp + 1;
            prob(ind_tmp) = prob(ind_tmp) + prob(ind_tmp-1);
        end
        tau(i) = ind_tmp;
        
%         if rand<prob(1)
%             tau(i) = 1;
%         else
%             tau(i) = 2;
%         end
        
        %     tau(i) = randsample(K,1,true,prob);
        T(tau(i)) = T(tau(i)) + 1;
    end
end

%% Considering the factor based on A

% tau_new = zeros(1,n);
% 
% tmp_A = zeros(1,n,n);
% tmp_A(1,:,:) = A';
% tmp_A = repmat(tmp_A,[K,1,1]);
% 
% tmp = (1 - tmp_A) .* repmat(f(:,tau(1:n),1),[1,1,n]) + tmp_A .* repmat(f(:,tau(1:n),2),[1,1,n]);
% 
% nv = (zeros(1,n,n) == 1);
% nv(1,:,:) = (eye(n) == 1);
% nv = repmat(nv,[K,1,1]);
% 
% tmp(nv) = [];
% tmp = reshape(tmp,[K,n-1,n]);
% 
% prob = squeeze(prod(tmp,2));
% 
% if (model == 1)
%     % factor of Pi_star
%     prob = prob .* repmat(Pi_star',[1,n]);
%     
%     % make sure all greater than zero
%     prob(prob<0) = 0;
%     % normalization
%     prob = prob./[sum(prob);sum(prob)];
% else
%     %% consider the factor gamma (K-by-K)
%     % f_gamma(i,j) indicates the factor when tau = i changes to tau_new = j
%     % we will devide all the prob by prod(gamma(theta + T))
%     f_gamma = ones(K,K);
%     for i = 1:K
%         for j = 1:K
%             if i ~= j
%                 f_gamma(i,j) = (gamma(theta(i) + T(i))/gamma(theta(i) + T(i) - 1))*(gamma(theta(j) + T(j))/gamma(theta(j) + T(j) + 1));
%             end
%         end
%     end
%     
%     tmp = zeros(K,n);
%     for i = 1:K
%         nv = (tau == i);
%         tmp(:,nv) = repmat(f_gamma(i,:)',[1,sum(nv)]);
%     end
%     
%     prob = prob.*tmp;
%     % make sure all greater than zero
%     prob(prob<0) = 0;
%     % normalization
%     prob = prob./[sum(prob);sum(prob)];
% end
% 
% %% Sampling tau_new from the multinomial distribution
% % r = rand(1,n);
% % nv = (zeros(1,n) == 0);
% % for i = 1:K
% %     if i>1
% %         prob(i,:) = prob(i,:) + prob(i-1,:);
% %     end
% %     nv1 = (r <= prob(i,:)) & nv;
% %     tau_new(nv1) = i;
% %     T(i) = sum(nv1);
% %     nv(nv1) = logical(zeros(1,T(i)));
% % end
% 
% r = rand(1,n);
% nv1 = (r<=prob(1,:));
% nv2 = ~nv1;
% tau_new(nv1) = 1;
% tau_new(nv2) = 2;
% T(1) = sum(nv1);
% T(2) = n - T(1);
