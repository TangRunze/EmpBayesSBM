% 
% nu_hat = [-0.7174, 0.0002; -0.6054, -0.0052];
% 
% % Generate valid initialial nu (K-by-d)
% flag = 0;
% nu = zeros(K,d);
% while (~flag) || (~Check_S(nu,K))
%     flag = 1;
%     nu = mvnrnd(nu_hat, Sigma_hat/n);
%     
%     while any(any((nu*nu' > ones(2,2))|(nu*nu' < zeros(2,2))))
%         nu = mvnrnd(nu_hat, Sigma_hat/n);
%     end
%     
%     % For the restriction nu1*nu1' < nu2*nu2' in S
%     if nu(1,:)*nu(1,:)' >= nu(2,:)*nu(2,:)'
%         tmp_nu = nu;
%         nu(1,:) = tmp_nu(2,:);
%         nu(2,:) = tmp_nu(1,:);
%     end
% end
% 
% beta = pi - atan(nu_hat(1,2)/nu_hat(1,1));
% [cos(beta),sin(beta);-sin(beta),cos(beta)]*nu_hat(1,:)'
% [cos(beta),sin(beta);-sin(beta),cos(beta)]*nu_hat(2,:)'

nu_hat = [-0.7174, 0.0002; -0.6054, -0.0052];
tmp = nu_hat';
tmp = tmp(:);
nu = fmincon(@(x) ProjectObj(x,tmp),tmp,[],[],[],[],-ones(4,1),ones(4,1),@ProjectCon);
