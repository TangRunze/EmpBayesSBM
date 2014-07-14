function [nuTilde, fNumerator] = nuUpdate(nVertex, tau, nu, muHat, ...
    sigmaHat, nuStar, sigmaStar, adjMatrix, fDenominator, nBlock, ...
    dimLatentPosition, scaleCovariance, modelType, isHomophily, ...
    isIdentifiable)
% Update latent positions nu

%% --- Generate valid nu_tilde ---
nuTilde = nuGenerator(nVertex, nBlock, dimLatentPosition, nuStar, ...
    sigmaStar, muHat, sigmaHat, scaleCovariance, modelType, isHomophily,...
    isIdentifiable, 0);

%% pre-calculation of f = <nu_i,nu_j>^p*(1-<nu_i,nu_j>)^(1-p)
% (nBlock x nBlock x 2)
fNumerator = fCalculator(nuTilde);

%% --- Calculate the probability of accepting ---
tmpMatrixNumerator = (1 - adjMatrix).*log(fNumerator(tau, tau, 1)) + ...
    adjMatrix.*log(fNumerator(tau, tau, 2));
tmpMatrixDenominator = (1 - adjMatrix).*log(fDenominator(tau, tau, 1)) +...
    adjMatrix.*log(fDenominator(tau, tau, 2));
tmpMatrix = tmpMatrixNumerator - tmpMatrixDenominator;
tmpMatrix = triu(tmpMatrix, 1);
logRatio = sum(sum(tmpMatrix));

if modelType == 4
    for iBlock = 1:nBlock
        logRatio = logRatio + ( - 1/2*(nu(iBlock, :) - muHat(iBlock, :))...
            *((sigmaHat(:, :, iBlock))\(nu(iBlock, :) - ...
            muHat(iBlock, :))'));
        logRatio = logRatio - ( - 1/2*(nuTilde(iBlock, :) - ...
            muHat(iBlock, :))*((sigmaHat(:, :, iBlock))...
            \(nuTilde(iBlock, :) - muHat(iBlock, :))'));
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

ratio = exp(logRatio);

%% --- Accept or reject the proposed status according to the ratio ---
isAccepted = rand;

if (isAccepted > ratio)
    nuTilde = nu;
    fNumerator = fDenominator;
end


