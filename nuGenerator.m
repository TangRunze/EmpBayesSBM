function nu = nuGenerator(K,d,nu_star,Sigma_star,mu_hat,Sigma_hat,...
    model,c,homophily,identifiability,init)

% Generate valid nu based on constraints in S. Init=1 means we are
% generating the initial nu, so c=1 leads to a normal (ASGE) initial
% instead of Flat.

flag = 0;
while (~flag) || ((c~=5) && (~CheckS(nu,homophily,identifiability)))
    flag = 1;
    if c == 1
        if init == 1
            % The initialization is still normal (ASGE) even if the
            % proposal & prior is Flat.
            nu = mvnrnd(mu_hat, Sigma_hat);
            while any(any((nu*nu' > ones(K,K))|(nu*nu' < zeros(K,K))))
                nu = mvnrnd(mu_hat, Sigma_hat);
            end
        else
            % When nu is uniformly distributed, it does not matter what is
            % the mean of the distribution. So all the methods will have
            % the same uniform distribution for the initial point.
            % Since the probability constraints in S (probability should
            % be between 0 & 1), the norm of each nu should be less or
            % equal to 1, i.e. each element of nu should be located inside
            % the unit circle.
            nu = 2*rand(K,d) - 1;
            while any(any((nu*nu' > ones(K,K))|(nu*nu' < zeros(K,K))))
                nu = 2*rand(K,d) - 1;
            end
        end
    elseif c == 2
        if (model == 2) || (model == 4)
            nu = mvnrnd(mu_hat, Sigma_hat);
            while any(any((nu*nu' > ones(K,K))|(nu*nu' < zeros(K,K))))
                nu = mvnrnd(mu_hat, Sigma_hat);
            end
        elseif model == 3
            nu = mvnrnd(mu_hat, Sigma_star/n);
            while any(any((nu*nu' > ones(K,K))|(nu*nu' < zeros(K,K))))
                nu = mvnrnd(mu_hat, Sigma_star/n);
            end
        elseif model == 1
            nu = mvnrnd(nu_star, Sigma_star/n);
            while any(any((nu*nu' > ones(K,K))|(nu*nu' < zeros(K,K))))
                nu = mvnrnd(nu_star, Sigma_star/n);
            end
        end
    elseif c == 3
        tmp = reshape(repmat(T0,[d*d,1]),[d,d,K]);
        if model == 2
            nu = mvnrnd(mu_hat, Sigma_hat./tmp);
            while any(any((nu*nu' > ones(K,K))|(nu*nu' < zeros(K,K))))
                nu = mvnrnd(mu_hat, Sigma_hat./tmp);
            end
        elseif model == 3
            nu = mvnrnd(mu_hat, Sigma_star./tmp/n);
            while any(any((nu*nu' > ones(K,K))|(nu*nu' < zeros(K,K))))
                nu = mvnrnd(mu_hat, Sigma_star./tmp/n);
            end
        elseif model == 1
            nu = mvnrnd(nu_star, Sigma_star./tmp/n);
            while any(any((nu*nu' > ones(K,K))|(nu*nu' < zeros(K,K))))
                nu = mvnrnd(nu_star, Sigma_star./tmp/n);
            end
        end
    elseif c == 4
        if model == 2
            nu = mvnrnd(mu_hat, Sigma_hat/n);
            while any(any((nu*nu' > ones(K,K))|(nu*nu' < zeros(K,K))))
                nu = mvnrnd(mu_hat, Sigma_hat/n);
            end
        elseif model == 3
            nu = mvnrnd(mu_hat, Sigma_star/(n^2));
            while any(any((nu*nu' > ones(K,K))|(nu*nu' < zeros(K,K))))
                nu = mvnrnd(mu_hat, Sigma_star/(n^2));
            end
        elseif model == 1
            nu = mvnrnd(nu_star, Sigma_star/(n^2));
            while any(any((nu*nu' > ones(K,K))|(nu*nu' < zeros(K,K))))
                nu = mvnrnd(nu_star, Sigma_star/(n^2));
            end
        end
    elseif c == 5
        if (model == 2) || (model ==3)
            nu = mu_hat;
        elseif model == 1
            nu = nu_star;
        end
    end
end
