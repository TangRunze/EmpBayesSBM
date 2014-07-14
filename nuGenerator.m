function nu = nugenerator(nVertex, nBlock, dimLatentPosition, nuStar, ...
    sigmaStar, muHat, sigmaHat, scaleCovariance, modelType, isHomophily,...
    isIdentifiable, initType)
% Generate valid nu based on constraints in S. initType = 1 means we are
% generating the initial nu, so scaleCovariance = 1 leads to a normal 
% (ASGE) initial instead of Flat.

flag = 0;
while (~flag) || ((scaleCovariance ~= 5) && ...
        (~checkconstraints(nu, isHomophily, isIdentifiable)))
    flag = 1;
    if scaleCovariance == 1
        if initType == 1
            % The initialization is still normal (ASGE) even if the
            % proposal & prior is Flat.
            nu = mvnrnd(muHat, sigmaHat);
            while any(any((nu*nu' > ones(nBlock, nBlock))|...
                    (nu*nu' < zeros(nBlock, nBlock))))
                nu = mvnrnd(muHat, sigmaHat);
            end
        else
            % When nu is uniformly distributed, it does not matter what is
            % the mean of the distribution. So all the methods will have
            % the same uniform distribution for the initial point.
            % Since the probability constraints in S (probability should
            % be between 0 & 1), the norm of each nu should be less or
            % equal to 1, i.e. each element of nu should be located inside
            % the unit circle.
            nu = 2*rand(nBlock, dimLatentPosition) - 1;
            while any(any((nu*nu' > ones(nBlock, nBlock)) | ...
                    (nu*nu' < zeros(nBlock, nBlock))))
                nu = 2*rand(nBlock, dimLatentPosition) - 1;
            end
        end
    elseif scaleCovariance == 2
        if (modelType == 2) || (modelType == 4)
            nu = mvnrnd(muHat, sigmaHat);
            while any(any((nu*nu' > ones(nBlock, nBlock)) | ...
                    (nu*nu' < zeros(nBlock, nBlock))))
                nu = mvnrnd(muHat, sigmaHat);
            end
        elseif modelType == 3
            nu = mvnrnd(muHat, sigmaStar/nVertex);
            while any(any((nu*nu' > ones(nBlock, nBlock)) | ...
                    (nu*nu' < zeros(nBlock, nBlock))))
                nu = mvnrnd(muHat, sigmaStar/nVertex);
            end
        elseif modelType == 1
            nu = mvnrnd(nuStar, sigmaStar/nVertex);
            while any(any((nu*nu' > ones(nBlock, nBlock)) | ...
                    (nu*nu' < zeros(nBlock, nBlock))))
                nu = mvnrnd(nuStar, sigmaStar/nVertex);
            end
        end
    elseif scaleCovariance == 3
        tmpMatrix = reshape(repmat(T0, ...
            [dimLatentPosition*dimLatentPosition, 1]), ...
            [dimLatentPosition, dimLatentPosition,nBlock]);
        if modelType == 2
            nu = mvnrnd(muHat, sigmaHat./tmpMatrix);
            while any(any((nu*nu' > ones(nBlock, nBlock)) | ...
                    (nu*nu' < zeros(nBlock, nBlock))))
                nu = mvnrnd(muHat, sigmaHat./tmpMatrix);
            end
        elseif modelType == 3
            nu = mvnrnd(muHat, sigmaStar./tmpMatrix/nVertex);
            while any(any((nu*nu' > ones(nBlock, nBlock)) | ...
                    (nu*nu' < zeros(nBlock, nBlock))))
                nu = mvnrnd(muHat, sigmaStar./tmpMatrix/nVertex);
            end
        elseif modelType == 1
            nu = mvnrnd(nuStar, sigmaStar./tmpMatrix/nVertex);
            while any(any((nu*nu' > ones(nBlock, nBlock)) | ...
                    (nu*nu' < zeros(nBlock, nBlock))))
                nu = mvnrnd(nuStar, sigmaStar./tmpMatrix/nVertex);
            end
        end
    elseif scaleCovariance == 4
        if modelType == 2
            nu = mvnrnd(muHat, sigmaHat/nVertex);
            while any(any((nu*nu' > ones(nBlock, nBlock)) | ...
                    (nu*nu' < zeros(nBlock, nBlock))))
                nu = mvnrnd(muHat, sigmaHat/nVertex);
            end
        elseif modelType == 3
            nu = mvnrnd(muHat, sigmaStar/(nVertex^2));
            while any(any((nu*nu' > ones(nBlock, nBlock)) | ...
                    (nu*nu' < zeros(nBlock, nBlock))))
                nu = mvnrnd(muHat, sigmaStar/(nVertex^2));
            end
        elseif modelType == 1
            nu = mvnrnd(nuStar, sigmaStar/(nVertex^2));
            while any(any((nu*nu' > ones(nBlock, nBlock)) | ...
                    (nu*nu' < zeros(nBlock, nBlock))))
                nu = mvnrnd(nuStar, sigmaStar/(nVertex^2));
            end
        end
    elseif scaleCovariance == 5
        if (modelType == 2) || (modelType ==3)
            nu = muHat;
        elseif modelType == 1
            nu = nuStar;
        end
    end
end
