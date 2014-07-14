function sigma = covariancecalculator(nu, rho)
% Input nu and rho, where nu is the latent positions and rho is block
% proportion, CovarianceCalculator returns a list of covariance matrices.

[nBlock, dimLatentPosition] = size(nu);
Delta = zeros(dimLatentPosition, dimLatentPosition);
for iBlock = 1:nBlock
    Delta = Delta + nu(iBlock,:)'*nu(iBlock,:)*rho(iBlock);
end

sigma = zeros(dimLatentPosition, dimLatentPosition, nBlock);
for iBlock = 1:nBlock
    tmpVariable1 = nu(iBlock, :)*nu';
    tmpVariable2 = tmpVariable1 - tmpVariable1.^2;
    B = zeros(dimLatentPosition, dimLatentPosition)';
    for jBlock = 1:nBlock
        B = B + nu(jBlock, :)'*nu(jBlock, :)*tmpVariable2(jBlock)*...
            rho(jBlock);
    end
    sigma(:, :, iBlock) = Delta\B/Delta;
end