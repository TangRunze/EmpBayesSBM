function Sigma = CovarianceCalculator(nu,rho)
% Input nu and rho, where nu is the latent positions and rho is block
% proportion, CovarianceCalculator returns a list of covariance matrices.

[K, d] = size(nu);
Delta = zeros(d,d);
for i = 1:K
    Delta = Delta + nu(i,:)'*nu(i,:)*rho(i);
end

Sigma = zeros(d,d,K);
for i = 1:K
    tmp1 = nu(i,:) * nu';
    tmp2 = tmp1 - tmp1.^2;
    B = zeros(d,d)';
    for j = 1:K
        B = B + nu(j,:)'*nu(j,:)*tmp2(j)*rho(j);
    end
    Sigma(:,:,i) = Delta\B/Delta;
end
