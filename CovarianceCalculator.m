function Sigma = CovarianceCalculator(nu,rho)

% nu is your nu and rho is your rho, it returns a list of covariance
% matrices

[K, d] = size(nu);
Delta = zeros(size(nu));
for i = 1:K
    Delta = Delta + nu(i,:)'*nu(i,:)*rho(i);
end

Sigma = zeros(d,d,K);
for i = 1:K
    tmp1 = nu(i,:) * nu';
    tmp2 = tmp1 - tmp1.^2;
    B = zeros(size(nu))';
    for j = 1:K
        B = B + nu(j,:)'*nu(j,:)*tmp2(j)*rho(j);
    end
    Sigma(:,:,i) = Delta\B/Delta;
end
