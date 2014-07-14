function f = fCalculator(nu)
% pre-calculation of f = <nu_i,nu_j>^p*(1-<nu_i,nu_j>)^(1-p) (KxKx2)

K = size(nu,1);
f = zeros(K,K,2);
for i = 1:K
    for j = i:K
        for p = 1:2
            f(i,j,p) = (nu(i,:)*nu(j,:)')^(p-1)*(1-nu(i,:)*nu(j,:)')^(2-p);
            f(j,i,p) = f(i,j,p);
        end
    end
end