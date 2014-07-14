function f = fCalculator(nu)
% pre-calculation of f = <nu_i,nu_j>^p*(1-<nu_i,nu_j>)^(1-p) (KxKx2)

nBlock = size(nu, 1);
f = zeros(nBlock, nBlock, 2);
for iBlock = 1:nBlock
    for jBlock = iBlock:nBlock
        for p = 1:2
            f(iBlock, jBlock, p) = (nu(iBlock, :)*nu(jBlock, :)')^ ...
                (p - 1)*(1 - nu(iBlock, :)*nu(jBlock, :)')^(2 - p);
            f(jBlock, iBlock, p) = f(iBlock, jBlock, p);
        end
    end
end