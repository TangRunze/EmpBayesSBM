
maxiter = 100000000;

for iter = 1:maxiter
    T2 = zeros(1,K);
    for i = 1:K
        T2(i) = sum(tau2 == i);
    end
    
    T1 = tabulate(tau2);
    T1 = T1(:,2)';
end
