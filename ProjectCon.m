function [c,ceq] = ProjectCon(x,K,d,homophily,identifiability)

%% Pre-calculation
nu_tmp = zeros(K,d);
for i = 1:K
    for j = 1:d
        nu_tmp(i,j) = x((i-1)*d+j);
    end
end
B = nu_tmp*nu_tmp';

%% Inequalities
c = [];
for i = 1:K
    for j = 1:K
        c = [c; -B(i,j); B(i,j)-1];
        if (homophily == 1)
            c = [c; -B(i,i)+B(i,j)];
        end
        if (identifiability == 1) && (i > j)
            c = [c; -B(i,i)+B(j,j)];
        end
    end
end

%% Equalities
ceq = 0;

end