function ac = CheckS(nu,homophily,identifiability)
% Check if nu*nu' satisfies the probability constraints
% 0 <= <nu_i, nu_j> <= 1
% If homophily = 1, then it also enforces homophily
% <nu_i,nu_j> <= <nu_i,nu_i>
% If identifiability = 1, then it also enforces
% <nu_i,nu_i> >= <nu_j,nu_j> for i > j

K = size(nu,1);
B = nu*nu';
ac = 1;
for i = 1:K
    for j = 1:K
        if (B(i,j) < 0) || (B(i,j) > 1) || ...
                (homophily == 1)&&(B(i,i) < B(i,j)) || ...
                (identifiability == 1)&&(i > j)&&(B(i,i) < B(j,j))
            ac = 0;
            return;
        end
    end
end