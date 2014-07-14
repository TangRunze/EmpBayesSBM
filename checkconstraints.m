function isValid = checkconstraints(nu, requireHomophily, requireIdentifiability)
% Check if nu*nu' satisfies the probability constraints
% 0 <= <nu_i, nu_j> <= 1
% If homophily = 1, then it also enforces homophily
% <nu_i,nu_j> <= <nu_i,nu_i>
% If identifiability = 1, then it also enforces
% <nu_i,nu_i> >= <nu_j,nu_j> for i > j

nBlock = size(nu, 1);
B = nu*nu';
isValid = 1;
for iBlock = 1:nBlock
    for jBlock = 1:nBlock
        isProbability = ((B(iBlock, jBlock) >= 0) && ...
            (B(iBlock, jBlock) <= 1));
        isHomophily = (B(iBlock, iBlock) >= B(iBlock, jBlock));
        isIdentifiable = ((iBlock > jBlock) && ...
            (B(iBlock, iBlock) >= B(jBlock, jBlock)));
        if (~isProbability || ...
                (requireHomophily == 1)&&(~isHomophily) || ...
                (requireIdentifiability == 1)&&(~isIdentifiable))
            isValid = 0;
            return;
        end
    end
end