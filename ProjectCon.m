function [c, ceq] = ProjectCon(x, nBlock, dimLatentPosition, ...
    isHomophily, isIdentifiable)
% Condition functions in projection optimization problem

%% Pre-calculation
tmpNu = zeros(nBlock, dimLatentPosition);
for iBlock = 1:nBlock
    for jBlock = 1:dimLatentPosition
        tmpNu(iBlock, jBlock) = x((iBlock - 1)*dimLatentPosition + jBlock);
    end
end
B = tmpNu*tmpNu';

%% Inequalities
c = [];
for iBlock = 1:nBlock
    for jBlock = 1:nBlock
        c = [c; - B(iBlock, jBlock); B(iBlock, jBlock) - 1];
        if (isHomophily == 1)
            c = [c; - B(iBlock, iBlock) + B(iBlock, jBlock)];
        end
        if (isIdentifiable == 1) && (iBlock > jBlock)
            c = [c; - B(iBlock, iBlock) + B(jBlock, jBlock)];
        end
    end
end

%% Equalities
ceq = 0;

end