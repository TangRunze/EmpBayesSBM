function [muHat, sigmaHat, tauHat, pTauHat] = ChangeOrder(muHat, ...
    sigmaHat, tauHat, pTauHat, nBlock)
% Reorder the blocks when they do not satisfy the identifiability
% constraints.

tmpVariable = diag(muHat*muHat');
[~, position] = sort(tmpVariable);
tmpMuHat = muHat(position, :);
tmpSigmaHat = sigmaHat(:, :, position);
tmpPTauHat = pTauHat(position, :);
tmpTauHat = tauHat;
for iBlock = 1:nBlock
    nv = (tauHat == position(iBlock));
    tmpTauHat(nv) = iBlock;
end
muHat = tmpMuHat;
sigmaHat = tmpSigmaHat;
tauHat = tmpTauHat;
pTauHat = tmpPTauHat;