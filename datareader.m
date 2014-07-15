function [nVertex, adjMatrix, muHat, sigmaHat, tauHat, pTauHat, tauStar]...
    = datareader(nBlock, dimLatentPosition, iGraph, hasLabel)
% Read and preprocess data

tauStar = 0;

if (exist(['data/graph' int2str(iGraph) '.mat'], 'file') == 0)
    if ((exist(['data/label' int2str(iGraph)], 'file') == 0) && ...
        (hasLabel == 1))
        error(['label' int2str(iGraph) 'does not exist!'])
    elseif (exist(['data/graph' int2str(iGraph)], 'file') == 0)
        error(['graph' int2str(iGraph) 'does not exist!'])
    else
        % Part 1: Read graph adjacency matrix.
        adjMatrix = textread(['data/graph' int2str(iGraph)]);
        nVertex = size(adjMatrix, 1);
        if (hasLabel == 1)
            tauStar = textread(['data/label' int2str(iGraph)]);
            if (length(tauStar) ~= nVertex)
                error(['label' int2str(iGraph) 'has wrong number!'])
            end
        end
        % Part 2: Obtain estimates from ASGE o GMM.
        xHat = asge(adjMatrix, dimLatentPosition);
        gm = fitgmdist(xHat, nBlock, 'Replicates', 10);
        tauHat = cluster(gm, xHat)';
        pTauHat = posterior(gm, xHat)';
        muHat = gm.mu;
        sigmaHat = gm.Sigma;
        
        % Save the data
        save(['data/graph' int2str(iGraph)], 'adjMatrix', 'tauHat', ...
            'pTauHat', 'muHat', 'sigmaHat', 'nVertex', 'tauStar');
    end
else
    % Read the existing data
    data = load(['data/graph' int2str(iGraph) '.mat']);
    adjMatrix = data.adjMatrix;
    muHat = data.muHat;
    sigmaHat = data.sigmaHat;
    tauHat = data.tauHat;
    pTauHat = data.pTauHat;
    nVertex = data.nVertex;
    tauStar = data.tauStar;
end
