function [adjMatrix, muHat, sigmaHat, tauHat, pTauHat] = ...
    datagenerator(nVertex, nBlock, dimLatentPosition, B, rho, ...
    epsilonInB, iGraph)
% Generate data if there does not exist one, otherwise read the
% existing data.

% Calculate the sizes
nVectorStar = nVertex*rho;
nVectorStarStart = cumsum(nVectorStar);
nVectorStarStart = [1, nVectorStarStart(1:(end-1)) + 1];
nVectorStarEnd = cumsum(nVectorStar);

if exist(['data/sim-n' num2str(nVertex) '-eps' num2str(epsilonInB) ...
        '-graph' int2str(iGraph) '.mat'], 'file') == 0
    
    disp(['Generating graph ' int2str(iGraph) '...'])
    
    % Part 1: Generate graph adjacency matrix.
    adjMatrix = zeros(nVertex);
    for iBlock = 1:nBlock
        for jBlock = iBlock:nBlock
            adjMatrix(nVectorStarStart(iBlock):nVectorStarEnd(iBlock), ...
                nVectorStarStart(jBlock):nVectorStarEnd(jBlock)) = ...
                binornd(1, B(iBlock,jBlock), nVectorStar(iBlock), ...
                nVectorStar(jBlock));
        end
    end
    adjMatrix = triu(adjMatrix, 1);
    adjMatrix = adjMatrix + adjMatrix';
    
    % Part 2: Obtain estimates from ASGE o GMM.
    xHat = asge(adjMatrix, dimLatentPosition);
    
    gm = fitgmdist(xHat, nBlock, 'Replicates', 10);
    
    tauHat = cluster(gm, xHat)';
    % pihat = gm.PComponents;
    pTauHat = posterior(gm, xHat)';
    muHat = gm.mu;
    sigmaHat = gm.Sigma;
    
    % Plot
    % cl_nv = false(K,n);
    % for i = 1:K
    %     cl_nv(i,:) = (idx == i);
    % end
    % scatter(Xhat(cl_nv(1,:),1),Xhat(cl_nv(1,:),2),10,'r+');
    % hold on
    % scatter(Xhat(cl_nv(2,:),1),Xhat(cl_nv(2,:),2),10,'bo');
    % hold on
    % scatter(Xhat(cl_nv(3,:),1),Xhat(cl_nv(3,:),2),10,'g.');
    % hold off
    % legend('Cluster 1','Cluster 2','Cluster 3','Location','NW')
    
    % Save the data
    save(['data/sim-n' num2str(nVertex) '-eps' num2str(epsilonInB) ...
        '-graph' int2str(iGraph) '.mat'], 'adjMatrix', 'tauHat', ...
        'pTauHat', 'muHat', 'sigmaHat');
else
    % Read the existing data
    data = load(['data/sim-n' num2str(nVertex) '-eps' ...
        num2str(epsilonInB) '-graph' int2str(iGraph) '.mat']);
    adjMatrix = data.adjMatrix;
    muHat = data.muHat;
    sigmaHat = data.sigmaHat;
    tauHat = data.tauHat;
    pTauHat = data.pTauHat;
end
