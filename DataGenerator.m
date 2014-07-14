function [A,mu_hat,Sigma_hat,tau_hat,p_tau_hat] = ...
    DataGenerator(n,K,d,B,rho,g)

ni = n*rho;
ni_start = cumsum(ni);
ni_start = [1, ni_start(1:(end-1)) + 1];
ni_end = cumsum(ni);


if exist(['data/graph' int2str(g)],'file') == 0
    
    disp(['Generating graph ' int2str(g) '...'])
    
    % Part 1: Generate graph adjacency matrix.
    A = zeros(n);
    for i = 1:K
        for j = i:K
            A(ni_start(i):ni_end(i),ni_start(j):ni_end(j)) = ...
                binornd(1,B(i,j),ni(i),ni(j));
        end
    end
    A = triu(A,1);
    A = A + A';
    
    % Part 2: Obtain estimates from ASGE o GMM.
    Xhat = asge(A,d);
    
    gm = fitgmdist(Xhat,K,'Replicates',10);
    
    tau_hat = cluster(gm,Xhat)';
    % pihat = gm.PComponents;
    p_tau_hat = posterior(gm,Xhat)';
    mu_hat = gm.mu;
    Sigma_hat = gm.Sigma;
    
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
    save(['data/graph' int2str(g)],'A','tau_hat','p_tau_hat',...
        'mu_hat','Sigma_hat');
else
    data = load(['data/graph' int2str(iter) '.mat']);
    A = data.A;
end
