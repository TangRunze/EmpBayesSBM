% this function generates the adjacency spectral graph embedding
function Xhat = asge(A,d)
    [V, D] = eigs(A,d,'LA');
    Xhat = V*sqrt(D);
end