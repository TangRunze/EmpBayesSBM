EmpBayesSBM
===========

Empirical Bayes Stochastic Block Model

EmpBayesSBM implements an empirical Bayes methodology for estimation of
block memberships of vertices in a random graph drawn from the stochastic
blockmodel.

There are two versions of EmpBayesSBM: simulation (ebsbm.m) & real data 
(ebsbmsim.m).

========== real data version ==========

EBSBM(nBlock, gStart, gEnd, scaleCovarianceStart, scaleCovarianceEnd, ...
    hasLabel)

nBlock selects the number of blocks in the stochastic blockmodel.

gStart/gEnd selects the starting/ending graph in the simulation.

scaleCovarianceStart/scaleCovarianceEnd selects the starting/ending scale
of the covariance matrices.
scaleCovarianceStart/scaleCovarianceEnd should be integers between 1 and 
5.
      scaleCovariance = 1:  0
      scaleCovariance = 2:  n
      scaleCovariance = 3:  n*n_k
      scaleCovariance = 4:  n^2
      scaleCovariance = 5:  Infinity

hasLabel = 1 means true labels are known. Otherwise truelabels are
unknown.
Default hasLabel = 0.

EXAMPLE: 
      EBSBM(3, 1, 10, 4, 5, 1)

EBSBM(nBlock, gStart, gEnd, scaleCovarianceStart, scaleCovarianceEnd, ...
    hasLabel, nCore)

nCore selects the number of cores in parallel programming.
nCore = 1 do NOT run the code parallel.
nCore should be integers between 1 and 12.
Default nCore = 1.

EXAMPLE:
      EBSBM(3, 1, 10, 4, 5, 1, 1)     NOT parallel
      EBSBM(3, 1, 10, 4, 5, 1, 12)    use 12 cores

EBSBM(nBlock, gStart, gEnd, scaleCovarianceStart, scaleCovarianceEnd, ...
    hasLabel, nCore, nBurnIn, nConverge)

nBurnIn selects the number of iterations in burn-in part.
Default nBurnIn = 19000.

nConverge selects the number of iterations for analysis after the burn-in
part.
Default nConverge = 500.

EXAMPLE:
      EBSBM(3, 1, 10, 4, 5, 1, 1, 10000, 1000)

EBSBM(nBlock, gStart, gEnd, scaleCovarianceStart, scaleCovarianceEnd, ...
    hasLabel, nCore, nBurnIn, nConverge, dimLatentPosition)

dimLatentPosition selects the dimension of latent positions.
dimLatentPosition should be integers between 1 and nBlock.
Default dimLatentPosition = nBlock.

EXAMPLE:
      EBSBM(3, 1, 10, 4, 5, 1, 1, 10000, 1000, 2)

EBSBM(nBlock, gStart, gEnd, scaleCovarianceStart, scaleCovarianceEnd, ...
    hasLabel, nCore, nBurnIn, nConverge, dimLatentPosition, ...
    isHomophily, isIdentifiable)

isHomophily controls the constraints. 
If isHomophily = 1, then it enforces homophily:
      <nu_i,nu_j> <= <nu_i,nu_i>  for any i, j
Default isHomophily = 1.

isIdentifiable controls the constraints.
If isIdentifiable = 1, then it enforces identifiability:
      <nu_i,nu_i> >= <nu_j,nu_j>  for any i > j
Default isIdentifiable = 0.

EXAMPLE:
      EBSBM(3, 1, 10, 4, 5, 1, 1, 10000, 1000, 2, 1, 1)

EBSBM(nBlock, gStart, gEnd, scaleCovarianceStart, scaleCovarianceEnd, ...
    hasLabel, nCore, nBurnIn, nConverge, dimLatentPosition, ...
    isHomophily, isIdentifiable, theta)

theta selects the hyperparameters for the prior distribution for rho.
Default theta = [1, 1, ..., 1] of length nBlock.

EXAMPLE:
      EBSBM(3, 1, 10, 4, 5, 1, 1, 10000, 1000, 2, 1, 1, [1, 3, 2])

EBSBM(nBlock, gStart, gEnd, scaleCovarianceStart, scaleCovarianceEnd, ...
    hasLabel, nCore, nBurnIn, nConverge, dimLatentPosition, ...
    isHomophily, isIdentifiable, theta, diagonalAugmentation)

diagonalAugmentation = 1 means using diagonal augmentation when
preprocessing the data.
Default diagonalAugmentation = 1.

EXAMPLE:
      EBSBM(3, 1, 10, 4, 5, 1, 1, 10000, 1000, 2, 1, 1, [1, 3, 2], 0)


========== simulation version ==========

EBSBMSIM(nVertex, nBlock, epsilonInB, gStart, gEnd, ...
    modelTypeStart, modelTypeEnd, scaleCovarianceStart,
    scaleCovarianceEnd)

nVertex selects the number of vertices in the graph.

nBlock selects the number of blocks in the stochastic blockmodel.

epsilonInB controls the true model. The probability matrix
      B = (0.5 - epsilonInB)*J + 2*epsilonInB*I
epsilonInB should be inside [0, 0.5].

gStart/gEnd selects the starting/ending graph in the simulation.

modelTypeStart/modelTypeEnd selects the starting/ending type of model.
modelTypeStart/modelTypeEnd should be integers between 1 and 4.
      modelType = 1:    Gold
      modelType = 2:    ASGE
      modelType = 3:    ASGE1
      modelType = 4:    Flat

scaleCovarianceStart/scaleCovarianceEnd selects the starting/ending scale
of the covariance matrices.
scaleCovarianceStart/scaleCovarianceEnd should be integers between 1 and 
5.
      scaleCovariance = 1:  0
      scaleCovariance = 2:  n
      scaleCovariance = 3:  n*n_k
      scaleCovariance = 4:  n^2
      scaleCovariance = 5:  Infinity

EXAMPLE: 
      EBSBMSIM(150, 3, 0.1, 1, 10, 1, 2, 4, 5)

EBSBMSIM(nVertex, nBlock, epsilonInB, gStart, gEnd, ...
    modelTypeStart, modelTypeEnd, scaleCovarianceStart, ...
    scaleCovarianceEnd, nCore)

nCore selects the number of cores in parallel programming.
nCore = 1 do NOT run the code parallel.
nCore should be integers between 1 and 12.
Default nCore = 1.

EXAMPLE:
      EBSBMSIM(150, 3, 0.1, 1, 10, 1, 2, 4, 5, 1)     NOT parallel
      EBSBMSIM(150, 3, 0.1, 1, 10, 1, 2, 4, 5, 12)    use 12 cores

EBSBMSIM(nVertex, nBlock, epsilonInB, gStart, gEnd, ...
    modelTypeStart, modelTypeEnd, scaleCovarianceStart, ...
    scaleCovarianceEnd, nCore, nBurnIn, nConverge)

nBurnIn selects the number of iterations in burn-in part.
Default nBurnIn = 19000.

nConverge selects the number of iterations for analysis after the burn-in
part.
Default nConverge = 500.

EXAMPLE:
      EBSBMSIM(150, 3, 0.1, 1, 10, 1, 2, 4, 5, 1, 10000, 1000)

EBSBMSIM(nVertex, nBlock, epsilonInB, gStart, gEnd, ...
    modelTypeStart, modelTypeEnd, scaleCovarianceStart, ...
    scaleCovarianceEnd, nCore, nBurnIn, nConverge, rho)

rho selects the block proportions.
rho should sum up to 1.
Default rho = [1/nBlock, 1/nBlock, ..., 1/nBlock] of length nBlock.

EXAMPLE:
      EBSBMSIM(150, 3, 0.1, 1, 10, 1, 2, 4, 5, 1, 10000, 1000, ...
          [1/2, 1/3, 1/6])

EBSBMSIM(nVertex, nBlock, epsilonInB, gStart, gEnd, ...
    modelTypeStart, modelTypeEnd, scaleCovarianceStart, ...
    scaleCovarianceEnd, nCore, nBurnIn, nConverge, rho, ...
    dimLatentPosition)

dimLatentPosition selects the dimension of latent positions.
dimLatentPosition should be integers between 1 and nBlock.
Default dimLatentPosition = nBlock.

EXAMPLE:
      EBSBMSIM(150, 3, 0.1, 1, 10, 1, 2, 4, 5, 1, 10000, 1000, ...
          [1/2, 1/3, 1/6], 2)

EBSBMSIM(nVertex, nBlock, epsilonInB, gStart, gEnd, ...
    modelTypeStart, modelTypeEnd, scaleCovarianceStart, ...
    scaleCovarianceEnd, nCore, nBurnIn, nConverge, rho, ...
    dimLatentPosition, isHomophily, isIdentifiable)

isHomophily controls the constraints. 
If isHomophily = 1, then it enforces homophily:
      <nu_i,nu_j> <= <nu_i,nu_i>  for any i, j
Default isHomophily = 1.

isIdentifiable controls the constraints.
If isIdentifiable = 1, then it enforces identifiability:
      <nu_i,nu_i> >= <nu_j,nu_j>  for any i > j
Default isIdentifiable = 0.

EXAMPLE:
      EBSBMSIM(150, 3, 0.1, 1, 10, 1, 2, 4, 5, 1, 10000, 1000, ...
          [1/2, 1/3, 1/6], 2, 1, 1)

EBSBMSIM(nVertex, nBlock, epsilonInB, gStart, gEnd, ...
    modelTypeStart, modelTypeEnd, scaleCovarianceStart, ...
    scaleCovarianceEnd, nCore, nBurnIn, nConverge, rho, ...
    dimLatentPosition, isHomophily, isIdentifiable, theta)

theta selects the hyperparameters for the prior distribution for rho.
Default theta = [1, 1, ..., 1] of length nBlock.

EXAMPLE:
      EBSBMSIM(150, 3, 0.1, 1, 10, 1, 2, 4, 5, 1, 10000, 1000, ...
          [1/2, 1/3, 1/6], 2, 1, 1, [1, 3, 2])
