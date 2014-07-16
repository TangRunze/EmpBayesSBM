
clear all

%% --- Simulation ---

% Parameter Setting

modelType = 2;
scaleCovariance = 2;
epsilon = 0.1;
gStart = 1;
gEnd = 100;

ind = [];
for iGraph = gStart:gEnd
    if exist(['./result/results-SBM-model' num2str(modelType) '-scale' ...
            num2str(scaleCovariance) '-sim-eps' num2str(epsilon) ...
            '-graph' num2str(iGraph) '.mat'])
        ind = [ind iGraph];
    end
end

maxIter = length(ind);

if ((modelType == 2) && (scaleCovariance == 2))
    errorAsge = zeros(1, maxIter);
    for iInd = 1:maxIter
        load(['./result/results-SBM-model' num2str(modelType) '-scale' ...
            num2str(scaleCovariance) '-sim-eps' num2str(epsilon) ...
            '-graph' num2str(ind(iInd)) '.mat']);
        % errorAsge(iInd) = errorRate;
        errorAsge(iInd) = error_rate;
    end
    mean(errorAsge)
    [mean(errorAsge) - 1.96*std(errorAsge)/sqrt(maxIter), ...
        mean(errorAsge) + 1.96*std(errorAsge)/sqrt(maxIter)]
    median(errorAsge)
end


