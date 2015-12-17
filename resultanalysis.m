
clear all

%% --- Simulation ---

% Parameter Setting
nVertex = 150;
modelTypeStart = 1;
modelTypeEnd = 3;
scaleCovarianceStart = 1;
scaleCovarianceEnd = 5;
epsilon = 0.1;
gStart = 1;
gEnd = 1000;

error = cell(3, 5);
errorMean = zeros(3, 5);
errorCI = zeros(3, 5, 2);
errorMedian = zeros(3, 5);

for modelType = modelTypeStart:modelTypeEnd
    for scaleCovariance = scaleCovarianceStart:scaleCovarianceEnd
        ind = [];
        for iGraph = gStart:gEnd
            if exist(['./results/results-SBM-model' num2str(modelType) ...
                    '-scale' num2str(scaleCovariance) '-sim-n' ...
                    num2str(nVertex) '-eps' num2str(epsilon) '-graph' ...
                    num2str(iGraph) '.mat'])
                ind = [ind iGraph];
            end
        end
        
        maxIter = length(ind);
        
        error{modelType, scaleCovariance} = zeros(1, maxIter);
        
        for iInd = 1:maxIter
            load(['./results/results-SBM-model' num2str(modelType) ...
                '-scale' num2str(scaleCovariance) '-sim-n' ...
                num2str(nVertex) '-eps' num2str(epsilon) '-graph' ...
                num2str(ind(iInd)) '.mat']);
            error{modelType, scaleCovariance}(iInd) = errorRate;
        end
        
        errorMean(modelType, scaleCovariance) = ...
            mean(error{modelType, scaleCovariance});
        errorCI(modelType, scaleCovariance, :) = ...
            [mean(error{modelType, scaleCovariance}) - ...
            1.96*std(error{modelType, scaleCovariance})/sqrt(maxIter), ...
            mean(error{modelType, scaleCovariance}) + ...
            1.96*std(error{modelType, scaleCovariance})/sqrt(maxIter)];
        errorMedian(modelType, scaleCovariance) = ...
            median(error{modelType, scaleCovariance});
    end
end

%% --- Plot ---

n = 1:5;

figure
set(gcf, 'Color', [1, 1, 1]);
plot(n, errorMean(1, :), 'b-', 'LineWidth', 2);
axis([0, 6, 0, 0.4]);
hold on;

plot(n, errorMean(2, :), 'r-', 'LineWidth', 2);
plot(n, errorMean(3, :), 'g-', 'LineWidth', 2);

hold all;
leg = legend('Model = Gold', 'Model = ASGE', 'Model = ASGE1');
legend('boxoff');

set(leg, 'location', 'eastoutside', 'FontSize', 11)
set(gca, 'ygrid', 'on')

plotshadedarea(n, squeeze(errorCI(1, :, :))', 'b');
plotshadedarea(n, squeeze(errorCI(2, :, :))', 'r');
plotshadedarea(n, squeeze(errorCI(3, :, :))', 'g');
hold off

xlabel('Scale of Covariance Matrices', 'FontSize', 12);
ylabel('Classification Error', 'FontSize', 12)

% set(gca,'XTick',[0:3]) % This automatically sets 
labels = {'', '0', 'n', 'n*n_k', 'n^2', 'Infinity', ''};
set(gca,'XTickLabel',labels)

%% --- Paired Test ---

% Alternative hypothesis: pair1 < pair2.

modelTypePair1 = 2;
scaleCovariancePair1 = 2;

modelTypePair2 = 2;
scaleCovariancePair2 = 3;

indPair = [];
for iGraph = gStart:gEnd
    if exist(['./results/results-SBM-model' num2str(modelTypePair1) ...
            '-scale' num2str(scaleCovariancePair1) '-sim-n' ...
            num2str(nVertex) '-eps' num2str(epsilon) '-graph' ...
            num2str(iGraph) '.mat']) && ...
            exist(['./results/results-SBM-model' num2str(modelTypePair2)...
            '-scale' num2str(scaleCovariancePair2) '-sim-n' ...
            num2str(nVertex) '-eps' num2str(epsilon) '-graph' ...
            num2str(iGraph) '.mat'])
        indPair = [indPair iGraph];
    end
end

maxIterPair = length(indPair);

errorPair = zeros(2, maxIterPair);

for iIndPair = 1:maxIterPair
    load(['./results/results-SBM-model' num2str(modelTypePair1) ...
        '-scale' num2str(scaleCovariancePair1) '-sim-n' ...
        num2str(nVertex) '-eps' num2str(epsilon) '-graph' ...
        num2str(ind(iIndPair)) '.mat']);
    errorPair(1, iIndPair) = errorRate;
    load(['./results/results-SBM-model' num2str(modelTypePair2) ...
        '-scale' num2str(scaleCovariancePair2) '-sim-n' ...
        num2str(nVertex) '-eps' num2str(epsilon) '-graph' ...
        num2str(ind(iIndPair)) '.mat']);
    errorPair(2, iIndPair) = errorRate;
end


% 1-sided sign-test
tmpStats = sum(errorPair(1, :) < errorPair(2, :));
pValue = 1 - binocdf(tmpStats - 1, maxIterPair, 0.5)

tmpStats = round((errorPair(1, :) - errorPair(2, :))*nVertex);
hist(tmpStats, -60:1:60)



[n, xout] = hist(tmpStats);
bar(xout, n, 1)




