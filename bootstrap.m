% Bootstrap Method
% This function takes N bootstrap samples with replacement of the
% original data and stores the calculation of the desired statistic
% in the returned array.
%
% bootstrap(data,N)
function [stat]=bootstrap(data1, data2, N, nConverge)

stat=zeros(N,1); %initialize array to be returned of bootstrap estimations


for k = 1:N
    iboot1 = randsample(nConverge, nConverge, true);
    iboot2 = randsample(nConverge, nConverge, true);
    stat(k) = psrf(data1(iboot1), data2(iboot2));
end
