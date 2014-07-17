function  parsave(fName, errorRate, tau, tauConfidence, tauResult)
% save data
    save(fName, 'errorRate', 'tau', 'tauConfidence', 'tauResult')
end