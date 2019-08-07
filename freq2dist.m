function distMat = freq2dist(heatmapMatrix,alpha,beta)
if ~exist('beta','var'); beta = 1.0; end
distMat = (1./heatmapMatrix).^alpha;
distMat(isinf(distMat)) = 0;
distMat = triu(distMat,1)+triu(distMat,1)';
distMat = beta*distMat;
end