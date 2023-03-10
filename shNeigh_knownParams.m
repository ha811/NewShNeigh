function [error_term, XYZ, shMat]=shNeigh_knownParams(heatmapMatrix,selbin,alpha,rho,sig,nmtype)
if ~exist('nmtype','var'); nmtype = 1; end
n = length(selbin);
M = neighborMatrix(n,sig,nmtype);
M=M(selbin,selbin);
distMat = freq2dist(heatmapMatrix,alpha,1.0);
shMat = distances(graph(distMat));
XYZ = MDSneigh(shMat,M,rho);
error_term = FitFreq(squareform(pdist(XYZ')),heatmapMatrix,-alpha);
end
