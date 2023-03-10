function [FreqMat,selbin]=matFilter(heatmapMatrix)
n=size(heatmapMatrix,1);
selbin=sum(heatmapMatrix>0)>2; %at least connect to two points
selbinIdx=find(selbin);
heatmapMatrix=heatmapMatrix(selbin,:);
heatmapMatrix=heatmapMatrix(:,selbin);

FreqMat=heatmapMatrix;
[num_parts,comps]=conncomp(digraph(FreqMat),'Type', 'weak');
%%make sure connected   
if num_parts>1
    compstat=tabulate(comps);
    [cm,ci]=max(compstat(:,2));
    finalbin=comps==ci;
    selbinIdx=selbinIdx(finalbin);
    FreqMat=FreqMat(finalbin,finalbin);
    selbin=false(n,1);
    selbin(selbinIdx)=true;
end