function [posMatrix,retVals, shMat]=ShNeigh(heatmapMatrix,method_type)
tic;
alpha_in = 1.0;
[rMat,selbin]=matFilter(heatmapMatrix);
selbinIdx=find(selbin);
n=size(rMat,1);
sig_in=n*0.023;
cover=(nnz(rMat)-n)/n/(n-1);
rho_in=max((1-cover)*sqrt(n),min(3,0.2*sqrt(n)));

FreqMat=rMat;
if method_type==1  %Gaussian neighboring matrix,fixed alpha,fixed rho,fixed sigma
    [error_term,XYZ, shMat]=shNeigh_knownParams(FreqMat,selbin,alpha_in,rho_in,sig_in,1);
    retVals.alpha=alpha_in;
    retVals.rho=rho_in;
    retVals.sig=sig_in;
elseif method_type==2    %Gaussian neighboring matrix,variable alpha,fixed sigma
    fMpar = @(x)shNeigh_knownParams(FreqMat,selbin,x(1),x(2),sig_in,1);
    [alpharho,fval,exitflag,output] = fminsearch(fMpar,[alpha_in rho_in],optimset('MaxFunEvals',1000));
    [error_term,XYZ]=shNeigh_knownParams(FreqMat,selbin,alpharho(1),alpharho(2),sig_in,1);
    retVals.alpha=alpharho(1);
    retVals.rho=alpharho(2);
    retVals.sig=sig_in;
end

runningTime = toc;

posMatrix=[];     
for i=1:size(XYZ,2)
    posMatrix= [posMatrix;selbinIdx(i),XYZ(1,i),XYZ(2,i),XYZ(3,i)];
end
retVals.runningTime=runningTime;
retVals.freqerr=error_term;

end
