function [XYZ,fval,exitflag] = MDSneigh(shMat,M,rho)
    exitflag = 1;
    n = length(shMat);
    A = -shMat.*shMat/2;
    H = eye(n)-ones(n,1)*ones(1,n)/n;
    gramMat = H*A*H;
    
    D = diag(sum(M));
    lapMat = D - M;
    w = rho*norm(gramMat,'fro')/max(0.000001,norm(lapMat,'fro'));
    BB= gramMat-w*lapMat;
 
    [V,lambdaMat]=eig(BB);
    lambda=diag(lambdaMat);
    if lambda(3) < 0
        exitflag=0;
    end
    XYZ0=[V(:,1)*sqrt(lambda(1)) V(:,2)*sqrt(lambda(2)) V(:,3)*sqrt(lambda(3))];
    XYZ = real(XYZ0');
    
    fval=0;
%     Dpred = squareform(pdist(XYZ'));
%     Apred = -Dpred.*Dpred/2;
%     gramPred = H*Apred*H;
%     fval = sum(sum((gramMat-gramPred).^2+w*M.*Dpred.*Dpred));
%    fval = trace((XYZ'*XYZ-BB)*(XYZ'*XYZ-BB))+trace(gramMat*gramMat-BB*BB);
    %fval = trace((XYZ'*XYZ-BB)*(XYZ'*XYZ-BB));
