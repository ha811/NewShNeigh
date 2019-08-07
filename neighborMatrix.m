function M = neighborMatrix(n,sig,type)
%generate neighbor matrix
if type==0
    M = spdiags(ones(n,2),[-1 1],n,n);
%     M(1:n+1:end) = 0;
else
    if sig == 0
        M = eye(n);
    else
        M = exp(-squareform(pdist((1:n)')).^2/2.0/sig/sig);
    end
end
end