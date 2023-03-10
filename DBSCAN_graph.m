function [IDX, isnoise]=DBSCAN_graph(X,MinPts,Minscore)

    C=0;
    
    n=size(X,1);
    
    % IDX labels the cluster index for each point falling into that
    % cluster, other points will be assigned as 0
    IDX=zeros(n,1);
    
    
    
    visited=false(n,1);
    isnoise=false(n,1);
    
    for i=1:n
        if ~visited(i)
            visited(i)=true;            
            Neighbors=find(X(i,:) > Minscore);
            Neighbors = transpose(Neighbors);
            if numel(Neighbors)<MinPts
                % X(i,:) is NOISE
                isnoise(i)=true;
            else
                C=C+1;
                ExpandCluster(i,Neighbors,C);
            end
            
        end
    
    end
    
    function ExpandCluster(i,Neighbors,C)
        IDX(i)=C;
        
        k = 1;
        while true
            j = Neighbors(k);
            
            if ~visited(j)
                visited(j)=true;
                Neighbors2=find(X(j,:) > Minscore);
                Neighbors2 = transpose(Neighbors2);
                Lia = ismember(Neighbors2,Neighbors,'rows');
                Neighbors2(Lia,:) = [];
                if numel(Neighbors2)>=MinPts
                    Neighbors=cat(1,Neighbors,Neighbors2);   
                end
            end
            if IDX(j)==0
                IDX(j)=C;
            end
            
            k = k + 1;
            if k > numel(Neighbors)
                break;
            end
        end
    end

end



