function [connected] = walk(source,connected,matrix) 

%This function finds all nodes connected to the source using adjacency
%matrix 

neighbors = find(matrix(source,:));

connected = [connected;source];
for kk = neighbors
    
    if sum(kk == connected)
        continue
    else
        connected = walk(kk,connected,matrix);
    end
    
    
end