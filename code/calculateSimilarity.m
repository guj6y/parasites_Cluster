function [simMx] = calculateSimilarity(res,con)

N = max([res;con]);

A = sparse(res,con,1,N,N);

%Maximum Similarity: need to compare all species to each other -gross! >.<

%Thought - this is obviously a way to measure competition.  Probably a
%reference about that?  Are species with high similiarity (or just
%prey/pred similarity) more or less likely to go extinct?

%NK's attempt to improve max sim with ridiculous 3d arrays:
coPreyMxs = zeros(N,N,N);
coPredMxs = zeros(N,N,N);
anyPreyMxs = zeros(N,N,N);
anyPredMxs = zeros(N,N,N);

%might be possible to reduce complexity with clever indexing.  Still, less
%than a second for up to about 250x250 matrices.  This might also be slower
%depending on the machine it is on - less RAM might mean this doesn't even
%work.  Probalby fine on any modern laptop.
%Method of Coralie & Perrine (significantly slower when I tested it):
%{
S_ij = zeros(N,N);
for ii = 1:(N-1)
    for jj = (ii+1):N
        %number Prey in common:
        co_prey = sum((A(:,ii)>0) & (A(:,jj))>0);
        %Predators in common:
        co_pred = sum((A(ii,:)>0) & (A(jj,:)>0));
        %Total links for both species:
        tot_prey = sum((A(:,ii)>0) | (A(:,jj))>0) +...
            sum((A(ii,:)>0) | (A(jj,:)>0));
        %Saving the similarity.
        S_ij(ii,jj) = (co_prey + co_pred)/tot_prey;
        S_ij(jj,ii) = (co_prey + co_pred)/tot_prey;
    end
end
%}


for ii = 1:N
    
    %co_prey_mxs(ii,jj,kk) = 1 means species kk (3rd dim) and species jj
    %(2nd dim) both eat species ii (1st dim).
    coPreyMxs(:,:,ii) = ((A(:,ii)*ones(1,N))>0)&(A>0);
    
    %co_pred_mxs(ii,jj,kk) = 1 means species kk (3rd dim) and species ii
    %(1st dim) both are eaten by species jj (2nd dim).
    coPredMxs(:,:,ii) = ((A(ii,:)'*ones(1,N))>0)'&(A>0);
    
    
    %any_prey_mxs(ii,jj,kk) = 1 means species kk (3rd dim) or species jj
    %(2nd dim) bot eat species ii (1st dim).
    anyPreyMxs(:,:,ii) = ((A(:,ii)*ones(1,N)>0)|(A>0));
    
    
    %any_pred_mxs(ii,jj,kk) = 1 means species kk (3rd dim) or species ii
    %(1st dim) is eaten by species jj (2nd dim).
    anyPredMxs(:,:,ii) = ((A(ii,:)'*ones(1,N)>0)'|(A>0));
    
    
    %Tried reducing the complexity by filling in at 'right angles' ( fill
    %the 3rd dimension, then a strip in one of the others at each time
    %step).  Ended up increasing computation time ever so slightly.  oh
    %well.  ¯\_(?)_/¯
    %{
    %co_prey_mxs(ii,jj,kk) = 1 means species kk (3rd dim) and species jj
    %(2nd dim) both eat species ii (1st dim).
    co_prey_mxs(:,ii:N,ii) = ((A(:,ii)*ones(1,(N-ii+1)))>0)&(A(:,ii:N)>0);
    co_prey_mxs(:,ii,(ii+1):N) = co_prey_mxs(:,(ii+1):N,ii);
    %co_pred_mxs(ii,jj,kk) = 1 means species kk (3rd dim) and species ii
    %(1st dim) both are eaten by species jj (2nd dim).
    co_pred_mxs(ii:N,:,ii) = ((A(ii,:)'*ones(1,(N-ii+1)))>0)'&(A(ii:N,:)>0);
    co_pred_mxs(ii,:,(ii+1):N) = co_pred_mxs((ii+1):N,:,ii)';
    
    %any_prey_mxs(ii,jj,kk) = 1 means species kk (3rd dim) or species jj
    %(2nd dim) bot eat species ii (1st dim).
    any_prey_mxs(:,ii:N,ii) = ((A(:,ii)*ones(1,(N-ii+1)))>0)|(A(:,ii:N)>0);
    any_prey_mxs(:,ii,(ii+1):N) = any_prey_mxs(:,(ii+1):N,ii);
    
    %any_pred_mxs(ii,jj,kk) = 1 means species kk (3rd dim) or species ii
    %(1st dim) is eaten by species jj (2nd dim).
    any_pred_mxs(ii:N,:,ii) = ((A(ii,:)'*ones(1,(N-ii+1)))>0)'|(A(ii:N,:)>0);
    any_pred_mxs(ii,:,(ii+1):N) = any_pred_mxs((ii+1):N,:,ii)';
    %}
end
%Had to fiddle a *bit* to get the dimensions right.. but this agrees with
%our french colleagues' method.
simMx = (permute(sum(coPreyMxs),[2,3,1]) +...
    permute(sum(coPredMxs,2),[1,3,2]))./(...
    (permute(sum(anyPreyMxs),[2,3,1])) ...
    + permute(sum(anyPredMxs,2),[1,3,2]))...
    -eye(N);


end