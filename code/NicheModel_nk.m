function [res,con,n_new,c_new,r_new]= ...
     NicheModel_nk(N, C)                
%------------------------------------------------------------------

%globalStream = RandStream.getGlobalStream;
%reset(globalStream);

tries=10000;                                                               
%error on the connectances
error_n=0.05;

%validation of the plain niche web:
ok_n=0; 

while ((N>0) && (tries>0) && (ok_n==0))
    tries=tries-1;
    %assign niche values from a uniform distribution
    n = rand(N,1);  
    
    
    
    
    %designate range for each species

    %parameters for beta distribution:
    alpha = 1;
    beta = (1-2*C)/(2*C); 
    r = betarnd(alpha,beta,N,1); 
    
    %vector of ranges: 
    r = r.*n;  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % set center of range, uniformly distributed in [r_i/2,n_i]; 
    c = min(1-r./2,rand(N,1).*(n-r./2)+r./2);

    %sort everything:
    [n_new, Indx] = sort(n);                                                
    %n_new: niche values in ascending order
    %Indx: indices of species in descending order 
    %(-> 1 is the index of the smallest niche range, 10 is the index of the
    %largest)

    %the smallest r to highest index species 
    r_new = r(Indx);                                                       %NK: Maybe not how I'd do this
    c_new = c(Indx); 

    r_new(1) = 0; %change the r of highest index species to 0
    %so we have a basal species in every web
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
    
    %lower border of niche range for every prey:
    preymins = c_new - r_new/2;
    
    %upper border of niche range for every predator:
    preymaxs = c_new + r_new/2;

    %fills the empty matrix with niche ranges:
    n_mx = n_new*ones(1,N); 
    %matrix with the lowest points of ranges in every column:
    preymins_mx = ones(N,1)*preymins'; 
    %same, with highest:
    preymaxs_mx = ones(N,1)*preymaxs';

    %Construct the web adjacency matrix;
    %if species in the row is in the diet range of the species in the 
    %column, it gets eaten (if species i is eaten by species j, set (i,j) =
    %1).
    
    web_mx=(n_mx>=preymins_mx)&(n_mx<=preymaxs_mx);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  if there is an isolated species                                    %
    %  or something not connected to a basal                               %
    %  or a disconnected group                                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %Make sure the graph is weakly connected (no isolated species)
    weak_comp = graphconncomp(sparse(web_mx),'Directed',1,'Weak',1);
    web_connected = (weak_comp == 1);
    
    %make sure that all species have a directed connection from a basal
    %species; only if connected.
    
    %'Walk' is where most of the time is spent.  At this point, however, it
    %takes about .04s to run on a 100 species, .12 connectance web.  IF you
    %need to, you could use the shortest paths matrix from floyd warshall
    %to figure out if everything goes to a basal.  Could maybe also figure
    %out if it is connected using the sp stuff - one function call for both
    %things might save some time.
    
    if web_connected
        unconnected_species = (1:N)';
        basal_species = unconnected_species(sum(web_mx)==0)';
        connected_species = [];

        for kk = basal_species
            connected_species = walk(kk,connected_species,web_mx);
        end

        unconnected_species(connected_species) = [];
        
        if isempty(unconnected_species)
            web_connected = 1;
        else
            web_connected = 0;
        end
    end
    
    
    
    
    if web_connected
        % indices of links (which element in web_mx?)
        links = sum(sum(web_mx));  
        % Actual connectance
        C_web = links/(N^2);  
        if (abs(C_web-C)*1.0/C) > error_n
            ok_n=0;
        else
            ok_n=1;
        end
    end  

end

    [res,con] = find(web_mx);

    
end