function [solOut] = integrateParasiteExperiments(p)
%This program runs the dynamical simulations for ATN with parasites.
%{
This is the model: (Boit et al 2012 without detritus or extra mortality
and! Also incorporating distribtuions from Brose 2006 and scalings
from WIlliams 2007: Homage to yodzis and Innes...
dB_i/dt = -f_m x_i B_i  maintenance
            + f_a x_i B_i\sum_j[y_{ij}F_{ij}(B) Benefit from eating
            - \sum_j (x_j Y_ji B_j F_ji(B))/e_{ji} Loss to being eaten

For basal species the first two terms are replaced by
r_iB_iG_i(B) -- (exudation fraction?)

F_ij as in Boit 2012, no interference terms:
 
w_ij Bj^q/(B_0 + \sum_resources w_il B_l^q)

G_i as in Boit 2012:

1- sum_produceers B_j/K

p is an object with all the parameters needed (don't want to calculate
properties in this function; also, want to have access to these after code
runs.
p. ...
             'S',S...
            ,'C',C...
            ,'res',res...
            ,'con',con...
            ,'K',K...
            ,'eij',eij...
            ,'wij',wij...
            ,'basal',basal...
            ,'r',r...
            ,'B0',B0...
            ,'h',h...
            ,'x',x...
            ,'halfSat',halfSat...
            ,'Tf',Tf ...
            ,'phi',.15...
            ,'extctThresh',1e-30...
            ,'AbsTol',1e-15...
            ,'RelTol',1e-5...
            ,'modelCode',modelCode...
            );


%}
%just use all of these and identify that as a potential drawback.
%Note that we don't need to identify concomittant links at all if we



includeFraction = p.modelCode(1)==2;
includeConcomittant = p.modelCode(2)==2;

res = p.res;
con = p.con;
S = p.S;
phi = p.para.*p.phi + ~p.para;

B0h = p.halfSat^p.h;
one = ones(p.S,1);
extctThresh = p.extctThresh;


wijmx = sparse(res,con,p.wij,p.S,p.S);


paraCon = p.para(con);

conM = p.M(con);
resM = p.M(res);

paraParaLink = paraCon & (conM<resM);

paraLinkPara = con(paraParaLink);
paraLinkHost = res(paraParaLink);

trophLink = ~paraParaLink;

trophicLinkCon = con(trophLink);
trophicLinkRes = res(trophLink);



yijParList = p.yij(paraParaLink);
yijTroList = p.yij(trophLink);


eijParList = p.eij(paraParaLink);
eijTroList = p.eij(trophLink);


wijParList = p.wij(paraParaLink);
wijTroList = p.wij(trophLink);

wijPar = sparse(paraLinkHost,paraLinkPara,wijParList,p.S,p.S);
wijTro = sparse(trophicLinkRes,trophicLinkCon,wijTroList,p.S,p.S);

x = p.x;
h = p.h;
K = p.K;
r = p.r;
para = p.para;
yij = p.yij;
eij = p.eij;
basal = p.basal;
%,one,B0h,p,wij,yij,includeFraction,eij,includeConcomittant,phi,wijTro,wijPar,yijTro,yijPar,eijTro,eijPar



%function [value,isterminal,direction] = extinction(~,B,S_,extctThresh_)
    function [value,isterminal,direction] = extinction(~,B)
        %This should be a smooth function! much easier on the solver.
        value = (log10(B)-log10(extctThresh));%*100; %/10
        value(B<0) = -1;
        isterminal = ones(S,1);
        direction = zeros(size(value));
        
    end


tnow = 0;
Tfinal = p.Tf;

try
    AbsTol = p.AbsTol;
    RelTol = p.RelTol;
catch
    AbsTol = extctThresh;
    RelTol = 1e-5;
end

options = odeset(p.options...
    ...,'Events',@(t,B) extinction(t,B,S,extctThresh)...
    ,'Events',@extinction...
    ,'AbsTol',AbsTol...
    ,'RelTol',RelTol...
    ,'NonNegative',uint8(1:40));
Tinterval = [tnow Tfinal];

B0 = p.B0;

sol = struct;

solOut.extctOrder = zeros(p.S,1);
B0s = zeros(p.S,p.S);
B0s(:,1) = B0;
solOut.extctTime = zeros(1,p.S);
extinct_sp = false(p.S,1);
countExtct = 1;
countSol = 0;



while tnow<Tfinal
    countSol = countSol+1;
    %could just use odextend
    sol = p.odeSolver(@(t,B) ATN_RHS(t,B,includeFraction,includeConcomittant,S,h,K,x,r,para,yij,yijTroList,yijParList,wijmx,wijTro,wijPar,eij,eijTroList,eijParList,res,con,paraLinkPara,paraLinkHost,trophicLinkRes,trophicLinkCon,paraParaLink,B0h,one,phi,basal),Tinterval,B0,options);
    %sol = p.odeSolver(@ATN_RHS,Tinterval,B0,options);
    tnow = sol.x(end);
    
    %Saving all solutions for diagnostic purposes.  Not necessary to save
    %the entire time series? The values at each extinction would be a good
    %compromise; can easily reconstruct particular parts of each sim. that
    %I want (but still, time-consuming to get the entire time series..).  I
    %just don't have the storage right now to make this work.
    
    solName = sprintf('sol%u',countSol);
    
    
    %Tfinal = tnow+1000; This *could* be a good idea
    if numel(sol.ie)>0
        Tfinal = max(Tfinal,ceil(tnow+p.Trelax));
    end
    %We need to catch species that get completely synchronized.  check:
    %correlation with extinct species?
    %relative error between extinct species?
    Tinterval = [tnow Tfinal];
    
    %The indices of the extinct species are here.
    extctIndices = sol.ie;
    %The array keeping track of extinct species is no longer up-to-date. Save the old array as extinct_sp0.
    extinct_sp0 = extinct_sp;
    %Update the new extpinct species list.
    extinct_sp(extctIndices) = true;
   
   %The new initial conditions are the final biomasses of this solution. 
    B0 = sol.y(:,end);
   
    %These are the biomasses of the dead species
    bDead = B0(sol.ie);
   
    %for each biomass of extinct species, 
    for ii = bDead'
        %Find the relative error between it's biomass and all other current biomasses.
        relErrori = abs(ii-B0)/ii;
        %whenever there is a species within relTol of the extinct species, mark it as extinct, too. The reason for this is that the solver can't actually tell these two bimoasses apart.  I was getting a problem where two species were pretty much locked together and going extinct at the same rate. One would hit the threshold first, then the second *should* have gone extinct in the very next timestep; but the solver pick up on it because it doesn't detect events in the first time step. So I had to manually kill that other species, too.
        extinct_sp(relErrori<RelTol)=true;
    end
   
    %Find the newly extinct species. 
    extctIndices = find(extinct_sp - extinct_sp0);
    
    %Zero out biomasses of the newly extinct species.
    B0(extinct_sp) = 0;
   
    %update the indices of the extinct species.
    sol.ie = extctIndices;
   
    %Mark the extinction event order; multiple species can go extinct at thsi time.
    solOut.extctOrder(sol.ie) = countExtct;
    
    %extctTime includes the initial condition, t=0.
    solOut.extctTime(countSol+1) = sol.x(end);
    
    %save the current solution structure.
    solOut.(solName) = sol;
    
    %update the current extinctions. multiple may go extinct at once. These should be snapshots of all the initial conditions.
    B0s(:,countSol+1) = B0;

    countExtct = countExtct+numel(extctIndices);
    
end
%delete the extra zeros; not all species go extinct.
solOut.extctTime(countSol+2:end) = [];
B0s(:,countSol+2:end) = [];

%save thesnapshots to the returned structure.
solOut.B0s = B0s;
solOut.n = countSol;

end


