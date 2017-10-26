%This code sets up parameters structures for a parallel implementation of
%the parasites experiment.

%Setting up the arrays and structures.

nWeb = 100;
models = fullfact([2,2,2,2]);
nModels = length(models);

kFrees = [1 2];   %Models(1): BSR exponents for free livers
kParas = [-3 -4];  %Models(2): BSR exponents for free livers
fracFrees = [0 1]; %Models(3): including fraction of free living time
fracParas = [0 1]; %Models(4): including concomittant links

fParAll0 = [0 0.025 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.40 0.45 0.50 0.75 1];
nfPar = numel(fParAll0);

nSims = nWeb*2+nWeb*(nModels)*(nfPar-1);

simParams = cell(nSims,1);
[simParams{:}] = deal(struct('web',0 ...
    ,'fPar',0 ...
    ,'kFree',0 ...
    ,'kPara',0 ...
    ,'fracFree',0 ...
    ,'fracPara',0 ...
    ,'modelCode',[0 0]...
    ,'para',zeros(40,0)...
    ,'LL',[]...
    ,'B0',zeros(40,1)...
    ,'gr',zeros(40,1)...
    ,'parOrder',zeros(40,1)...
    ,'TL',zeros(40,1)...
    ));

%%% Parameters of the Niche Webs
S = 40;
C = .15;

simNo = 0;
for ii = 1:nWeb %Generate new niche webs each time this code is run:
    webBad = true;

    while webBad
        [res, con,~,~,~] = NicheModel_nk(S,C);
        simMx = calculateSimilarity(res,con);
        mx = sparse(res,con,1,S,S);
        webBad = max(max(simMx))==1;
    end
    
    %Need to decide if each species is a basal.
    basal = false(S,1);
    
    for kk = 1:S
        if sum(con==kk)==0
            basal(kk) = true;
        end
    end
    
    B0 = .95*rand(S,1)+.05;
    gr = basal.*(randn(S,1)*.1+1);
    
    SList = 1:S;
    idxPar = datasample(SList(~basal),sum(~basal),'Replace',false);
    
    %short-weighted trophic level
    A = full(sparse(res,con,1,S,S));
    nonCannA = A;
    nonCannA(diag(true(S,1)))=0;

    %Formula for the prey-averaged trophic level.
    paTL_mx = sparse(nonCannA)*(diag(1./sum(sparse(nonCannA))));
    paTL = (speye(S)-paTL_mx')\ones(S,1);

    for  model = models'
        
        kFree = kFrees(model(1));
        kPara = kParas(model(2));
        
        fracFree = fracFrees(model(3));
        fracPara = fracParas(model(4));

        for fPar = fParAll0
            %omit extraneous simulations.
            if (fPar==0)&&((fracFree~=0)||(fracPara~=0)||(kPara~=-3))
                continue
            end
            
            %Number of parasites at this level. went with round to keep as
            %close to what was advertised as possible.
            nBasal = sum(basal);
            nFree = S-nBasal;
            para = false(S,1);
            nPar = round(fPar*nFree);
            para(idxPar(1:nPar)) = true;
            
            simNo = simNo+1;
            simParams{simNo}.web = ii;
            simParams{simNo}.fPar = fPar;
            simParams{simNo}.kFree = kFree;
            simParams{simNo}.kPara = kPara;
            simParams{simNo}.fracFree = fracFree;
            simParams{simNo}.fracPara = fracPara;
            simParams{simNo}.modelCode = [fracFree fracPara];
            simParams{simNo}.para = para;
            simParams{simNo}.LL = [res con];
            simParams{simNo}.B0 = B0;
            simParams{simNo}.gr = gr;
            simParams{simNo}.parOrder = idxPar;
            simParams{simNo}.patl = paTL;
            
        end
        
    end
end
save('../setup_parameters.mat','simParams')

TS = zeros(40,1000,nSims);
extcts = cell(nSims,1);

%Web Structure parameters
S = 40;
C = 0.15;

%Metabolic parameters for body sizes, metabolic rates
axFree = .314;
axPara = .314;
metabScaling = 0.25;

%%% Parameters needed in dynamical equations
%maximum production rate of consumers relative to metabolic rate
yFree = 8;
yPara = 8;

%shared carrying capacity for basal species
K = 5;

%Hill coefficient and half-saturation density for the functional response
h=1.2;
halfSat = 0.5;

%Simulation time is: Max(Tf,t_(last Extinction)+Trelax)
Tf = 10000;
Trelax = 2000;

%%% Simulation parameters
%Specifying threshold for species extinction and tolerance for solvers.
%These are necessarily related; must have abstol smaller that extctThresh
%so that we can have confidence in the species actually going extinct;
%otherwise we are guaranteed no accuracy in the species that goes extinct!
extctThresh = 1e-10;
AbsTol = 1e-16;
RelTol = 1e-6;

%This is the structure that passes into the integrate parasites function.
%contains all parameters needed to run the simulation.

params = struct(... 
    ... Web Parameters:
    'S',S...
    ,'C',C...
    ... DS Parameters:
    ,'K',K...
    ...FR Parameters:
    ,'halfSat',halfSat...
    ,'phi',.15...
    ,'h',h...
    ... Link Level properties:
    ,'res',[]...
    ,'con',[]...
    ,'eij',[]...
    ,'wij',[] ...
    ,'yij',[]...
    ... Species level Properties:
    ,'basal',zeros(40,1)...
    ,'r',zeros(40,1)...
    ,'B0',zeros(40,1)...
    ,'x',zeros(40,1)...
    ,'M',zeros(40,1)...
    ,'modelCode',[0 0]...
    ,'para',zeros(40,1)...
    ... Simulation Parameter:
    ,'Tf',Tf ...
    ,'Trelax',Trelax ...
    ... Solver Parameters:
    ,'extctThresh',extctThresh ...
    ,'AbsTol',AbsTol...
    ,'RelTol',RelTol...
    ,'odeSolver',@ode45... .
    ,'options',odeset()...
    );
%This automatically makes the parallel pool use ever core possible. on ocelote, their
%default is (for some reason?) 12.

nCores = feature('numcores');
parpool('local',nCores)

parfor  (simNo = 1:nSims, nCores)
   
    p = params;
    simParam = simParams{simNo};
    
    %%% Define properties for this model
    kFree = simParam.kFree;
    kPara = simParam.kPara;
    modelCode = simParam.modelCode;
    
    
    LL = simParam.LL;
    
    res = LL(:,1);
    con = LL(:,2);
    
    patl = simParam.patl;
    basal = simParam.gr>0;
    
    nBasal = sum(basal);
    nFree = S-nBasal;
    
    %Assimilation efficiency is given in Brose et al as .45 and .85 for
    %herbivores and carnivores, respectively.  This would correspond to
    %the links we talkin bout.  Easy enough to figure:
    eij = zeros(size(res));
    eij(basal(res)) = .45;
    eij(~basal(res)) = .85;  %Future: add an eij for parasitic links
    wij = true(size(res));
    
    
    %Set the properties that we have; web-level properties, independent
    %of parasite identity.
    p.B0 = simParam.B0;
    
    p.res = res;
    p.con = con;
    p.eij = eij;
    p.wij = wij;
    p.basal = basal;
    p.r = simParam.gr;
    p.modelCode = modelCode;
    
    %Pre-allocate final data.
    yFinals = zeros(S,1);
    yMeans = zeros(S,1);
    yStds = zeros(S,1);
    
    ZFree = 10^kFree;
    
    %%% Define the parameters that change with parasite identities.
    ax = zeros(S,1);
    y = zeros(S,1);
    M = zeros(S,1);
    para = simParam.para;
    free = ~para;
    %assimilation rates
    y(free) = yFree;
    y(para) = yPara;
    %scaling constant
    ax(free) = axFree;
    ax(para) = axPara;
    %bodymass
    M(free) = ZFree.^(patl(free)-1);
    M(para) = 10.^(kPara + kFree*(patl(para)-2));
    %metabolic rate
    x = ax.*M.^(-0.25);
    x(basal) = 0;
    
    %Save parameters to parameter structre.
    p.x = x;
    p.para = para;
    p.M = M;
    p.yij = y(con);
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Integrate!
    sol = integrateParasiteExperiments(p);    
    %Save the extinction times and B0 snapshots. Get the intial conditions and the final conditions as well. Doesn't matter,r eally! 
    extcts{simNo} = [sol.extctTime;sol.B0s];
    endSol = sol.(sprintf('sol%u',sol.n));
    xEnd = endSol.x(end);
    xRange = linspace(xEnd,xEnd-999,1000);
    TS(:,:,simNo) = deval(endSol,xRange);
    %sol is too large to save. a single run can be as much as 77MB. Which
    %obviously doesn't scale well for 21000 runs.
    
    
end

save('../out.mat','extcts','TS','-v7.3');

