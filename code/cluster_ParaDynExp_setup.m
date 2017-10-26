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
[simParams{:}] = deal(struct('web',struct() ...
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
fprintf('Setting up parameters ...\n')
for ii = 1:nWeb
    fprintf('Generating Web no. %4u ...',ii)
    %Generate new niche webs each time this code is run:
    webBad = true;
    while webBad
        [res, con,~,~,~] = NicheModel_nk(S,C);
        simMx = calculateSimilarity(res,con);
        mx = sparse(res,con,1,S,S);
        webBad = max(max(simMx))==1;
    end
    fprintf('Done.')
    
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
    fprintf(repmat('\b',1,32));
end
fprintf('Done.\n')
save('../simulationDataStructures.mat','simParams')

