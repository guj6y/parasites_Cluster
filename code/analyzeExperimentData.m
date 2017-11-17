%THis script computes statistics for the data that were extracted from the
%text output files from the experiment data.

%Modified from something that runs on Lucille2, meant to process data from the cluster, 
%on a chivo machine from math department.

%Someday I might setup the experiment to read in a file with the
%parameters/ experimental setup....  Then that file works for this, too.
%Goals, I guess.
kFree = [1 2];
kPara = [-3,-4];


%Load the data:
%"Rule of thumb to avoid poofing variables into the workspace"
%These are the names of the variables as I saved them over in
%readParasiteExperimentData.m:
var1 = 'simParams';
var2 = 'finalBiomasses';
var3 = 'meanBiomasses';
var4 = 'stdBiomasses';
%A little info about these things:
%It isn't at all clear what you load here. So:
%out_4d.mat contains the three 4-d arrays, finalBiomass, meanBiomass, and stdBiomass.
%setup_parameters.mat contains the cell array simParams. each element of the cell array is a structure. This is the structure:
%{
          web: 1                 %This is the index of tyhe niche web used
         fPar: 0                 %This is the fraction of parasites 
        kFree: 1                 %This is exponent  on free-living BSR
        kPara: -3                %This is exponent on parasite BSR
     fracFree: 0                 %This is host refuge identifier   
     fracPara: 0                 %This is the concomittant identifier
    modelCode: [0 0]             %This is redundant.
         para: [40×1 logical]    %This tells which species are parasites
           LL: [250×2 double]    %This is the link list
           B0: [40×1 double]     %THis is the initial Biomass.
           gr: [40×1 double]     %Growth rates of basal species
     parOrder: [1×33 double]     %This is the order of parasites.
           TL: [40×1 double]     %This is nothing.
         patl: [40×1 double]     %This is the prey-averaged trophic level.
%}

%setup parameters.
T = load('../simParams.mat');
webStructures = T.(var1);

%The results
T = load('../out.mat');
meanBiomass = T.(var2);
finalBiomass= T.(var3);
stdBiomass = T.(var4);

clear T

%So much for poofing! This has the basic sort of sim data that I gotta have.
load '../metaSimData.mat'

extctBinArray = finalBiomass>0;

%Want to get an array that lets me extract parasites.  A little bit
%challenging to work through the 4D array, but it works out.
paraSp = false(S,nFPar,1,nWeb);

%Binary array that tells me which species are basal.  (generality = 0).
%This array is 40x1x100.
basalSp = false(S,1,nWeb);

nBasal = zeros(nFPar,1,nWeb);
nFree = zeros(nFPar,1,nWeb);
nPara = zeros(nFPar,1,nWeb);

allLinkTypes = zeros(nWeb,nFPar,nModel,4,2);

bsrAll = cell(nWeb,1);
typesAll = cell(nWeb,1);

bsrFree = zeros(S,2);
bsrPara = zeros(S,2);

%The order of body size configurations in the third loop below:
bsConf = fullfact([2,2]);
bsAll = zeros(S,nFPar,4,S);

nSims = numel(webStructures);
%For each web...
for ii = 1:nWeb
    simNo = nSims/nWeb*(ii-1)+1;
    %Load the parasite orders
    par_ii = webStructures{simNo}.parOrder;    
    % Load the link list
    LL = webStructures{simNo}.LL;
    
    %L = number of links
    [L,~] = size(LL);

    %Pre-allocate arrays
    bsrAll{ii} = zeros(L,nFPar,4);
    bs_kk = zeros(S,nFPar);
    typesAll{ii} = zeros(L,nFPar);
    
    %For each fraction of parasites...
    
    %Need to calculate body sizes for each web.  Do it for all parasites
    %and all free-livers, then can switch as necessary. Note that this is MUCH nicer with     %as of r2016b... chivo has r2016a, fml.
    patl = webStructures{simNo}.patl;
    bsFree = permute(repmat(10.^(repmat(kFree,[S,1]).*(repmat(patl,[1,2])-1)),[1,1,nFPar]),[1 3 2]);
    bsPara = permute(repmat(10.^(repmat(kPara,[S,1]) + repmat(kFree,[S,1]).*(repmat(patl,[1,2])-2)),[1,1,nFPar]),[1 3 2]);
    
    basalSp(:,1,ii) = webStructures{simNo}.gr>0;
    nBasal(:,1,ii) = sum(basalSp(:,1,ii));
    nFree(:,1,ii) = S-nBasal(:,1,ii);
    nPara(:,1,ii) = round(fParAll'.*nFree(:,1,ii));
    nFree(:,1,ii) = nFree(:,1,ii) - nPara(:,1,ii);
    
    for jj = 1:nFPar
        %Complicated here.  wanted to avoid the double loop structure, but
        %it runs in O(.01s) so it doesn't really matter.  The first
        %dimension represents the species, so I pick the appropriate number
        %of parasites in the list and update thsoe indices, for that
        %fraction of parasites, for that web.
        paraSp(par_ii(1:nPara(jj,1,ii)),jj,:,ii) = true;
        
        %Each type of link; f->p, f->f ,etc.  Set up so that each type gets
        %a number 1,2,3, or 4; deal with each individually in the loop
        %below
        typesAll{ii}(:,jj) = paraSp(LL(:,1),jj,:,ii)*2 + paraSp(LL(:,2),jj,:,ii) + 1;
        
        
        for kk = 1:4
            allLinkTypes(ii,jj,:,kk,1) = sum(typesAll{ii}(:,jj)==kk);
            allLinkTypes(ii,jj,:,kk,2) = sum(repmat(typesAll{ii}(:,jj)==kk,[1 1 nModels]).*(extctBinArray(LL(:,1),jj,:,ii).*extctBinArray(LL(:,2),jj,:,ii)));
            
        end
        
    end

    for kk = 1:4
%         kFree = [1 2];
%         kPara = [-3,-4];

        freeSize = bsFree(:,:,bsConf(kk,1));
        paraSize = bsPara(:,:,bsConf(kk,2));
        bs_kk(paraSp(:,:,1,ii)) = paraSize(paraSp(:,:,1,ii));
        bs_kk(~paraSp(:,:,1,ii)) = freeSize(~paraSp(:,:,1,ii));
        
        bsAll(:,:,kk,ii) = bs_kk;
        bsrAll{ii}(:,:,kk) = bs_kk(LL(:,2),:)./bs_kk(LL(:,1),:);
    end
end

%Transform the cell arrays to normal arrays.  We lose identiy of individual
%webs by doing this.
bsrAll = cell2mat(bsrAll);
typesAll = cell2mat(typesAll);


% for ii = 1:3
%     h = figure;
%     bsConfiguration = 4;
%     histogramFPar = 1+5*(ii-1);
%     [n, edges] = histcounts(log10(bsrAll(:,histogramFPar,bsConfiguration)));
%     n1 = histcounts(log10(bsrAll(typesAll(:,histogramFPar)==1,histogramFPar,bsConfiguration)),edges);
%     n2 = histcounts(log10(bsrAll(typesAll(:,histogramFPar)==2,histogramFPar,bsConfiguration)),edges);
%     n3 = histcounts(log10(bsrAll(typesAll(:,histogramFPar)==3,histogramFPar,bsConfiguration)),edges);
%     n4 = histcounts(log10(bsrAll(typesAll(:,histogramFPar)==4,histogramFPar,bsConfiguration)),edges);
%     bar(edges(2:end),[n1',n2',n3',n4'],1,'stacked')
% %      hold on
% %      histogram('binEdges',edges,'BinCounts',n1)
% %      histogram('binEdges',edges,'BinCounts',n2)
% %      histogram('binEdges',edges,'BinCounts',n3)
% %      histogram('binEdges',edges,'BinCounts',n4)
% %      hold off
%     legend('Predation','Parasitism','Consumption of Parasites','Hyper-Parasitism',...
%         'Location','NorthWest')
%     xlabel('$\log_{10}$Body Size Ratio','Interpreter','LaTex')
%     ylabel('Frequency','Interpreter','LaTex')
%     saveas(h,...
%         sprintf('../../../Papers/ParasiteSwitching/floats/results/bsrStackHist%u',ii),'png')
%     close(h)
% end


fracLinkChange = allLinkTypes(:,:,:,:,2)./allLinkTypes(:,:,:,:,1);

paraSp = repmat(paraSp, [1,1,nModel,1]);
livingPara = sum(paraSp.*finalBiomass);

basalSp = permute(repmat(basalSp,[1,nFPar,1,nModel]),[1 2 4 3]);
freeSp = ~(paraSp|basalSp);


freePersistence = permute(squeeze(sum((finalBiomass>0)&freeSp)),[3,1,2])./...
    permute(repmat(nFree,[1,nModel,1]),[3,1,2]);

paraPersistence = permute(squeeze(sum((finalBiomass>0)&paraSp)),[3,1,2])./...
    permute(repmat(nPara,[1,nModel,1]),[3,1,2]);

basalPersistence = permute(squeeze(sum((finalBiomass>0)&basalSp)),[3,1,2])./...
    permute(repmat(nBasal,[1,nModel,1]),[3,1,2]);

%Getting net extinctions per parasite added
extctDiff = [nan(S,1,nModel,nWeb) diff(extctBinArray,1,2)];
netExtctPerPara = permute(squeeze(sum(extctDiff))./repmat([nan(1,1,nWeb); -diff(nPara)],[1,nModel,1]),[3 1 2]);


%Getting all persistence
persistenceAll = permute(squeeze(mean(finalBiomass>0)),[3 1 2]);
persistencePara = permute(squeeze(sum(finalBiomass>0&paraSp)./sum(paraSp)),[3 1 2]);
persistenceFree = permute(squeeze(sum(finalBiomass>0&freeSp)./sum(freeSp)),[3 1 2]);
persistenceBasal = permute(squeeze(sum(finalBiomass>0&basalSp)./sum(basalSp)),[3 1 2]);

%Try this
offset = repmat(persistenceAll(:,1,1),[1,nFPar,nModel]);
persistenceAllCentered = persistenceAll - offset;

%This is for getting the relative persistence change table.  Had to copy it
%over by hand, sucked.
meanPersistence = squeeze(mean(persistenceAll,1));
relPersistenceChange = (repmat(meanPersistence(1,:),[2,1]) - meanPersistence([6,nFPar],:))./repmat(meanPersistence(1,:),[2,1]);

CVBiomass = stdBiomass./meanBiomass;
log10CVBiomass = log10(CVBiomass);
log10CVBiomass((stdBiomass==0)|(meanBiomass==0)) = nan;

meanCVBiomass = permute(squeeze(mean(log10CVBiomass ,'omitnan')),[3 1 2]);

meanCVBiomassFree = permute(squeeze(...
    sum(log10CVBiomass.*freeSp ,'omitnan')./sum(finalBiomass>0&freeSp))...
                            ,[3 1 2]);
                        
meanCVBiomassPara  = permute(squeeze(...
    sum(log10CVBiomass.*paraSp ,'omitnan')./sum(finalBiomass>0&paraSp))...
                            ,[3 1 2]);
                        
meanCVBiomassBasal  = permute(squeeze(...
    sum(log10CVBiomass.*basalSp,'omitnan')./sum(finalBiomass>0&basalSp))...
                            ,[3 1 2]);


sumMeanBiomass = permute(squeeze(sum(meanBiomass)),[3 1 2]);
sumMeanParaBiomass = permute(squeeze(sum(meanBiomass.*paraSp)),[3 1 2]);
sumMeanBasalBiomass = permute(squeeze(sum(meanBiomass.*basalSp)),[3 1 2]);
sumMeanFreeBiomass = permute(squeeze(sum(meanBiomass.*(~basalSp))),[3 1 2]);
fracParaMeanFreeBiomass = sumMeanParaBiomass./sumMeanFreeBiomass;

meanParaBiomass = permute(squeeze(sum(meanBiomass.*paraSp)),[3 1 2]);
meanTotalBiomass = permute(squeeze(sum(meanBiomass)),[3 1 2]);

fracParaBiomass = meanParaBiomass./meanTotalBiomass;
fracFinalPara = permute(squeeze(sum((meanBiomass&paraSp)./repmat(sum((meanBiomass>0)&~basalSp),[40,1,1,1]))),[3 1 2]);
fracFinalFree = permute(squeeze(sum((meanBiomass&freeSp)./repmat(sum((meanBiomass>0)&~basalSp),[40,1,1,1]))),[3 1 2]);
                        
                        
diffPersistenceAll = diff(persistenceAll,1,2);
diffPersistenceAll = [nan(nWeb,1,nModel) diffPersistenceAll];

%Actual fraction of consumers as parasites in each web
fParVar = permute(repmat(nPara./(nFree+nPara),[1,nModel,1]),[3 1 2]);
fPar40 =  permute(repmat(nPara./40,[1,nModel,1]),[3 1 2]);
webVar = repmat((1:100)',[1,nFPar,16]);

fParNom = repmat(fParAll,100,1,16);
freeVar = permute(repmat(models(:,1)-1,[1,100,nFPar]),[2 3 1]);
paraVar = permute(repmat(models(:,2)-1,[1,100,nFPar]),[2 3 1]);
refVar = permute(repmat(models(:,3)-1,[1,100,nFPar]),[2 3 1]);
conVar = permute(repmat(models(:,4)-1,[1,100,nFPar]),[2 3 1]);
persistenceAll(persistenceAll==0) = nan;

%{
groupingAll = reshape(1 + freeVar + 2*paraVar + 4*refVar + 8*conVar,[],1);
offset(:,1,:) = [];
offset = reshape(offset,[],1);

varAll = [reshape(cat(4,persistenceAll,persistenceAllCentered,...
                     round((1-persistenceAll)*40),...
                     round((persistenceAllCentered + 1)*40),...
                     ... 
                     fParVar,fParVar-mean(fParVar(:,2:nFPar,:),2),fParNom,...
                     fPar40,...
                     freeVar,paraVar,refVar,conVar,...
                     webVar)...
                 ,[],13),groupingAll];

varAll(isinf(varAll)) = nan;
varAll(fParVar==0,:) = [];

varAllTab = array2table(varAll,...
    'VariableNames',{'Persistence',...                     1
                     'Persistence_Centered',...            2
                     'nExtct',...                          3
                     'nExtct_Centered',...                 4
                     'Parasite_Fraction_notCentered',...   5
                     'Parasite_Fraction',...               6
                     'NomFraction',...                     7
                     'Fraction_40',...                     8 
                     'Small_Free',...                      9
                     'Big_Para',...                        10 
                     'Refuge',...                          11 
                     'Concomittant',...                    12 
                     'Web',...                             13
                     'modelID'});                       ...14    
varAllTab.Web = categorical(varAllTab.Web);
varAllTab.modelID = categorical(varAllTab.modelID);
dataUsed = [1 6];


modelSpec0 = 'Persistence ~ Parasite_Fraction';
modelCSpec0 = 'Persistence_Centered ~ Parasite_Fraction';
modelNSpec0 = 'nExtct ~ Parasite_Fraction';
addModels = ' + Small_Free + Big_Para + Refuge + Concomittant';
addQuad = ' + Parasite_Fraction^2';
addInter = ' + Small_Free:Big_Para + Refuge:Concomittant';
addWebs = ' + (1|Web)';
addParaCondWeb = ' + (Parasite_Fraction|Web)';
addFracRefCInter = ' + Parasite_Fraction^2*Refuge';
addFracRefInter = ' + Parasite_Fraction^2*Refuge';
addFracRefNInter = ' + Parasite_Fraction^2*Refuge';



%Splitting up training and testing data.  Fair comparison of predictive
%power?  
fTrain = 0.7;
nTrain = round(fTrain*nWeb);
trainData = randsample(nWeb,nTrain);
trainBin = false(nWeb,1);
trainBin(trainData) = true;

%{
n_ii = 160;


modelBin = varAll(varAll(:,13)==1,9:12);
grouping = sum(modelBin.*[1 2 4 8],2);
modelAllBin = repmat(modelBin,nTrain,1);
groupingAll = sum(modelAllBin.*[1 2 4 8],2);
webId = reshape(repmat(1:nTrain,160,1),[],1);

hatMxAll = cell(1,nTrain);
YHatAll = cell(1,nTrain);
residAll = cell(1,nTrain);

count = 0;
for ii =  trainData'
    count = count+1;
        
    data_ii = varAll(varAll(:,13)==ii,[1,6]);
    Y = data_ii(:,1);
    fPar_ii = data_ii(:,2);
    
    X = [ones(n_ii,1)...
         ,fPar_ii...     
         ,modelBin(:,1)... All that should matter here is free liver
         ,fPar_ii.*modelBin ... Posibly matters
         ...,fPar_ii.^2 ... After peaking at data, this looks good. 
         ...,fPar_ii.^2.*modelBin ... Also might be important.
         ... model interactions?
         ,fPar_ii.*modelBin(:,1).*modelBin(:,2) ...
         ,fPar_ii.*modelBin(:,1).*[modelBin(:,1), modelBin(:,2)].*modelBin(:,3) ...
         ]; ... Maybe, but can I justify?
         
       [n,p] = size(X);
    
    if count == 1
        allCoef = zeros(p,nTrain);
        anovaAll = zeros(nTrain,4);
        nAll = zeros(nTrain,1);
        residAll = nan(n,nTrain);
        YHatAll = nan(n,nTrain);
        YAll =  nan(n,nTrain);
        XAll = nan(n,p,nTrain);
    end
    
    %There are a few missing observations
    missing_no = isnan(Y);
    X(missing_no,:) = [];
    Y(missing_no) = [];
    
    [n,p] = size(X);
    
    
    XtX = X'*X;
    b = XtX\(X'*Y);
    H = X*(XtX\X');
    Yhat = H*Y;
    resid = Y-Yhat;
    SSR = Y'*(H-ones(n)/n)*Y;
    SSE = Y'*(eye(n)-H)*Y;
    SST = Y'*(eye(n)-ones(n)/n)*Y;
    Rsq = 1- (SSE/(n-p))/(SST/(n-1));
    
    anovaAll(count,:) = [SSR SSE SST Rsq];
    nAll(count) = n;
    allCoef(:,count) = b;
    hatMxAll{count} = H;
    
    YHatAll(~missing_no,count) = Yhat;
    YAll(~missing_no,count) = Y;
    XAll(~missing_no,:,count) = X;
    residAll(~missing_no,count) = resid;
    
end

residAllVec = reshape(residAll,[],1);
YHatVec = reshape(YHatAll,[],1);
scatter(groupingAll,residAllVec,'o');
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%  MAKING STUFF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
S = 40;
C = .15;

close all




SList = 1:S;

%The abundance level web directory
dataDir = sprintf('/Volumes/Buster/Research/Parasitism/Dynamics/JulySimulations');


fParAll = linspace(0,0.5,11);

codeDir = pwd;
startWeb = 1;
endWeb = 100;
nWeb = endWeb-startWeb+1;

%Considered doing this, just want to save the code snippet here.  counts
%the number of files in the directory with the given filter.  (look up
%doc...)
%nWeb = numel(dir([dataDirSC,'/web*']))


%Assuming that model 0000 is null (it is and always should be!), we get the 
%total number of parasite fractions (nCol) here (automatically!):
model = models(1,:);


%Getting the webs; filtered by abundance then connectance
dataDirS = sprintf('%s/S%03u',dataDir,S);

%The abundance,connectance level web directory
dataDirSC = sprintf('%s/C%03u',dataDirS,round(100*C));


ZFrees = false(nModels,2);
ZParas = false(nModels,2);
freeLivings = false(nModels,2);
concs = false(nModels,2);

ZFrees(:,1) = models(:,1) == 1; %Small Free
ZFrees(:,2) = models(:,1) == 2; %big Free
ZParas(:,1) = models(:,2) == 1; %Big Para
ZParas(:,2) = models(:,2) == 2; % SMall Para
freeLivings(:,1) = models(:,3) == 1; %Don't include free LIving
freeLivings(:,2) = models(:,3) == 2; %INclude Free Living
concs(:,1) = models(:,4) == 1; %Don't include concomittant
concs(:,2) = models(:,4) == 2; %INclude concomittant
AXW = fullfact([2,2]);



CI = true;
sigLevel = .05;
nullChoice = ZFrees(:,2)+1;

xmin = min(fParAll)-.05;
xmax = max(fParAll)+.05;
ymin = -0.05;
ymax = 1.05;

middlePercentage = .5;
lowerQuantile = (1-middlePercentage)/2;
upperQuantile = 1 - lowerQuantile;

AXW = fullfact([2,2]);
ZParaColors = {'r','b'};
ZFreeMarks = {'o','x'};

marksMatrix = {{'k^','MarkerFaceColor','k'},{'x','Color',[.75,0,.75]};
    {'v','MarkerFaceColor',[0,.5,.5]}, {'+','Color',[0 0.75 0]}};

legendEntries = {'(Z_{Free},Z_{Para})=(10,1e-3)'...
    ,'(Z_{Free},Z_{Para})=(10,1e-4)'...
    ,'(Z_{Free},Z_{Para})=(100,1e-3)'...
    ,'(Z_{Free},Z_{Para})=(100,1e-4)'...
    };

param = struct('xVals',fParAll...
    ,'makePlot',true...
    ,'savePlot',false...
    ,'saveData',true...
    ,'statFcn',@mean...
    ,'statFcnOpt','omitnan'...
    ,'errBars',struct...
    ,'plotParam',struct...
    ,'dataParam',struct...
    ,'modelCodes',models...
    ,'dataAnalyzed','persistenceAll');

param.errBars = struct('errBarsFcn',@calculateConfidenceIntervals...
    ,'param',struct);

param.errBars.param = struct('alpha',.05);

param.errBars.flag = true;

param.plotParam = struct('marksMatrix',{marksMatrix} ...
    ,'legendEntries',{legendEntries} ...
    ,'title1',sprintf('No Free-Living Stage')...
    ,'title2',sprintf('Free-Living Stage')  ...
    ,'y2label1',sprintf('No Concomittant') ...
    ,'y2label2',sprintf('Concomittant') ...
    ,'ylabel',sprintf('Fraction of All Species %s'...
    ,param.dataAnalyzed) ...
    ,'xlabel',sprintf('Fraction of Consumers as Parasites') ...
    ,'plotRange',[xmin,xmax,ymin,ymax] ...
    ,'filename','AllSpeciesPersistenceFig'...
    ,'format','jpg');

param.dataParam =...
    struct('headers','x,y,yp,ym' ...
    ,'filenames',{{'AllSpecies00';...
    'AllSpecies0F';...
    'AllSpeciesC0';...
    'AllSpeciesCF';}}...
    ,'filePath','../../../Papers/ParasiteSwitching/data/'...
    );

param.regressionParam = ...
    struct('doRegression',true...
           ,'cmd',@fitlm...
           ,'saveReg',true...
           ,'regOpts',struct);

processExperimentOutput(persistenceAllCentered,param);
% 
% param.dataAnalyzed = 'diff-persistence';
% param.plotParam.ylabel = 'Sequential persistence difference';
% processExperimentOutput(diffPersistenceAll,param);

param.dataAnalyzed = 'net-extct-per-para';
param.plotParam.ylabel = 'Net new extinctions per parasite added';
processExperimentOutput(netExtctPerPara,param);

param.errBars.flag = false;
param.dataAnalyzed = 'frac-more-extct';
param.plotParam.ylabel = 'Fraction of runs with Non-negative persistence Change';
processExperimentOutput(netExtctPerPara<=0,param);

param.dataAnalyzed = 'frac-less-extct';
param.plotParam.ylabel = 'Fraction of runs with Persistence Decrease';
processExperimentOutput(netExtctPerPara>0,param);

param.errBars.flag = true;
param.dataAnalyzed = 'frac-cff';
param.plotParam.ylabel = 'Persistence of ff links';
processExperimentOutput(fracLinkChange(:,:,:,1),param);

param.dataAnalyzed = 'frac-cfp';
param.plotParam.ylabel = 'Persistence of fp links';
processExperimentOutput(fracLinkChange(:,:,:,2),param);

param.dataAnalyzed = 'frac-cpf';
param.plotParam.ylabel = 'Persistence of pf links';
processExperimentOutput(fracLinkChange(:,:,:,3),param);

param.dataAnalyzed = 'frac-cpp';
param.plotParam.ylabel = 'Persistence of pp links';
processExperimentOutput(fracLinkChange(:,:,:,4),param);

param.dataAnalyzed = 'free-persistence';
param.plotParam.ylabel = 'Persistence of Free Livers';
processExperimentOutput(freePersistence,param);

param.dataAnalyzed = 'para-persistence';
param.plotParam.ylabel = 'Persistence of Parasites';
processExperimentOutput(paraPersistence,param);

param.dataAnalyzed = 'basal-persistence';
param.plotParam.ylabel = 'Persistence of Basal Species';
processExperimentOutput(basalPersistence,param);

% param.dataAnalyzed = 'hamming-null';
% param.plotParam.ylabel = 'Hamming Distance from Null';
% processExperimentOutput(hammingNullAll,param);

%}
%Tried glme.  Better in the no random effects case, but inferior when
%constant random effects for webs are included.  Might need to come back
%and mention the results here, but I am hoping to avoid it since I'm not
%an expert here.
%{
glme0 = fitlme(varAllTab,modelNSpec0,...
                'Distribution','binomial',...
                'binomialSize',40);
glme0o = fitlme(varAllTab,modelNSpec0,...
                'Distribution','binomial',...
                'binomialSize',40,...
                'offset',offset);

modelNSpec = strcat(modelNSpec0,addModels);
glme1 = fitlme(varAllTab,modelNSpec,...
                'Distribution','binomial',...
                'binomialSize',40);
glme1o = fitlme(varAllTab,modelNSpec,...
                'Distribution','binomial',...
                'binomialSize',40,...
                'offset',offset);

modelNSpec = strcat(modelNSpec0,addModels,addQuad);
glme2 = fitlme(varAllTab,modelNSpec,...
                'Distribution','binomial',...
                'binomialSize',40);
glme2o = fitlme(varAllTab,modelNSpec,...
                'Distribution','binomial',...
                'binomialSize',40,...
                'offset',offset);

modelNSpec = strcat(modelNSpec0,addModels,addQuad,addInter);
glme3 = fitlme(varAllTab,modelNSpec,...
                'Distribution','binomial',...
                'binomialSize',40);
glme3o = fitlme(varAllTab,modelNSpec,...
                'Distribution','binomial',...
                'binomialSize',40,...
                'offset',offset);

modelNSpec = strcat(modelNSpec0,addModels,addQuad,addWebs);
glme4 = fitlme(varAllTab,modelNSpec,...
                'Distribution','binomial',...
                'binomialSize',40);
glme4o = fitlme(varAllTab,modelNSpec,...
                'Distribution','binomial',...
                'binomialSize',40,...
                'offset',offset);
            
modelNSpec = strcat(modelNSpec0,addModels,addWebs);
glme5 = fitlme(varAllTab,modelNSpec,...
                'Distribution','binomial',...
                'binomialSize',40);
glme5o = fitlme(varAllTab,modelNSpec,...
                'Distribution','binomial',...
                'binomialSize',40,...
                'offset',offset);
            
modelNSpec = strcat(modelNSpec0,addModels,addFracRefNInter,addWebs);
glme6 = fitlme(varAllTab,modelNSpec,...
                'Distribution','binomial',...
                'binomialSize',40);
glme6o = fitlme(varAllTab,modelNSpec,...
                'Distribution','binomial',...
                'binomialSize',40,...
                'offset',offset);
%}
%{
lme0 = fitlme(varAllTab,modelSpec0);
lme0c = fitlme(varAllTab,modelCSpec0);

modelSpec = strcat(modelSpec0,addModels);
modelCSpec = strcat(modelCSpec0,addModels);
lme1 = fitlme(varAllTab,modelSpec);
lme1c = fitlme(varAllTab,modelCSpec);

modelSpec = strcat(modelSpec0,addModels,addQuad);
modelCSpec = strcat(modelCSpec0,addModels,addQuad);
lme2 = fitlme(varAllTab,modelSpec);
lme2c = fitlme(varAllTab,modelCSpec);

modelSpec = strcat(modelSpec0,addModels,addQuad,addInter);
modelCSpec = strcat(modelCSpec0,addModels,addQuad,addInter);
lme3 = fitlme(varAllTab,modelSpec);
lme3c = fitlme(varAllTab,modelCSpec);

modelSpec = strcat(modelSpec0,addModels,addQuad,addWebs);
modelCSpec = strcat(modelCSpec0,addModels,addQuad,addWebs);
lme4 = fitlme(varAllTab,modelSpec);
lme4c = fitlme(varAllTab,modelCSpec);

modelSpec = strcat(modelSpec0,addModels,addWebs);
modelCSpec = strcat(modelCSpec0,addModels,addWebs);
lme5 = fitlme(varAllTab,modelSpec);
lme5c = fitlme(varAllTab,modelCSpec);

modelSpec = strcat(modelSpec0,addModels,addFracRefInter,addWebs);
modelCSpec = strcat(modelCSpec0,addModels,addFracRefCInter,addWebs);
lme6 = fitlme(varAllTab,modelSpec);
lme6c = fitlme(varAllTab,modelCSpec);

modelCSpec = 'Persistence ~ 1 + Parasite_Fraction^2*(Small_Free + Big_Para + Refuge + Concomittant) + (1|Web)';
lme7c = fitlme(varAllTab,modelCSpec);
modelSpec = 'Persistence ~ 1 + Parasite_Fraction^2*(Small_Free*Big_Para + Refuge*Concomittant) + (1|Web)';
lme7 = fitlme(varAllTab,modelSpec);
%}
%{
modelsHighLow = 2*models-3;

freeDiffs = persistenceAll(:,2:11,modelsHighLow(:,1)==1)-persistenceAll(:,2:11,(modelsHighLow(:,1)==-1));
freeDiffsVec = reshape(freeDiffs,[],1);

freeSmall = modelsHighLow(modelsHighLow(:,1)==-1,1);
paraBig = modelsHighLow(modelsHighLow(:,1)==-1,2);
ref = modelsHighLow(modelsHighLow(:,1)==-1,3);
con = modelsHighLow(modelsHighLow(:,1)==-1,4);

plotStylePara = cell(6,1);
plotStylePara(paraBig==-1) = deal({'b'});
plotStylePara(paraBig==1) = deal({'r'});

plotStyleRef =  cell(6,1);
plotStyleRef(ref==-1) = deal({'o'});
plotStyleRef(ref==1) = deal({'x'});


plotStyleCon =  cell(6,1);
plotStyleCon(con==-1) = deal({'-'});
plotStyleCon(con==1) = deal({'--'});


meanMx = squeeze(mean(freeDiffs,'omitnan'));
hold on
for ii = 1:8
   plotStyle = strcat(plotStylePara{ii},plotStyleRef{ii},plotStyleCon{ii});
   plot(linspace(0.05,0.5,10),meanMx(:,ii),plotStyle)
    
end
plot(linspace(0.05,0.5,10),zeros(1,10),'k-')
plot(linspace(0.05,0.5,10),ones(1,10).*linspace(-0.075,0.075,7)','k--')
hold off

n = sum(isfinite(freeDiffsVec));
tStar = tinv(1-.05/4/2,n-1);  

meanFreeDiff = mean(freeDiffsVec,'omitnan');
varFreeDiff = var(freeDiffsVec,0,'omitnan');
ME = tStar*sqrt(varFreeDiff/n);
freeDiffCI = meanFreeDiff+[-ME,ME];
tFree = meanFreeDiff/ME*tStar;

paraDiffs = persistenceAll(:,2:11,modelsHighLow(:,2)==1)-persistenceAll(:,2:11,modelsHighLow(:,2)==-1);
paraDiffsVec = reshape(paraDiffs,[],1);
meanParaDiff = mean(paraDiffsVec,'omitnan');
varParaDiff = var(paraDiffsVec,0,'omitnan');

ME = tStar*sqrt(varParaDiff/n);
paraDiffCI = meanParaDiff+[-ME,ME];
tPara = meanParaDiff/ME*tStar;

freeSmall = modelsHighLow(modelsHighLow(:,2)==-1,1);
paraBig = modelsHighLow(modelsHighLow(:,2)==-1,2);
ref = modelsHighLow(modelsHighLow(:,2)==-1,3);
con = modelsHighLow(modelsHighLow(:,2)==-1,4);

plotStyleFree = cell(6,1);
plotStyleFree(freeSmall==-1) = deal({'b'});
plotStyleFree(freeSmall==1) = deal({'r'});

plotStyleRef =  cell(6,1);
plotStyleRef(ref==-1) = deal({'o'});
plotStyleRef(ref==1) = deal({'x'});


plotStyleCon =  cell(6,1);
plotStyleCon(con==-1) = deal({'-'});
plotStyleCon(con==1) = deal({'--'});


meanMx = squeeze(mean(paraDiffs,'omitnan'));

hold on
for ii = 1:8
   plotStyle = strcat(plotStyleFree{ii},plotStyleRef{ii},plotStyleCon{ii});
   plot(linspace(0.05,0.5,10),meanMx(:,ii),plotStyle)
    
end
plot(linspace(0.05,0.5,10),zeros(1,10),'k-')
plot(linspace(0.05,0.5,10),ones(1,10).*linspace(-0.075,0.075,7)','k--')
hold off


refDiffsVec = reshape(persistenceAll(:,2:11,modelsHighLow(:,3)==1)-persistenceAll(:,2:11,modelsHighLow(:,3)==-1),[],1);
meanRefDiff = mean(refDiffsVec,'omitnan');
varRefDiff = var(refDiffsVec,0,'omitnan');
ME = tStar*sqrt(varRefDiff/n);
refDiffCI = meanRefDiff+[-ME,ME];
tRef = meanRefDiff/ME*tStar;

conDiffsVec = reshape(persistenceAll(:,2:11,modelsHighLow(:,4)==1)-persistenceAll(:,2:11,modelsHighLow(:,4)==-1),[],1);
meanConDiff = mean(conDiffsVec,'omitnan');
varConDiff = var(conDiffsVec,0,'omitnan');
ME = tStar*sqrt(varConDiff/n);
conDiffCI = meanConDiff+[-ME,ME];
tCon = meanConDiff/ME*tStar;
[freeDiffCI;paraDiffCI;refDiffCI;conDiffCI]
[tFree;tPara;tRef;tCon]
%}
%{
%do regression from .05-.50 parasites. that is, index 3 to 12.

%Let's try to reduce variability by looking at delta extinctions. so,
%get this by subtracting the new persistence from the extinction AFTER THE
%FIRST PERTURBATION!  This is a little odd and might make for akward
%interpretations but I think this is really what I want; should (hopefully 
%reduce variability in the data by controlling for the effect of the
%initial perturbation/extinction level..  Huge difference between web with
%0 and web with 1 parasite; treat that as separate, distinct difference.
%increasing the number of parasites is omething that is more reasonable to
%compare.

%gives change in extinction: negative means fewer extinctions, positive
%means more.
%This is relative change
deltaExtinctions = 40*(-persistenceAll(:,2,:) + persistenceAll(:,3:end,:));%./persistenceAll(:,1,:);
bigJumps = sum(squeeze(max(abs(deltaExtinctions),[],2)>10),2)>0;
bigJumps = false(size(bigJumps));
deltaOmitted = deltaExtinctions;
deltaOmitted(bigJumps,:,:) = [];



n = numel(deltaOmitted);
lmData4dArray = cat(4,deltaOmitted...    1
                     ,fParNom(~bigJumps,3:end,:)...  2
                     ,fParVar(~bigJumps,3:end,:)...  3
                     ,fPar40(~bigJumps,3:end,:)...   4 = nX
                    );

nX= 4;
lmData2dArray = reshape(lmData4dArray,n,nX);
varNames = {'D_Extct',...   1
            'x1'...         2
            ,'x2'...        3
            ,'x3'...        4
            }; 
lmDataTab = array2table(lmData2dArray,'VariableNames',varNames);
lmDataTab.x4 = categorical(lmDataTab.x1);
lmDataTab.free = reshape(freeVar(~bigJumps,3:end,:)==1,[],1);
lmDataTab.para = reshape(paraVar(~bigJumps,3:end,:)==1,[],1);
lmDataTab.ref= reshape(refVar(~bigJumps,3:end,:)==1,[],1);
lmDataTab.con = reshape(conVar(~bigJumps,3:end,:)==1,[],1);

lmData2dArray = table2array(lmDataTab(:,[1:4,6:end]));
predMatrix = [ones(length(lmData2dArray),1),lmData2dArray];
predNames =  lmDataTab.Properties.VariableNames(2:end);
X1 = x2fx(lmData2dArray(:,[2,5:8]),'interactions',2:5);
X = x2fx(lmData2dArray(:,[3,5:8]),'interactions',2:5);
X = X(:,2:end);
%variable names:
subsetVarNames = {'constant','fracPar','free','para','ref','con','free:fracPar','para:fracPar','ref:fracPar','con:fracPar','free:para','free:ref','free:con','para:ref','para:con','ref:con'};
X3 = x2fx(lmData2dArray(:,[4,5:8]),'interactions',2:5);
Y = lmData2dArray(:,1);

%Borrowing from Loren's art of matlab to find best subsets:

nModels = 2^15-1;
indx = dec2bin(1:nModels)=='1';
best5Radj = zeros(5,15,2);
best5BIC = inf(5,15,2);
for ii = 1:nModels   
    ii
    lm = fitlm(X,Y,'linear','PredictorVars',indx(ii,:));
    
    nVar = sum(indx(ii,:));
    Rsqs = sortrows([squeeze(best5Radj(:,nVar,:));[lm.Rsquared.Adjusted ii]],-1);
    BICs = sortrows([squeeze(best5BIC(:,nVar,:));[lm.ModelCriterion.BIC,ii]]);
    
    best5Radj(:,nVar,:) = Rsqs(1:5,:);
    best5BIC(:,nVar,:) = BICs(1:5,:);
    
end



%So now we try building models.  
lm = struct('fNom1',struct...
           ,'fVar1',struct...
           ,'fTot1',struct...
           ,'fNom2',struct...
           ,'fVar2',struct...
           ,'fTot2',struct);

%%%%%%%%%%%%%%Using nominal parasite fraction, 
%%%%%% linear model
stepwiselmOpts = {'PredictorVars',predNames([1 5:8])...
                  ,'Criterion','BIC'...
                  ,'ResponseVar','D_Extct'...
                  ,'upper','linear'};
lm.fNom1 = struct('bkwd',stepwiselm(lmDataTab,'D_Extct~(x1*free*para*ref*con)'...
                              ,stepwiselmOpts{:}...
                             )...
                 ,'frwd',stepwiselm(lmDataTab,'constant'...
                             ,stepwiselmOpts{:}...
                             )...
                );
%%%%%% quadratic model
stepwiselmOpts = {'PredictorVars',predNames([1 5:8])...
                  ,'Criterion','BIC'...
                  ,'ResponseVar','D_Extct'...
                  ,'upper','quadratic'};
lm.fNom2 = struct('bkwd',stepwiselm(lmDataTab,'D_Extct~(x1^2*free*para*ref*con)'...
                              ,stepwiselmOpts{:}...
                             )...
                 ,'frwd',stepwiselm(lmDataTab,'constant'...
                             ,stepwiselmOpts{:}...
                             )...
                );
%%%%%%%%%%%%%%Using actual parasite fraction (consumers, 
%%%%%% linear model            
stepwiselmOpts = {'PredictorVars',predNames([2 5:8])...
                  ,'Criterion','BIC'...
                  ,'ResponseVar','D_Extct'...
                  ,'upper','linear'};
lm.fVar1 = struct('bkwd',stepwiselm(lmDataTab,'D_Extct~(x2*free*para*ref*con)'...
                              ,stepwiselmOpts{:}...
                             )...
                 ,'frwd',stepwiselm(lmDataTab,'constant'...
                             ,stepwiselmOpts{:}...
                             )...
                );
%%%%%% quadratic model     
stepwiselmOpts = {'PredictorVars',predNames([2 5:8])...
                  ,'Criterion','BIC'...
                  ,'ResponseVar','D_Extct'...
                  ,'upper','quadratic'};
lm.fVar2 = struct('bkwd',stepwiselm(lmDataTab,'D_Extct~(x2^2*free*para*ref*con)'...
                              ,stepwiselmOpts{:}...
                             )...
                 ,'frwd',stepwiselm(lmDataTab,'constant'...
                             ,stepwiselmOpts{:}...
                             )...
                );
            
%%%%%%%%%%%%%%Using actual parasite fraction (all sp, 
%%%%%% linear model            
stepwiselmOpts = {'PredictorVars',predNames([3 5:8])...
                  ,'Criterion','BIC'...
                  ,'ResponseVar','D_Extct'...
                  ,'upper','linear'};
lm.fTot1 = struct('bkwd',stepwiselm(lmDataTab,'D_Extct~(x3*free*para*ref*con)'...
                              ,stepwiselmOpts{:}...
                             )...
                 ,'frwd',stepwiselm(lmDataTab,'constant'...
                             ,stepwiselmOpts{:}...
                             )...
                );
%%%%%% quadratic model     
stepwiselmOpts = {'PredictorVars',predNames([3 5:8])...
                  ,'Criterion','BIC'...
                  ,'ResponseVar','D_Extct'...
                  ,'upper','quadratic'};
lm.fTot2 = struct('bkwd',stepwiselm(lmDataTab,'D_Extct~(x3^2*free*para*ref*con)'...
                              ,stepwiselmOpts{:}...
                             )...
                 ,'frwd',stepwiselm(lmDataTab,'constant'...
                             ,stepwiselmOpts{:}...
                             )...
                );
lmSummary=table();

lmSummary.AdjRsquared = [lm.fNom1.bkwd.Rsquared.Adjusted;
                         lm.fNom1.frwd.Rsquared.Adjusted;
                         lm.fNom2.bkwd.Rsquared.Adjusted;
                         lm.fNom2.frwd.Rsquared.Adjusted;
                         lm.fVar1.bkwd.Rsquared.Adjusted;
                         lm.fVar1.frwd.Rsquared.Adjusted;
                         lm.fVar2.bkwd.Rsquared.Adjusted;
                         lm.fVar2.frwd.Rsquared.Adjusted;
                         lm.fTot1.bkwd.Rsquared.Adjusted;
                         lm.fTot1.frwd.Rsquared.Adjusted;
                         lm.fTot2.bkwd.Rsquared.Adjusted;
                         lm.fTot2.frwd.Rsquared.Adjusted];
                     
lmSummary.Predictor = {'fracNom','fracNom','fracNom','fracNom'...
                      ,'fracVar','fracVar','fracVar','fracVar'...
                      ,'fracTot','fracTot','fracTot','fracTot'}';
                  
lmSummary.HighestTerm = {'linear','linear','quadratic','quadratic'...
                        ,'linear','linear','quadratic','quadratic'...
                        ,'linear','linear','quadratic','quadratic'}';
                    
lmSummary.Selection = {'Backward','Forward','Backward','Forward'...
                      ,'Backward','Forward','Backward','Forward'...
                      ,'Backward','Forward','Backward','Forward'}';
                  
lmSummary.NumberofTerms = [lm.fNom1.bkwd.NumCoefficients;
                           lm.fNom1.frwd.NumCoefficients;
                           lm.fNom2.bkwd.NumCoefficients;
                           lm.fNom2.frwd.NumCoefficients;
                           lm.fVar1.bkwd.NumCoefficients;
                           lm.fVar1.frwd.NumCoefficients;
                           lm.fVar2.bkwd.NumCoefficients;
                           lm.fVar2.frwd.NumCoefficients;
                           lm.fTot1.bkwd.NumCoefficients;
                           lm.fTot1.frwd.NumCoefficients;
                           lm.fTot2.bkwd.NumCoefficients;
                           lm.fTot2.frwd.NumCoefficients];
                       
lmSummary.Criterion = [lm.fNom1.bkwd.ModelCriterion.BIC;
                       lm.fNom1.frwd.ModelCriterion.BIC;
                       lm.fNom2.bkwd.ModelCriterion.BIC;
                       lm.fNom2.frwd.ModelCriterion.BIC;
                       lm.fVar1.bkwd.ModelCriterion.BIC;
                       lm.fVar1.frwd.ModelCriterion.BIC;
                       lm.fVar2.bkwd.ModelCriterion.BIC;
                       lm.fVar2.frwd.ModelCriterion.BIC;
                       lm.fTot1.bkwd.ModelCriterion.BIC;
                       lm.fTot1.frwd.ModelCriterion.BIC;
                       lm.fTot2.bkwd.ModelCriterion.BIC;
                       lm.fTot2.frwd.ModelCriterion.BIC];

lmFinal = fitlm(lmDataTab,'D_Extct ~ x2*(free+ref) + para');
% effect of adding .025 parasites (i.e. 1 para): (Perturbtion Idea)
%web26 failed at .025 parasites for one of the models, so just delete it.  1%
%ilure not bad here?  Also interesting to regress perturbations against
%#basal species or #free livers, or tl of parasite, etc. (or all?)
k = 4;


omitted = isnan(sum(persistenceAll(:,2,:),3));
n = 100-sum(omitted);
%sumPerturb is the number of new extinctions after adding .05 prasites.
perturb = 40*(persistenceAll(~omitted,1,:)-persistenceAll(~omitted,2,:));
sumPerturb = squeeze(sum(perturb));
varPerturb = squeeze(var(perturb));
MSE = sum(varPerturb*(n-1))/(16*n);
%The contrasts:
mainContrasts = models*2-3;
%
%Compute the main efects:
mainEffects = sum(mainContrasts.*sumPerturb)/(2^(k-1)*n);
tMainEffects = mainEffects./(MSE/(2^(k-2)*n));

interaction2Idx = nchoosek([1 2 3 4],2);
interaction2Contrasts = mainContrasts(:,interaction2Idx(:,1)).*mainContrasts(:,interaction2Idx(:,2));
interaction2Effects = sum(interaction2Contrasts.*sumPerturb)/(2^(k-1)*n);
tInteraction2Effects = interaction2Effects./(MSE/(2^(k-2)*n));

interaction3Idx = nchoosek([1 2 3 4],3);
interaction3Contrasts = mainContrasts(:,interaction3Idx(:,1)).*mainContrasts(:,interaction3Idx(:,2)).*mainContrasts(:,interaction3Idx(:,3));
interaction3Effects = sum(interaction3Contrasts.*sumPerturb)/(2^(k-1)*n);
tInteraction3Effects = interaction3Effects./(MSE/(2^(k-2)*n));

interaction4Contrasts = prod(mainContrasts,2);
interaction4Effects = sum(interaction4Contrasts.*sumPerturb)/(2^(k-1)*n);
tInteraction4Effects = interaction4Effects./(MSE/(2^(k-2)*n));
%Thees are kind of tricky to interpret; but here goes:
%The first number is the number of additional extinctions from using larger
%free-liver body sizes after adding .05 free-livers compared to the web
%with no parasites.  negative indicates fewer extinctions.

% Then, do ols on data for .1-.5
z = [reshape(freeVar(~bigJumps,1,:)==1,[],1),...
reshape(paraVar(~bigJumps,1,:)==1,[],1),...
reshape(refVar(~bigJumps,1,:)==1,[],1),...
reshape(conVar(~bigJumps,1,:)==1,[],1)];
%}
%}

