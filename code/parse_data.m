%It isn't at all clear what you load here. So:
%out.mat contains extct and TS. 
%  extct is a cell array with nSims cells. Each cell contains a 41xn array. (1,:) slices the event times (including start and end). (2:end,ii) slices the biomasses of all species at the ii-th extinction event.
%  TS is a 3-d array of size 40x1000xnSims. This is the last 1000 time steps of the biomasses for all species. This should be useful for analyzing the equilibrium state (check for oscillations, fft,etc.)
%
%simParams
%.mat contains the cell array simParams. each element of the cell array is a structure. This is the structure:
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
load '../out.mat'
load '../simParams.mat'
load '../metaSimData.mat'

%The goal of this code is to convert the contents of out.mat to complicated multi-dimensional arrays that I can plot.
finalBiomasses = nan(S,nWeb,numel(fParAll0),nFacts(1),nFacts(2),nFacts(3),nFacts(4));
meanBiomasses = nan(S,nWeb,numel(fParAll0),nFacts(1),nFacts(2),nFacts(3),nFacts(4));
stdBiomasses = nan(S,nWeb,numel(fParAll0),nFacts(1),nFacts(2),nFacts(3),nFacts(4));

nSims = numel(extcts);

fParAll = fParAll0;

for ii = 1:nSims
    
    webNo = simParams{ii}.web;
    fact1Level = simParams{ii}.kFree == kFrees;
    fact2Level = simParams{ii}.kPara == kParas;
    fact3Level = simParams{ii}.fracFree;
    fact4Level = simParams{ii}.fracPara;
    fParLevel = simParams{ii}.fPar == fParALll0;
 
    meanBiomasses(:,webNo,fParLevel,fact1Level,fact2Level,fact3Level,fact4Level) = mean(TS(:,:,ii),2);
    finalBiomasses(:,webNo,fParLevel,fact1Level,fact2Level,fact3Level,fact4Level) = TS(:,end,ii);
    stdBiomasses(:,webNo,fParLevel,fact1Level,fact2Level,fact3Level,fact4Level) = std(TS(:,:,ii),0,2);

end

save('../out.mat','meanBiomasses','finalBiomasses','stdBiomasses');
