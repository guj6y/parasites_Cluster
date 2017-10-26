%It isn't at all clear what you load here. So:
%out.mat contains extct and TS. 
%  extct is a cell array with nSims cells. Each cell contains a 41xn array. (1,:) slices the event times (including start and end). (2:end,ii) slices the biomasses of all species at the ii-th extinction event.
%  TS is a 3-d array of size 40x1000xnSims. This is the last 1000 time steps of the biomasses for all species. This should be useful for analyzing the equilibrium state (check for oscillations, fft,etc.)
%
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
load '../out.mat'
load '../setup_parameters.mat'

%The goal of this code is to convert the contents of out.mat to 4d arrays that my old code is capable of plotting. The 4d arrays need to be 40x14x16x100 (or some permutation of that...)
finalBiomasses = zeros(40,14,16,100);
meanBiomasses = zeros(40,14,16,100);
stdBiomasses = zeros(40,14,16,100);

nSims = numel(simParams);

fParAll = [0 0.025 0.05 0.1 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.75 1.0];

%Plotting code assumes that the models go in this stupid order:
models = circshift(fullfact([2,2,2,2])-1,[0 1]);

modelNo = bin2dec(num2str(models));

for ii = 1:nSims
    nWeb = simParams{ii}.web;
    model_ii = bin2dec(num2str([(simParams{ii}.kFree == 2) (simParams{ii}.kPara == -4) simParams{ii}.modelCode])); 
    modelIdx = find(modelNo == model_ii);
    fParIdx = find(simParams{ii}.fPar == fParAll);

    simParams{ii}.d2Slice = fParIdx;
    simParams{ii}.d3Slice = modelIdx;
    simParams{ii}.d4Slice = ii;
    
    meanBiomasses(:,fParIdx,modelIdx,nWeb) = mean(TS(:,:,ii),2);
    finalBiomasses(:,fParIdx,modelIdx,nWeb) = TS(:,end,ii);
    stdBiomasses(:,fParIdx,modelIdx,nWeb) = std(TS(:,:,ii),0,2);
end


