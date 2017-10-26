fig = figure('Position', [0,0,1440,900],'Visible','off');
alpha = 0.05;
dataToPlotCode  = 'perAll';
plottingData = struct('perAll',persistenceAll...
                      ,'perBasal',persistenceBasal...
                      ,'perFree',persistenceFree...
                      ,'perPara',persistencePara...
                      ,'cvAll',meanCVBiomass...
                      ,'cvBasal',meanCVBiomassBasal...
                      ,'cvFree',meanCVBiomassFree...
                      ,'cvPara',meanCVBiomassPara...
                      ,'fracParaBio',fracParaMeanFreeBiomass...
                      ,'totalBio',sumMeanBiomass...
                      ... ,'slopes_e',slopes_e,...
                      ... ,'slopes',slopes,...
                      ... ,'rsquared',rSquares,...
                      ... ,'rsquared_e',rSquares_e...
                      );

ylabels = struct('perAll','Persistence of All Species',...
                'perBasal','Persistence of Basal Species',...
                'perFree','Persistence of Free-living Species',...
                'perPara','Persistence of Parasitic Species',...
                'cvAll','CV of all species',...
                'cvBasal','CV of basal species',...
                'cvFree','CV of free-living species',...
                'cvPara','CV of parasitic species',...
                'fracParaBio','Fraction of Consumer Biomass as Parasites',...
                'totalBio','Total Average Biomass',...
                'slopes_e','Average Slope of Abundance Body Size Curves (ln)',...
                'slopes','Average Slope of Abundance Body Size Curves (log_{10})',...
                'rsquared','Average r-squared of Abundance Body Size Curves (log_{10})',...
                'rsquared_e','Average r-squared of Abundance Body Size Curves (ln)');

dataToPlot = plottingData.(dataToPlotCode);
%dataToPlot = fracParaMeanFreeBiomass;
includeFStats = false;
subplotTitles = {'Free Liver BSR (Z_f)','Parasite BSR (Z_p)','Host Refuge','Concomittant Consumption'};
SS = zeros(2,2);
DF = SS;
MS = SS;
legendEntries = {'Z_f = 10','Z_f=100';
                 'Z_p = 10^{-3}','Z_p = 10^{-4}';
                 'Without Refuge','With Refuge';
                 'Without Concomittant Losses','With Concomittant Losses'};
xmax = -inf;
xmin = inf;
ymax = -inf;
ymin = inf;


printData = false;

for ii = 1:4
    subplot(2,2,ii)
    
    dataOff = dataToPlot(:,:,models(:,ii)==1);
    dataOn = dataToPlot(:,:,models(:,ii)==2);
    
    dataOff = reshape(permute(dataOff,[1 3 2]),[],nFPar);
    dataOn = reshape(permute(dataOn,[1 3 2]),[],nFPar);
    
    
    
    nOff = sum(isfinite(dataOff));
    nOn = sum(isfinite(dataOn));
    
    meanOff = mean(dataOff,'omitnan');
    meanOn = mean(dataOn,'omitnan');
    
    I = sum(~isnan(meanOn));
    
    %Anova (hard coded.. why not?)
    SS(1,1) = sum(sum((dataOff-meanOff).^2,'omitnan'),'omitnan');
    DF(1,1) = sum(nOff) - I;
    SS(2,1) = sum((meanOff-mean(dataOff(:))).^2.*nOff,'omitnan');
    DF(2,1)= I - 1;
    
    SS(1,2) = sum(sum((dataOn-meanOn).^2,'omitnan'),'omitnan');
    DF(1,2) = sum(nOn) - I;
    SS(2,2) = sum((meanOn-mean(dataOn(:))).^2.*nOn,'omitnan');
    DF(2,2)= I - 1;
    
    MS = SS./DF;
    F = MS(2,:)./MS(1,:);
    
    stdOff = std(dataOff,[],'omitnan');
    stdOn = std(dataOn,[],'omitnan');
    
    stdErrOff = stdOff./sqrt(nOff);
    stdErrOn = stdOn./sqrt(nOn);
    
    
    tCritOff = tinv(1-alpha/(2*2*I),nOff-1);
    tCritOn = tinv(1-alpha/(2*2*I),nOn-1);
    
    hold on
    errorbar(fParAll,meanOff,stdErrOff.*tCritOff,'^-');
    errorbar(fParAll,meanOn,stdErrOn.*tCritOn,'o-');
    hold off
    
    
    title(subplotTitles{ii});
    if includeFStats
    legend(sprintf('off (1), F=%.2e (P=%.3f)',F(1),fcdf(F(1),DF(2,1),DF(1,1),'upper')),...
        sprintf('on (2), F=%.2e (P=%.3f)',F(1),fcdf(F(2),DF(2,2),DF(1,2),'upper')))
    else
        legend(legendEntries{ii,:})
    end
    if mod(ii,2)==1
        ylabel(ylabels.(dataToPlotCode));
    end
    
    if ii >2
        xlabel('Fraction of free livers as parasites')
    end
    xl = xlim;
    yl = ylim;
    
    xmin = min(xmin,xl(1));
    xmax = max(xmax,xl(2));
    
    ymin = min(ymin,yl(1));
    ymax = max(ymax,yl(2));
    grid on
    nanRowsOn = isnan(meanOn)|isnan(stdErrOn);
    if printData 
    fidOn = fopen(sprintf('../../../Papers/ParasiteSwitching/data/%s-effect-on-subplot-%u',dataToPlotCode,ii),'w');
    fprintf(fidOn,'x,y,yp,ym\n');
    fprintf(fidOn,'%.5f,%.5f,%.5f,%.5f\n',...
        [fParAll(~nanRowsOn)',meanOn(~nanRowsOn)',stdErrOn(~nanRowsOn)',stdErrOn(~nanRowsOn)']');
    fclose(fidOn);
    
    nanRowsOff = isnan(meanOff)|isnan(stdErrOff);
    fidOff = fopen(sprintf('../../../Papers/ParasiteSwitching/data/%s-effect-off-subplot-%u',dataToPlotCode,ii),'w');
    fprintf(fidOff,'x,y,yp,ym\n');
    fprintf(fidOff,'%.5f,%.5f,%.5f,%.5f\n',...
        [fParAll(~nanRowsOff)',meanOff(~nanRowsOff)',stdErrOff(~nanRowsOff)',stdErrOff(~nanRowsOff)']');
    fclose(fidOff);
end
    %refline(1,0);
    
end

for ii = 1:4
    subplot(2,2,ii);
    axis([xmin xmax ymin ymax])
end
format = 'jpg';
figFilename = sprintf('../%s.%s',dataToPlotCode,format);
saveas(fig,figFilename,format)
