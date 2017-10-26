close all
paraNo = 1;

bsAll = cat(3,repmat(bsAll(:,:,1,:),1,1,4,1),...
    repmat(bsAll(:,:,3,:),1,1,4,1),...
    repmat(bsAll(:,:,2,:),1,1,4,1),...
    repmat(bsAll(:,:,4,:),1,1,4,1));


abundance = meanBiomass./bsAll;

bsIntegral = round(log10(bsAll));
bsIntegral_e = round(log(bsAll));

abundance = permute(abundance,[2 3 4 1]);
bsIntegral = permute(bsIntegral,[2 3 4 1]);
bsIntegral_e = permute(bsIntegral_e,[2 3 4 1]);

webNo = 1:100;

inters = zeros(nFPar,nModels,nWeb);
slopes = zeros(nFPar,nModels,nWeb);
rSquares = zeros(nFPar,nModels,nWeb);


inters_e = zeros(nFPar,nModels,nWeb);
slopes_e = zeros(nFPar,nModels,nWeb);
rSquares_e = zeros(nFPar,nModels,nWeb);

for ii = 1:nWeb
    ii
    for jj = 1:nModels
        for kk = 1:nFPar
            
            [m,n] = grpstats(squeeze(abundance(kk,jj,ii,:)),squeeze(bsIntegral(kk,jj,ii,:)),{'mean','numel'});
            [m_e,n_e] = grpstats(squeeze(abundance(kk,jj,ii,:)),squeeze(bsIntegral_e(kk,jj,ii,:)),{'mean','numel'});
            
            x = unique(bsIntegral(kk,jj,ii,:));
            x_e = unique(bsIntegral_e(kk,jj,ii,:)); 
            y = log10(m.*n);
            y_e = log(m_e.*n_e);
            
            y(~isfinite(y)) = nan;
            y_e(~isfinite(y_e)) = nan;
            
            lm = fitlm(x,y);
            lm_e = fitlm(x_e,y_e);
            
            inters(kk,jj,ii) = lm.Coefficients.Estimate(1);
            slopes(kk,jj,ii) = lm.Coefficients.Estimate(2);
            rSquares(kk,jj,ii) = lm.Rsquared.Ordinary;
            
            inters_e(kk,jj,ii) = lm_e.Coefficients.Estimate(1);
            slopes_e(kk,jj,ii) = lm_e.Coefficients.Estimate(2);
            rSquares_e(kk,jj,ii) = lm_e.Rsquared.Ordinary;
        end
    end
end

slopes = permute(slopes,[3 1 2]);
slopes_e = permute(slopes_e,[3 1 2]);

rSquares = permute(rSquares, [3 1 2]);
rSquares_e = permute(rSquares_e, [3 1 2]);

            
