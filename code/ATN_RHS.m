function [dB] = ATN_RHS(~,B,includeFraction,includeConcomittant,S,h,K,x,r,para,yij,yijTroList,yijParList,wijmx,wijTro,wijPar,eij,eijTroList,eijParList,res,con,paraLinkPara,paraLinkHost,trophicLinkRes,trophicLinkCon,paraCon,B0h,one,phi,basal)


if includeFraction == 0
    %{
    %This is a tricky thing. We want to avoid the complex numbers
    %that can arise from taking a weird root (1.2 in f.p.
    %representation is a rational number with even
    %denominator,which causes imaginary numbers & mucks
    %intermediate calculations of candidate time steps in ode45).
    %One option is to express 1.2 properly as a rational number.
    %This is faster and is equivalent,since we should theoretically
    %always get a positive number, anyway.
    %}
    Bh = abs(B).^(h);

    xB = x.*B;
    dens = B0h + wijmx'*Bh;
    yF = Bh(res)./(dens(con)).*yij.*xB(con);

    assim = (one'*sparse(res,con,yF,S,S))';
    loss = (one'*sparse(res,con,yF./eij,S,S)')';
    
    dB = basal.*r.*B.*(1-sum(B.*basal)/K);  %Basal logistic growth
    dB = dB - xB;                           %Metabolic losses
    dB = dB + assim;                        %Assimilation
    dB = dB - loss;                         %Loss due to predation.
    
    if includeConcomittant == 1
        
        %total losses due to trophic consumption of host(i.e.
        %non-parasitic consumer)

        totLossHost = loss - (one'*sparse(paraLinkPara,paraLinkHost,yF(paraCon)./eij(paraCon),S,S))';
        totConPara = assim.*para;
        
        %this is so horrible!
        fhpVec = yF(paraCon)./totConPara(paraLinkPara).*B(paraLinkPara)./B(paraLinkHost);
        fhp = sparse(paraLinkHost,paraLinkPara,fhpVec,S,S);

        fhp(isnan(fhp)) = 0;
        concomittantModifier = fhp'*totLossHost;
        
        dB = dB - concomittantModifier.*para;
    end
    
else
    %This model separates out trophic and parasitic links for
    %parasites. It probably takes about twice as long
    Bphi = B.*phi;
    BPhih = abs(Bphi).^(h);
    xBp = x.*Bphi;
    xB1p = x.*(B-Bphi);
    %need functional responses for each type of interaction:
    %trophic or parasitic.
    trophDen = B0h + wijTro'*BPhih;
    paraDen = B0h + wijPar'*BPhih;

    xyFTro = BPhih(trophicLinkRes)./(trophDen(trophicLinkCon)).*yijTroList.*xBp(trophicLinkCon);
    xyFPar = BPhih(paraLinkHost)./(paraDen(paraLinkPara)).*yijParList.*xB1p(paraLinkPara);
    
    assimTro = (one'*sparse(trophicLinkRes,trophicLinkCon,xyFTro,S,S))';
    assimPar = (one'*sparse(paraLinkHost,paraLinkPara,xyFPar,S,S))';
    
    lossTro = (one'*sparse(trophicLinkRes,trophicLinkCon,xyFTro./eijTroList,S,S)')';
    lossPar = (one'*sparse(paraLinkHost,paraLinkPara,xyFPar./eijParList,S,S)')';
    
    dB = r.*B.*(1-basal'*B/K)...    Basal logistic growth
        - B.*x...                   Metabolic losses
        +assimTro...                trophic consumption
        +assimPar...                parasitic consumption
        -lossTro...                 Loss due to predation.
        -lossPar;                   %Loss due to parasitism.
    
    
    if includeConcomittant == 1
        
        %Total losses due to trophic interactions for each host is lossTro.
        
        fhpVec = xyFPar./assimPar(paraLinkPara).*(B(paraLinkPara)-Bphi(paraLinkPara))./Bphi(paraLinkHost);
        fhp = sparse(paraLinkHost,paraLinkPara,fhpVec,S,S);
        fhp(isnan(fhp)) = 0;
        concomittantModifier = fhp'*lossTro;
        
        %Subtract it off.
        dB = dB - concomittantModifier;
    end
end


end
