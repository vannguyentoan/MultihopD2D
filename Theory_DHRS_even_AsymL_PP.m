function [OP] = Theory_DHRS_even_AsymL_PP(KK,Nk,LL,LdB,yS,xPU,yPU,xPB,yPB,rho,alpha,beta,eta,PdB,IdB,Rth)  
% number of cluster is even (chan)(K-1 = 2,4,6,8...)
% --> number of hop odd (K = 3,5,7,9...)
% --> divide into(K-1)/2 dual hop and first hop
%
Ip          = 10.^(IdB/10);
PP          = 10.^(PdB/10);
% KK          = length(NumRelay);
LmW         = 10^(LdB/10);
% KK          = length(NumRelay);
d_Ik     = zeros(1,KK);     % khoang cach tu relay - PU (Interference link)
for aa = 1 : KK
    d_Ik(aa) = sqrt(((aa-1)/KK-xPU)^2 + (yS-yPU)^2);
end
LD_Ik    = d_Ik.^(-beta)*LmW;   %% LdB = 30 dB, --> L(W) = 10^(LdB/10)=10^3;

d_Ek      = zeros(1,KK);    % khoang cach tu relay - Beacon (Energy link)
for aa = 1 : KK
    d_Ek(aa) = sqrt(((aa-1)/KK-xPB)^2 + (yS-yPB)^2);
end
LD_Ek    = d_Ek.^(-beta)*LmW;   % Lambda = L/(d^beta/d0) : simple pathloss model

d_Dk      = 30/KK;           % khoang cach tu relay i_th - relay (i+1)_th (Data link)                  
LD_Dk    = d_Dk^(-beta)*LmW;

kappa = KK*eta*alpha/(1-alpha);

gmTH = 2^(KK*Rth/(1-alpha)) - 1;

% OMG_I = (1+ii-ii*rho^2)*LD_Ik(kk)/(ii+1);

SSum = 0;
for vv = 1 : (KK-1)/2+1
    kk = 2*vv;
    if vv == 1
        % calculate for the first hop      
      
        FSub3 = 0;
        for ii = 0 : LL-1
             for nn = 1 : Nk
                 FSub3 = FSub3 + fpsi(Nk,nn)*fvartheta(LL,ii)/2.*fchi(LD_Dk,nn*gmTH,OMG_Ik(ii,rho,LD_Ik(vv)),Ip);
             end
        end
        %     
        Hop = 1 - FSub3;
    else
        % calculate for first term       

        DSub2 = 0;
        for ii = 0 : LL-1
                 DSub2 = DSub2 + fvartheta(LL,ii)/2.*fchi(LD_Dk,gmTH,OMG_Ik(ii,rho,LD_Ik(kk-2)),Ip);
        end

        FirstTerm = 1 - DSub2;

        % calculate for second term       

        DSSub3 = 0;
        for ii = 0 : LL-1
             for nn = 1 : Nk
                 DSSub3 = DSSub3 + fpsi(Nk,nn)*fvartheta(LL,ii)/2.*fchi(LD_Dk,nn*gmTH,OMG_Ik(ii,rho,LD_Ik(kk-1)),Ip);
             end
        end
    %     
        SecondTerm = 1 - DSSub3;
        Hop = (1 - (1-FirstTerm).*(1-SecondTerm)).^Nk;
    end
    SSum = SSum + Hop;
end
OP = ones(1,length(PP)).*SSum;
% OP = SSum;
semilogy(PdB,OP,'k-'); grid on;hold on; 
end




