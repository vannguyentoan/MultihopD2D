function [OP] = Theory_DHRS_even_PP(KK,MM,Nk,LL,LdB,yS,xPU,yPU,xPB,yPB,rho,alpha,beta,eta,PdB,IdB,Rth)  
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

Prdt = 1;
for vv = 1 : (KK-1)/2+1
    kk = 2*vv;
    if vv == 1
        % calculate for first hop      
        FSub1 = 0;
        for nn = 1 : Nk
             FSub1 = FSub1 + 2*fpsi(Nk,nn).*fXi(nn*gmTH,LD_Dk,MM,kappa,PP,LD_Ek(vv));
        end
        
        FSub2 = 0;
        for ii = 0 : LL-1
                 FSub2 = FSub2 + fvartheta(LL,ii).*fXi(Ip,OMG_Ik(ii,rho,LD_Ik(vv)),MM,kappa,PP,LD_Ek(vv));
        end
        
        FSub3 = 0;
        for ii = 0 : LL-1
             for nn = 1 : Nk
                 FSub3 = FSub3 + fpsi(Nk,nn).*fvartheta(LL,ii).*...
                     fLambda(nn*gmTH,MM,kappa,PP,LD_Dk,LD_Ek(vv),Ip,OMG_Ik(ii,rho,LD_Ik(vv)));
             end
        end
        T31 = 1 - 1/gamma(MM).*(FSub1 + FSub2 - FSub3);

         FSum2 = 0;
         for tt = 0 : MM-1
             for nn = 0: Nk
                 for ii = 0 : LL-1
                     FSum2 = FSum2 + fvartheta(LL,ii)*fpsi(Nk,nn)*(-1).*fvarpi(Ip,kappa,PP,OMG_Ik(ii,rho,LD_Ik(vv)),LD_Ek(vv))/factorial(tt)...
                     .*fLambda(nn*gmTH,tt-1,kappa,PP,LD_Dk,LD_Ek(vv),Ip,OMG_Ik(ii,rho,LD_Ik(vv)));
                 end
             end
         end
        T41 = FSum2;
        %     
        Hop = 1 - T31 - T41;
    else
        % Calculate for next dual hop
        % calculate for first term       
        DSub1 = 2.*fXi(gmTH,LD_Dk,MM,kappa,PP,LD_Ek(kk-2));

        DSub2 = 0;
        for ii = 0 : LL-1
                 DSub2 = DSub2 + fvartheta(LL,ii).*fXi(Ip,OMG_Ik(ii,rho,LD_Ik(kk-2)),MM,kappa,PP,LD_Ek(kk-2));
        end
        
        DSub3 = 0;
        for ii = 0 : LL-1
                DSub3 = DSub3 + fvartheta(LL,ii).*fLambda(gmTH,MM,kappa,PP,LD_Dk,LD_Ek(kk-2),Ip,OMG_Ik(ii,rho,LD_Ik(kk-2)));
        end
        
        T12km1 = 1 - 1/gamma(MM).*(DSub1 + DSub2 - DSub3);
        % Calculate for I2

         DSub4 = 0;
         for tt = 0 : MM-1
             for ii = 0 : LL-1
                 DSub4 = DSub4 + fvartheta(LL,ii).*fvarpi(Ip,kappa,PP,OMG_Ik(ii,rho,LD_Ik(kk-2)),LD_Ek(kk-2))./factorial(tt).*...
                     (fXi(Ip,OMG_Ik(ii,rho,LD_Ik(kk-2)),tt-1,kappa,PP,LD_Ek(kk-2))...
                     -fLambda(gmTH,tt-1,kappa,PP,LD_Dk,LD_Ek(kk-2),Ip,OMG_Ik(ii,rho,LD_Ik(kk-2))));
             end
         end 
       %     
        T22km1 = DSub4;

        % calculate for second term       
        DSSub1 = 0;
        for nn = 1 : Nk
             DSSub1 = DSSub1 + 2*fpsi(Nk,nn).*fXi(nn*gmTH,LD_Dk,MM,kappa,PP,LD_Ek(kk-1));
        end

        DSSub2 = 0;
        for ii = 0 : LL-1
                 DSSub2 = DSSub2 + fvartheta(LL,ii).*fXi(Ip,OMG_Ik(ii,rho,LD_Ik(kk-1)),MM,kappa,PP,LD_Ek(kk-1));
        end

        DSSub3 = 0;
        for ii = 0 : LL-1
             for nn = 1 : Nk
                 DSSub3 = DSSub3 + fpsi(Nk,nn)*fvartheta(LL,ii)...
                     .*fLambda(nn*gmTH,MM,kappa,PP,LD_Dk,LD_Ek(kk-1),Ip,OMG_Ik(ii,rho,LD_Ik(kk-1)));
             end
        end
        T32k = 1 - 1/gamma(MM).*(DSSub1 + DSSub2 - DSSub3);
         % Calculate for I2

         DSSum2 = 0;
         for tt = 0 : MM-1
             for nn = 0: Nk
                 for ii = 0 : LL-1
                     DSSum2 = DSSum2 + fvarpi(Ip,kappa,PP,OMG_Ik(ii,rho,LD_Ik(kk-1)),LD_Ek(kk-1))./factorial(tt)...
                     .* fvartheta(LL,ii)*fpsi(Nk,nn)*(-1).*fLambda(nn*gmTH,tt-1,kappa,PP,LD_Dk,LD_Ek(kk-1),Ip,OMG_Ik(ii,rho,LD_Ik(kk-1)));
                 end
             end
         end 
        T42k = DSSum2;
        
        Hop = 1 - (1 - (1 - T12km1 -T22km1).*(1 - T32k - T42k)).^Nk;
    end
    Prdt = Prdt.*Hop;
end
OP = 1 - Prdt;
semilogy(PdB,OP,'k-'); grid on;hold on; 
end




