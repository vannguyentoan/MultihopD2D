function [OP] = Theory_RS_PP(KK,MM,LL,LdB,yS,xPU,yPU,xPB,yPB,rho,alpha,beta,eta,PdB,IdB,Rth)
% hop by hop RS
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

% % OMG_I = (1+ii-ii*rho^2)*LD_Ik(kk)/(ii+1);
% function Xi = fXi(xx,yy,zz,kappa,PP,LD)
%     Xi = (xx/kappa./PP./yy/LD).^(zz/2).*besselk(zz,2*sqrt(xx/kappa./PP./yy/LD));
% end
% function Lambda = fLambda(xx,zz,kappa,PP,LD_D,LD_E,Ip,OMG_I)
%     Lambda = (xx/kappa./PP./LD_D/LD_E + Ip./kappa./PP./OMG_I/LD_E).^(zz/2).*besselk(zz,2.*sqrt(xx/kappa./PP./LD_D/LD_E + Ip./kappa./PP./OMG_I/LD_E));
% end
% %%% OMG_I
% function OMG = OMG_Ik(ii,rho,LD_Ik)
%     OMG = (1+ii-ii*rho^2)*LD_Ik/(ii+1);
% end
% %%% varpi
% function varpi = fvarpi(Ip,kappa,PP,OMG_Ik,LD_Ek)
%     varpi = Ip./kappa./PP./OMG_Ik/LD_Ek;
% end
% %%% varpi
% function vartheta = fvartheta(xx,yy)
%     vartheta = nchoosek(xx-1,yy)*2*xx*(-1)^yy/(yy+1);
% end
%%% psi
% function psi = fpsi(xx,yy)
%     psi = nchoosek(xx,yy)*(-1)^(yy-1);
% end
%
Prdt = 1;
for kk = 1 : KK
    % calculate for OP      
    Sum1Sub1 = 2*fXi(gmTH,LD_Dk,MM,kappa,PP,LD_Ek(kk));

    Sum1Sub2 = 0;
    for ii = 0 : LL-1
             Sum1Sub2 = Sum1Sub2 + fvartheta(LL,ii).*fXi(Ip,OMG_Ik(ii,rho,LD_Ik(kk)),MM,kappa,PP,LD_Ek(kk));
    end
    Sum1Sub3 = 0;
    for ii = 0 : LL-1
         Sum1Sub3 = Sum1Sub3 + fvartheta(LL,ii).*fLambda(gmTH,MM,kappa,PP,LD_Dk,LD_Ek(kk),Ip,OMG_Ik(ii,rho,LD_Ik(kk)));
    end
    %
    Sum1 = 1 - 1/gamma(MM)*(Sum1Sub1 + Sum1Sub2 - Sum1Sub3);
     % Calculate for I2
     
     Sum2 = 0;
     for tt = 0 : MM-1
         for ii = 0 : LL-1
             Sum2 = Sum2 + fvartheta(LL,ii)/factorial(tt).*fvarpi(Ip,kappa,PP,OMG_Ik(ii,rho,LD_Ik(kk)),LD_Ek(kk)).*(fXi(Ip,OMG_Ik(ii,rho,LD_Ik(kk)),tt-1,kappa,PP,LD_Ek(kk))...
                 - fLambda(gmTH,tt-1,kappa,PP,LD_Dk,LD_Ek(kk),Ip,OMG_Ik(ii,rho,LD_Ik(kk))));
         end
     end
%     
    Prdt = Prdt.*(1 - Sum1 - Sum2);
end
OP = 1 - Prdt;
semilogy(PdB,OP,'k-'); grid on;hold on; 
end




