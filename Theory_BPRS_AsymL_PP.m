function [OP] = Theory_BPRS_AsymL_PP(KK,Nk,LL,LdB,yS,xPU,yPU,xPB,yPB,rho,alpha,beta,eta,PdB,IdB,Rth)  
% Best path relay selection
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
for kk = 1 : KK
    Hop = 0;
    for ii = 0 : LL-1
            Hop = Hop + fvartheta(LL,ii)/2.*...
                 fchi(LD_Dk,gmTH,OMG_Ik(ii,rho,LD_Ik(kk)),Ip);
    end   
    SSum = SSum + 1 - Hop;
end
% OP = ones(1,length(PP)).*(1 - SSum).^(Nk^(KK-1));
OP = ones(1,length(PP)).*(SSum).^(Nk^(KK-1));

semilogy(PdB,OP,'k-'); grid on;hold on; 
end




