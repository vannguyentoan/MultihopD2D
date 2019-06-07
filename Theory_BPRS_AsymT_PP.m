function [OP] = Theory_BPRS_AsymT_PP(KK,MM,Nk,xPU,yPU,xPB,yPB,alpha,beta,eta,PdB,IdB,Rth)  
% Best path relay selection
%
Ip          = 10.^(IdB/10);
PP          = 10.^(PdB/10);
% KK          = length(NumRelay);
d_Ii     = zeros(1,KK); % khoang cach tu relay - PU (Interference link)
for aa = 1 : KK
    d_Ii(aa) = sqrt(((aa-1)/KK-xPU)^2 + yPU^2);
end
LD_Ik    = d_Ii.^(-beta);%% LdB = 30 dB, --> L(W) = 10^(LdB/10)=10^3;

d_Ek      = zeros(1,KK);  % khoang cach tu relay - Beacon (Energy link)
for aa = 1 : KK
    d_Ek(aa) = sqrt(((aa-1)/KK-xPB)^2 + yPB^2);
end
LD_Ek    = d_Ek.^(-beta); % Lambda = L/(d^beta/d0) : simple pathloss model

d_Dk      = 1/KK;      % khoang cach tu relay i_th - relay (i+1)_th (Data link)                  
LD_Dk    = d_Dk^(-beta);

kappa = KK*eta*alpha/(1-alpha);

gmTH = 2^(KK*Rth/(1-alpha)) - 1;

% OMG_I = (1+ii-ii*rho^2)*LD_Ik(kk)/(ii+1);

SSum = 0;
for kk = 1 : KK
    % calculate for first term       
    Hop = 1 - 2/gamma(MM).*fXi(gmTH,LD_Dk,MM,kappa,PP,LD_Ek(kk));        
    SSum = SSum + Hop;
end
OP = SSum.^(Nk^(KK-1));

OP1 = ones(1,length(Ip)).*OP; % for IdB
semilogy(PdB,OP,'k-'); grid on;hold on; 
end




