function [OP] = Theory_RS_AsymT(KK,MM,xPU,yPU,xPB,yPB,alpha,beta,eta,PdB,IdB,Rth)
% hop by hop RS
Ip          = 10.^(IdB/10);
PP          = 10.^(PdB/10);
% KK          = length(NumRelay);
d_Ii     = zeros(1,KK); % khoang cach tu relay - PU (Interference link)
for aa = 1 : KK
    d_Ii(aa) = sqrt(((aa-1)*10/KK-xPU)^2 + yPU^2);
end
LD_Ik    = d_Ii.^(-beta)*(10^3);%% LdB = 30 dB, --> L(W) = 10^(LdB/10)=10^3;

d_Ek      = zeros(1,KK);  % khoang cach tu relay - Beacon (Energy link)
for aa = 1 : KK
    d_Ek(aa) = sqrt(((aa-1)*10/KK-xPB)^2 + yPB^2);
end
LD_Ek    = d_Ek.^(-beta)*(10^3); % Lambda = L/(d^beta/d0) : simple pathloss model

d_Dk      = 10/KK;      % khoang cach tu relay i_th - relay (i+1)_th (Data link)                  
LD_Dk    = d_Dk^(-beta)*(10^3);

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

Sum1 = 0;
for kk = 1 : KK
    % calculate for OP      
    Sum1Sub1 = 1 - 2/gamma(MM)*fXi(gmTH,LD_Dk,MM,kappa,PP,LD_Ek(kk));
    Sum1 = Sum1 + Sum1Sub1;
end
% OP_PdB = Sum1;
OP_IdB = ones(1,length(IdB)).*Sum1;
semilogy(IdB,OP_IdB,'k-'); grid on;hold on; 
end




