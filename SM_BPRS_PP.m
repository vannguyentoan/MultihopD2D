function [OP] = SM_BPRS_PP(KK,MM,Nk,LL,LdB,yS,xPU,yPU,xPB,yPB,rho,alpha,beta,eta,PdB,IdB,Rth,bit_frame)    
Ip          = 10.^(IdB/10);
PP          = 10.^(PdB/10);
OP          = zeros(1,length(PdB));
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

% %%%%
% d_En      = zeros(1,KK);  % khoang cach tu relay - Beacon (Energy link)
% for aa = 1 : KK
%     d_En(aa) = sqrt(((aa-1)*10/KK-xPB)^2 + yPB^2);
% end
% LD_En    = d_En.^beta/(10^3);
% d_Dn      = 10/KK;      % khoang cach tu relay i_th - relay (i+1)_th (Data link)                  
% LD_Dn    = d_Dn^(beta)/(10^3);
% %%%%%%%

kappa = KK*eta*alpha/(1-alpha);

for bb = 1 : length(PdB)
    fprintf('Running %d per %d \n',bb,length(PdB));
    for bitnum   =  1 : bit_frame(bb)
        SNR_e2e = zeros(1,KK);
        CC_ii = zeros(1,Nk^(KK-1));
        for ik = 1 : Nk^(KK-1)
            for kk = 1 : KK
                gmEk_sumk = 0;
                for dd = 1 : MM
                    h_Ei = sqrt(LD_Ek(kk)/2)*(randn(1,1)+1i*randn(1,1));
                    gmEk_sumk  = gmEk_sumk + abs(h_Ei)^2;
                end
                %  L primary users (imperfect channel gain)
                h_Ii = zeros(1,LL);
                for ii = 1 : LL      
                  h_Ii(1,ii)       = sqrt(LD_Ik(kk)/2)*(randn(1,1)+1i*randn(1,1));
                end
                Eps    = sqrt(LD_Ik(kk)/2)*(randn(1,1)+1i*randn(1,1));
                h_I    = rho*max(h_Ii) + sqrt(1-rho^2)*Eps;   % imperfect channel
                gm_Ik  = abs(h_I)^2;

                P_EHk  = kappa*PP(bb)*gmEk_sumk; % Power EH relay
                P_Ik   = Ip/gm_Ik;      % Interference Constraint
                Pk     = min(P_EHk,P_Ik);
                  % Relay selection
                h_Di = sqrt(LD_Dk/2)*(randn(1,1)+1i*randn(1,1));                   
                SNR_e2e(kk) = Pk*abs(h_Di)^2;
            end
            CC_ii(1,ik)  = (1-alpha)/KK*log2(1 + min(SNR_e2e));
        end
        CC_ib = max(CC_ii);
        if (CC_ib < Rth)            
            OP(bb) = OP(bb) + 1;
        end
    end
end
OP  = OP./bit_frame;
semilogy(PdB,OP,'o'); grid on;hold on; 
end




