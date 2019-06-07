function [OP] = SM_RS_PP(KK,MM,LL,LdB,yS,xPU,yPU,xPB,yPB,rho,alpha,beta,eta,PdB,IdB,Rth,bit_frame)  
% hop by hop RS
OP          = zeros(1,length(PdB));
Ip          = 10.^(IdB/10);
PP          = 10.^(PdB/10);

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
PkmW1 = 0; PkmW2 = 0; PkmW3 = 0;

for aa = 1 : length(PdB)
    fprintf('Running %d per %d \n',aa,length(PdB));
    for bitnum   =  1 : bit_frame(aa)
        SNR_e2e = zeros(1,KK);
        for kk = 1 : KK
            % M PBs
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
            
            P_EHk  = kappa*PP(aa)*gmEk_sumk; % Power EH relay
            P_Ik   = Ip/gm_Ik;      % Interference Constraint
            Pk     = min(P_EHk,P_Ik);
            %%%%
            switch kk
                case 1
                   PkmW1 = PkmW1 + Pk;
                case 2
                   PkmW2 = PkmW2 + Pk;
                otherwise
                   PkmW3 = PkmW3 + Pk;
            end
            %%%%
            % Relay selection

           h_Di = sqrt(LD_Dk/2)*(randn(1,1)+1i*randn(1,1));                   
           SNR_e2e(kk)   = Pk*abs(h_Di)^2;         
        end
        CC_ii  = (1-alpha)/KK*log2(1 + min(SNR_e2e));
        if (CC_ii < Rth)            
            OP(aa) = OP(aa) + 1;
        end
    end
end
PkmW1 = PkmW1/bit_frame
PkmW2 = PkmW2/bit_frame
PkmW3 = PkmW3/bit_frame
OP  = OP./bit_frame;
semilogy(PdB,OP,'o'); grid on;hold on; 
end




