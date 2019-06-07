function [OP] = SM_DHRS_even_PP(KK,MM,Nk,LL,LdB,yS,xPU,yPU,xPB,yPB,rho,alpha,beta,eta,PdB,IdB,Rth,bit_frame)  
% number of cluster is even (chan)(K-1 = 2,4,6,8...)
% --> number of hop odd (K = 3,5,7,9...)
% --> divide into(K-1)/2 dual hop and first hop
%
Ip          = 10.^(IdB/10);
PP          = 10.^(PdB/10);
OP          = zeros(1,length(PdB));
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

for aa = 1 : length(PdB)
    fprintf('Running %d per %d \n',aa,length(PdB));
    for bitnum   =  1 : bit_frame(aa)
        SNR_e2e = zeros(1,(KK-1)/2+1);
        for bb = 1 : (KK-1)/2+1 % include first hop and (K-1)/2 dual-hop
            kk = 2*bb;
            if bb == 1
                % the first hop
                gmEk_sum1 = 0;
                for mm = 1 : MM
                    h_Ei = sqrt(LD_Ek(bb)/2)*(randn(1,1)+1i*randn(1,1));
                    gmEk_sum1 = gmEk_sum1 + abs(h_Ei)^2;
                end
                %  L primary users (imperfect channel gain)
                h_Ii = zeros(1,LL);
                for ii = 1 : LL      
                  h_Ii(1,ii)  = sqrt(LD_Ik(bb)/2)*(randn(1,1)+1i*randn(1,1));
                end
                Eps     = sqrt(LD_Ik(bb)/2)*(randn(1,1)+1i*randn(1,1));
                h_I     = rho*max(h_Ii) + sqrt(1-rho^2)*Eps;   % imperfect channel 
                gm_I1   = abs(h_I)^2;

                P_EH1    = kappa*PP(aa)*gmEk_sum1; % Power EH relay
                P_I1     = Ip/gm_I1;      % Interference Constraint
                Pk1      = min(P_EH1,P_I1);

                gm_h1 = zeros(1,Nk);
                for nn = 1 : Nk
                   h_Di = sqrt(LD_Dk/2)*(randn(1,1)+1i*randn(1,1));                   
                   gm_h1(1,nn)   = Pk1*abs(h_Di)^2;
                end
                SNR_max = max(gm_h1); % Relay selection 
            else
                % the next dual hop
                SNR_candidate = zeros(1,Nk);
                for jj = 1 : Nk
                    % hop 2bb-2

                    gmEk_sumA = 0;
                    for mm = 1 : MM
                        h_Ei = sqrt(LD_Ek(kk-2)/2)*(randn(1,1)+1i*randn(1,1));
                        gmEk_sumA = gmEk_sumA + abs(h_Ei)^2;
                    end

                    %  L primary users (imperfect channel gain)
                    h_Ii = zeros(1,LL);
                    for ii = 1 : LL      
                      h_Ii(1,ii)  = sqrt(LD_Ik(kk-2)/2)*(randn(1,1)+1i*randn(1,1));
                    end
                    Eps    = sqrt(LD_Ik(kk-2)/2)*(randn(1,1)+1i*randn(1,1));
                    h_I    = rho*max(h_Ii) + sqrt(1-rho^2)*Eps;   % imperfect channel 
                    gm_IA  = abs(h_I)^2;
                    P_EHA  = kappa*PP(aa)*gmEk_sumA; % Power EH relay
                    P_IA   = Ip/gm_IA;      % Interference Constraint
                    PkA    = min(P_EHA,P_IA);
                    h_Di   = sqrt(LD_Dk/2)*(randn(1,1)+1i*randn(1,1));                   
                    gm_hA  = PkA*abs(h_Di)^2;
                    
                    % hop 2bb-1
                    gmEk_sumB = 0;
                    for mm = 1 : MM
                        h_Ei = sqrt(LD_Ek(kk-1)/2)*(randn(1,1)+1i*randn(1,1));
                        gmEk_sumB = gmEk_sumB + abs(h_Ei)^2;
                    end
                    %  L primary users (imperfect channel gain)
                    h_Ii = zeros(1,LL);
                    for ii = 1 : LL      
                      h_Ii(1,ii)   = sqrt(LD_Ik(kk-1)/2)*(randn(1,1)+1i*randn(1,1));
                    end
                    Eps    = sqrt(LD_Ik(kk-1)/2)*(randn(1,1)+1i*randn(1,1));
                    h_I    = rho*max(h_Ii) + sqrt(1-rho^2)*Eps;   % imperfect channel 
                    gm_IB  = abs(h_I)^2;
                    P_EHB  = kappa*PP(aa)*gmEk_sumB; % Power EH relay
                    P_IB   = Ip/gm_IB;      % Interference Constraint
                    PkB    = min(P_EHB,P_IB);
                    
                    gm_hk = zeros(1,Nk);
                    for nn = 1 : Nk
                       h_Di = sqrt(LD_Dk/2)*(randn(1,1)+1i*randn(1,1));                   
                       gm_hk(1,nn)   = PkB*abs(h_Di)^2;
                    end            
                    gm_hB  = max(gm_hk);
                    
                    SNR_candidate(1,jj) = min(gm_hA, gm_hB);
                end
                SNR_max = max(SNR_candidate); %% max-min criteria
            end
            SNR_e2e(bb)= SNR_max;
        end
        CC_ii  = (1-alpha)/KK*log2(1 + min(SNR_e2e));
        if (CC_ii < Rth)            
            OP(aa) = OP(aa) + 1;
        end
    end
end
OP  = OP./bit_frame;
semilogy(PdB,OP,'o'); grid on;hold on; 
end




