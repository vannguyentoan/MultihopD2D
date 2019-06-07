% % % Fig 3 : OP vs PdB (ve theo so PRs)
clear all;
clc;
yS             = 0;
xPB            = 15;
yPB            = 5;
xPU            = 20;
yPU            = -5;
LdB            = 30; %% loss 30 dB
beta           = 2.7;
Rth            = 1;
eta            = 0.8;

rho            = 1; % imperfect CSI : 0.6 and 1
MM             = 3;   % number of power beacons
LL             = 3; 
KK             = 3; % K le for even, K chan for odd   K: so hops
Nk             = 2;
IdB            = 23;%-10:2:40; 23
PdB            = 22; %20;  0:5:45
% PdB            = 10:5:30; %20; 
alpha          = 0.2;
bit_frame      = 5*10^5*ones(1,length(PdB));

%%%############# Simulation ###############################
% % % 
% Best path RS
% OP_SM_BP = SM_BPRS_PP(KK,MM,Nk,LL,LdB,yS,xPU,yPU,xPB,yPB,rho,alpha,beta,eta,PdB,IdB,Rth,bit_frame);
% Dual-hop relay selection
% EVEN
% SM_DHRS_even_PP(KK,MM,Nk,LL,LdB,yS,xPU,yPU,xPB,yPB,rho,alpha,beta,eta,PdB,IdB,Rth,bit_frame);
% % % RS
SM_RS_PP(KK,MM,LL,LdB,yS,xPU,yPU,xPB,yPB,rho,alpha,beta,eta,PdB,IdB,Rth,bit_frame);
% save('OP_SM_BP1.mat','OP_SM_BP');

%%%############################################
% % % % %%% Theory
% % % % % %%% random Scheme
% Theory_RS_PP(KK,MM,LL,LdB,yS,xPU,yPU,xPB,yPB,rho,alpha,beta,eta,PdB,IdB,Rth);
% % % % % % % Theory_RS_AsymT(KK,MM,xPU,yPU,xPB,yPB,alpha,beta,eta,PdB,IdB,Rth);
% Theory_RS_AsymL(KK,LL,LdB,yS,xPU,yPU,xPB,yPB,rho,alpha,beta,eta,PdB,IdB,Rth);
% % % % % % % 
% % % % 
% % % % % % Dual-hop relay selection
% % % % 
% % % % % %%%%% EVEN
% Theory_DHRS_even_PP(KK,MM,Nk,LL,LdB,yS,xPU,yPU,xPB,yPB,rho,alpha,beta,eta,PdB,IdB,Rth);
% % % % % % % % Theory_PP_DHRS_even_AsymT(KK,MM,Nk,xPU,yPU,xPB,yPB,alpha,beta,eta,PdB,IdB,Rth);
% Theory_DHRS_even_AsymL_PP(KK,Nk,LL,LdB,yS,xPU,yPU,xPB,yPB,rho,alpha,beta,eta,PdB,IdB,Rth);
% % % % % % % 
% % % % % % % % % %%% Best path RS
% Theory_BPRS_PP(KK,MM,Nk,LL,LdB,yS,xPU,yPU,xPB,yPB,rho,alpha,beta,eta,PdB,IdB,Rth);
% % % % % % % Theory_PP_BPRS_AsymT(KK,MM,Nk,xPU,yPU,xPB,yPB,alpha,beta,eta,PdB,IdB,Rth);
% Theory_BPRS_AsymL_PP(KK,Nk,LL,LdB,yS,xPU,yPU,xPB,yPB,rho,alpha,beta,eta,PdB,IdB,Rth);




