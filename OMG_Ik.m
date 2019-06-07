%%% OMG_I
function OMG = OMG_Ik(ii,rho,LD_Ik)
    OMG = (1+ii-ii*rho^2)*LD_Ik/(ii+1);
end