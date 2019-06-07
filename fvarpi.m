%%% varpi
function varpi = fvarpi(Ip,kappa,PP,OMG_Ik,LD_Ek)
    varpi = Ip./kappa./PP./OMG_Ik/LD_Ek;
end