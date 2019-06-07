function Lambda = fLambda(xx,zz,kappa,PP,LD_D,LD_E,Ip,OMG_I)
    Lambda = (xx/kappa./PP./LD_D/LD_E + Ip./kappa./PP./OMG_I/LD_E).^(zz/2).*besselk(zz,2.*sqrt(xx/kappa./PP./LD_D/LD_E + Ip./kappa./PP./OMG_I/LD_E));
end