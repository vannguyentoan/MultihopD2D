% OMG_I = (1+ii-ii*rho^2)*LD_Ik(kk)/(ii+1);
function Xi = fXi(xx,yy,zz,kappa,PP,LD)
    Xi = (xx/kappa./PP./yy/LD).^(zz/2).*besselk(zz,2*sqrt(xx/kappa./PP./yy/LD));
end