function psi = fpsi(xx,yy)
    psi = nchoosek(xx,yy)*(-1)^(yy-1);
end