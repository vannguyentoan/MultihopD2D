%%% varpi
function vartheta = fvartheta(xx,yy)
    vartheta = nchoosek(xx-1,yy)*2*xx*(-1)^yy/(yy+1);
end