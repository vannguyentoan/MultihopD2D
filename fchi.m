function chi = fchi(xx,yy,zz,Ip)
    chi = xx.*Ip./(yy.*zz + xx.*Ip);
end