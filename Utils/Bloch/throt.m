function Rth=throt(phi,theta)

Rz = zrot(-theta);
Rx = xrot(phi);
Rth = inv(Rz)*Rx*Rz;


