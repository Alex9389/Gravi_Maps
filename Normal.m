%Listing A.4:
% This function computes the normal field potential coefficients
function [J]= Normal (GM,AX,Nmax)
  GMS  = 0.3986005e15;
  AXS  = 6378137.0;
  JJ   = 0.108262982131e-2;
  FINV = 298.257;
  FLTN = 1.0/ FINV;
  E2   = FLTN*(2.0-FLTN);
  J(1) = GMS/GM;
  J(3) = -0.484169650276e-3;
  J(5) =  0.790314704521e-6;
  J(7) = -0.168729437964e-8;
  J(9) =  0.346071647263e-11;
  J(11)= -0.265086254269e-14;
  J(14) = 0;