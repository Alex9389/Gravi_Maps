%Listing A.5:
% *****************************************************%
% %
% This function computes the Legendre function and its%
% first order derivatives using recursive formulae %
% %
% INPUT %
% phi= latitude of the desired point %
% Nmax = maximum desired degree %
% Mmax = maximum desired order %
% %
% OUTPUT %
% %
% pnm= Normalized associated Legendre function %
% dP= first order derivative of the Normalized %
% associated Legendre function %
% %
% REFERENCES %
% Hwang , C. and M.J. Lin , (1998) , Fast Integration %
% of low orbiter 's trajctory perturbed by the %
% Earth non - sphericity , Journal of Geodesy , %
% vol 72:578 -585 %
% %
% Borre , Kai (2004) , Geoid Undulations computed by%
% EGM96 , report , Aalborg University %
% %
% %
% By Mehdi Eshagh and Ramin Kiamehr 2006 %
% Division of Geodesy %
% Royal Institute of Technology %
% Stckholm , Sweden %
% Email : eshagh@kth .se %
% %
% *****************************************************%

function [pnm,dP]= Pnm(phi ,Nmax, Mmax)
  
  nrow = Nmax+1; np1 = Mmax+1;
  phii = phi*pi/180;
  
  x=sin(phii);
  y=cos(phii);
  
  pnm(1,1) =1.0;
  pnm(2,1)= sqrt(3.0)*x;
  pnm(2,2)= sqrt(3.0)*y;
  pnm(3,2)= sqrt(5.0)*pnm(2,2)*x;
  
  for i=3:np1
    n=i-1;
    pnm(i,i)= sqrt((n +0.5)/n)*pnm(i-1,i-1)*y;
  end
  k=np1-1;
 
  for i=3:k
    n=i-1;
    pnm(i+1,i)= sqrt(2.0* n+3)*pnm(i,i)*x;
  end
  nm1=np1-2;
  
  for j=1:nm1
    m=j-1;
    k=j+2;
    for i=k:np1
      n=i-1;
      c =(2.0*n+1.0)/(n-m)/(n+m);
      c1= c*(2.0*n-1.0) ;
      c2= c*(n+m-1)/(2*n -3)*(n-m-1);
      c1= sqrt(c1);
      c2= sqrt(c2);
      pnm(i,j)=c1*x*pnm(i -1,j)-c2*pnm(i -2,j);
    end
  end
  
  for n=1:Nmax-1
    for m=1:n
      dP(n,m)= sqrt(((n -1)-(m -1))*((n -1)+(m-1)+1)*(1+kron(m-1,m-1))*pnm(n,m+1)-(m-1)*tan(phi)*pnm(n,m));
    end
  end