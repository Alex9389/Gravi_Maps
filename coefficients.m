%Listing A.3:
% ***************************************************************%
% %
% This function computes the coefficients of the Legendre %
% functions %
% %
% INPUT %
% the maximum desired degree of geopotential model to be used %
% plus 2 -it is suggested to introduce higher values than the%
% maximum degree of the model. %
% %
% OUTPUT %
% %
% all of the coefficients of Legendre functions needed for %
% computing the gravity gradients %
% %
% REFERENCE %
% Petrovskaya , M.S. and A.N. Vershkov (2006) , Non/ Singular %
% expressions for the gravity gradients in the local %
% north - oriented and orbital referencse frames . Journal of %
% Geodesy , Vol 80, 117 -127% %
% %
% %
% by %
% Mehdi Eshagh and Ramin Kiamehr 2006 %
% Division of Geodesy %
% Royal Institute of Technology %
% Stockholm , Sweden %
% Email : eshagh@kth .se %
% %
% Modified by Marshall Rogers 2008 %
% %
% ***************************************************************%

function [a, b, c, d, g, h, beta, psi, mu, eta] = coefficients (N)
% Проверка полученного N из файла
if (N<10)  
	 error('ERROR \nchange file spherical harmonic model, wherein N greater 10')
end
% Смотри Appendix B. Legendre Function Coefficients used for EGM96 Gradient Calculations
% Здесь n-1 и m-1, потому что массивы в Matlab начинаются с первого элемента, а не нулевого
for n=1:N
  for m=1:n
  
    if ((m-1) == 0 || (m-1) == 1)
	  a(n,m) = 70;
	  b(n,m) = (n-1+(m -1)+1) * (n-1+(m-1)+2) / (2*(m-1)+1);
	  c(n,m) = sqrt(1+kron((m-1),0)) * sqrt((n-1)^2-((m-1)+1)^2) *...
               sqrt((n-1)-(m-1)) * sqrt(n-1+(m-1)+2) / 4;
    
    elseif (2 <= (m-1) <= (n-1))
      a(n,m) = sqrt(1+kron((m-1),2)) * sqrt((n-1)^2-((m-1)-1)^2) *...
               sqrt((n-1)+(m-1)) * sqrt(n-1-(m-1)+2) / 4;
      
      b(n,m) = ((n-1)^2 + (m-1)^2 + 3*(n-1) + 2) / 2;
      
      c(n,m) = sqrt((n-1)^2-((m-1)+1)^2) * sqrt((n-1)-(m-1)) *...
               sqrt((n-1)+(m -1)+2) / 4;
     
      d(n,m) = -(m-1)/4/(m -1) * sqrt((2*(n-1)+1) / (2*(n-1)-1)) * sqrt(1+kron(m-1,2)) *...
               sqrt((n-1)^2-((m-1)-1)^2) * sqrt(n-1+(m -1) ) * sqrt(n-1+(m -1)-2);
      
      g(n,m) = (m-1)/2 * sqrt((2*(n-1)+1) / (2*(n-1)-1)) * sqrt(n-1+(m-1)) * sqrt(n-1-(m-1));
     
      h(n,m) = 1/4 * sqrt((2*(n-1)+1)/(2*(n-1)-1)) * sqrt((n-1)^2-((m-1)+1)^2) *...
               sqrt(n-1-(m-1)) * sqrt(n-1-(m-1)-2);
    end
    
    if ((m-1) == 1)
      d(n,m) = 0;
      g(n,m) = 1/4 * sqrt((2*(n-1)+1)/(2*(n-1)-1)) * sqrt(n) * sqrt(n-2) * (n+1);
      h(n,m) = 1/4 * sqrt((2*(n-1)+1)/(2*(n-1)-1)) * sqrt(n-4) * sqrt(n-3) *...
	           sqrt(n-2) * sqrt(n+1);
    end
    
    if ((m-1) == 0)
      beta(n,m) = 0;
      psi(n,m)  = -(n+1) * sqrt((n-1)*n/2);
    elseif (1 <= (m-1) <= (n-1))
      beta(n,m) =  (n+1)/2 * sqrt(1+kron((m-1),1)) * sqrt(n-1+(m-1)) * sqrt(n-1-(m-1)+1);
      psi(n,m)  = -(n+1)/2 * sqrt(n-1-(m-1)) * sqrt(n-1+(m-1)+1);
    end
    
    if (0 <= (m-1) <= (n-1))
      mu(n,m) =  -(n+1)/2 * sqrt((2*(n-1)+1)/(2*(n-1)-1)) *...
                 sqrt(1+kron((m-1),1)) * sqrt(n-1+(m-1)) * sqrt(n-1+(m-1)-1);
      eta(n,m) = -(n+1)/2 * sqrt((2*(n-1)+1)/(2*(n-1)-1)) *...
	             sqrt(n-1-(m-1)) * sqrt(n-1-(m-1)-1);
    end
  end
end