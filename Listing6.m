%Listing A.6:
% ******************************* %
% Gradients due to KC -10 Tanker  %
% ( used to generate white noise )%
% ******************************* %
close all
clear all
clc
% z positive "downward"
G =6.67E-11;
p =132; % density contrast for tanker
Eotvos =1E-9; % Eotvos conversion
a =[-3.05+25 3.05+25]; % Point 1, Permuted for 7 remaining points
b =[-5-24.4 -5+24.4]; % Point 1, Permuted for 7 remaining points
c =[-5-3.05 -5+3.05]; % Point 1, Permuted for 7 remaining points
alt =0;
Ng  =0;
Eg  =0;
sumTxx =0;

for ctr1 =1:2
  for ctr2 =1:2
    for ctr3 =1:2
      sumTxx = sumTxx +(( -1) ^( ctr1 )*( -1) ^( ctr2 )*( -1) ^( ctr3 )).*atan2 (((Ng -b( ctr2 )).*( alt -c( ctr3 ))) ,((Eg -a( ctr1 ))...
                                                               .*(( Eg -a( ctr1 )) .^2+( Ng -b( ctr2 )) .^2+( alt -c( ctr3 )).^2) .^(.5) ));
    end
  end
end
Txx =G*p*sumTxx/Eotvos;
sumTxy =0;

for ctr1 =1:2
  for ctr2 =1:2
    for ctr3 =1:2
      sumTxy = sumTxy +(( -1) ^( ctr1 + ctr2 + ctr3 )).* log ((( alt -c(...
      ctr3 ))+(( Eg -a( ctr1 )) .^2+( Ng -b( ctr2 )) .^2+( alt -c( ctr3)) .^2) .^(.5) ));
    end
  end
end
Txy =-G*p*sumTxy;
sumTxz =0;

for ctr1 =1:2
  for ctr2 =1:2
    for ctr3 =1:2
      sumTxz = sumTxz +(( -1) ^( ctr1 )*( -1) ^( ctr2 )*( -1) ^( ctr3 )).*...
      log ((( Ng -b( ctr2 ))+((Eg -a( ctr1 )) .^2+( Ng -b( ctr2 )).^2+( c( ctr3 )-alt) .^2) .^(.5) ));
    end
  end
end
Txz =-G*p*sumTxz/Eotvos;
sumTyy =0;

for ctr1 =1:2
  for ctr2 =1:2
    for ctr3 =1:2
      sumTyy = sumTyy +(( -1) ^( ctr1 )*( -1) ^( ctr2 )*( -1) ^( ctr3 )).*atan2 (((Eg -a( ctr1 )).*( alt -c( ctr3 ))) ,((Ng -b( ctr2 ))...
                                                               .*(( Eg -a( ctr1 )) .^2+( Ng -b( ctr2 )) .^2+( c( ctr3 )-alt ).^2) .^(.5) ));
    end
  end
end
Tyy =G*p*sumTyy/Eotvos;
sumTyz =0;

for ctr1 =1:2
  for ctr2 =1:2
    for ctr3 =1:2
      sumTyz = sumTyz +(( -1) ^( ctr1 )*( -1) ^( ctr2 )*( -1) ^( ctr3 )).*...
      log ((( Eg -a( ctr1 ))+((Eg -a( ctr1 )) .^2+( Ng -b( ctr2 )).^2+( c( ctr3 )-alt) .^2) .^(.5) ));
    end
  end
end
Tyz =-G*p*sumTyz/Eotvos;
Tzz =-(Tyy+Txx);