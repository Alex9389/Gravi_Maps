%Listing A.7:
% ****************************%
% Terrain Avoidance Scenario  %
% ****************************%
close all
clear all
clc
vi =[50 0 0];
speed = norm (vi);
G =6.67E-11;
p =2670; % density contrast ground

% p =1000; % density contrast water tower
spacing =1;
N = [0: spacing :1550]; % setup grid
E = [0: spacing :1550]; % setup grid
Eotvos = 1E-9; % Eotvos conversion
update_rate   = 1;
filter_cutoff =.2;

%% 25m Cubic Object - dimensions permuted for each obstacle size, water tower gradients calculated similarly
length1 =25; % object base ( assumes square )
c=[0 25]; % obstacle height
t1 =0: spacing / speed :( round ( length (N)/2) -( length1 / spacing )/2 -1)/speed; % Truth signal time array

t2 =0:1/ update_rate :( round ( length (N)/2) -length1 /2 -1)/ speed; % time array for sensor with 1Hz sampling rate
a =[1000 1000+ length1 ]; % putting obstacle on grid
b=[round ( length (N)/2) -length1 /2 round ( length (N)/2)+ length1 /2]; % putting obstacle on grid

% obstacle height
alt =15; % 10m below obstacle top
ipos =[0 (length(N)/2) alt];
[Eg ,Ng]= meshgrid(E,N);

%Txx
sumTxx =0;
for ctr1 =1:2
  for ctr2 =1:2
    for ctr3 =1:2
      sumTxx = sumTxx +(( -1) ^( ctr1 )*( -1) ^( ctr2 )*( -1) ^( ctr3 )).*atan2 (((Ng -b( ctr2 )).*( alt -c( ctr3 ))) ,((Eg -a( ctr1 ))...
                                                              .*(( Eg -a( ctr1 )) .^2+( Ng -b( ctr2 )) .^2+( alt -c( ctr3 )).^2) .^(.5) ));
    end
  end
end
Txx =G*p* sumTxx / Eotvos ;

% Txy
sumTxy =0;
for ctr1 =1:2
  for ctr2 =1:2
    for ctr3 =1:2
      sumTxy = sumTxy +(( -1) ^( ctr1 + ctr2 + ctr3 )).* log ((( alt -c(...
      ctr3 ))+(( Eg -a( ctr1 )) .^2+( Ng -b( ctr2 )) .^2+( alt -c(ctr3)) .^2) .^(.5) ));
    end
  end
end
Txy =-G*p* sumTxy ;
Txy (727 ,1051) = -0.000001; % Prevent infinite values
Txy (727 ,1001) =  0.000001; % Prevent infinite values
Txy (777 ,1051) =  0.000001; % Prevent infinite values
Txy (777 ,1001) = -0.000001; % Prevent infinite values
Txy =Txy/Eotvos;

% Txz
sumTxz =0;
for ctr1 =1:2
  for ctr2 =1:2
    for ctr3 =1:2
      sumTxz = sumTxz +(( -1) ^( ctr1 )*( -1) ^( ctr2 )*( -1) ^( ctr3 )).*...
      log ((( Ng -b( ctr2 ))+((Eg -a( ctr1 )) .^2+( Ng -b( ctr2 )).^2+( c( ctr3 )-alt) .^2) .^(.5) ));
    end
  end
end
Txz =G*p*sumTxz/Eotvos;

% Tyy
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

% Tyz
sumTyz =0;
for ctr1 =1:2
  for ctr2 =1:2
    for ctr3 =1:2
      sumTyz = sumTyz +(( -1) ^( ctr1 )*( -1) ^( ctr2 )*( -1) ^( ctr3 )).*...
      log ((( Eg -a( ctr1 ))+((Eg -a( ctr1 )) .^2+( Ng -b( ctr2 )).^2+( c( ctr3 )-alt) .^2) .^(.5) ));
    end
  end
end
Tyz =G*p*sumTyz/Eotvos;

% Tzz
Tzz =-(Tyy+Txx);
% Simulation
sim (' terrain_avoidance_sim ', 0:.01:30) ; %run sim for 30s