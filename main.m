% **************************************************%
% Gravity Gradient Signal and " Truth " Calculation %
% **************************************************%
close all
clear all
clc

G =6.67E-11; % Universal Gravity Const
p =2670; % Average terrain density Средняя плотность местности???
Eotvos=1E-9; %use to convert units to Eotvos

data_of_map = hmap('..\..\MAPS\new_map.mtw');
x0 = data_of_map(4);
y0 = data_of_map(5);
x1 = data_of_map(6);
y1 = data_of_map(7);
b0 = data_of_map(8)/pi*180;
l0 = data_of_map(9)/pi*180;
b1 = data_of_map(10)/pi*180;
l1 = data_of_map(11)/pi*180;
fprintf ('Size MAP:\n');
fprintf ('   lat: [%g %g]\n', b0, b1);
fprintf ('  long: [%g %g]\n', l0, l1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ПОЛЬЗОВАТЕЛЬСКАЯ НАСТРОЙКА %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Задание размера области, минимум 600*600 элементов для 3 секундной карты
%(иначе нельзя проинтерполировать 30 минутную карту (EGM96) до 3 секундной)
size_map_lat  = 2000; %размер карты по широте (в элементах) 
size_map_long = 2000; %размrер карты по долготе (в элементах)
step_map_sec  = 3;    %шаг карты (в секундах)
alt = 25000;          %высота,м
cut_map = 100; % ЗАПЛАТКА: почему-то на краях карты наблюдаются ошибки вычисления
% градиента, потому их обрезаем. Установить 0, если алгоритм исправлен!

M= input ('Rough or Smooth Terrain [R/S]? ','s');
if ((M== 'R')||(M== 'r'))
  % Начальные координаты для пострения участка с грубым рельефом
  start_lat_grad = 41.851401;
  start_long_grad  = 41.272401;
  M = 'R';
elseif ((M== 'S')||(M== 's'))
  % Начальные координаты для пострения участка с гладким рельефом
  start_lat_grad = 40.351432;  
  start_long_grad  = 40.272432;
  M = 'S';
else
  error('Wrong Answer !');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%для записи в файл данных
size_map_lat_r  = size_map_lat;
size_map_long_r = size_map_long;
%для алгоритма
size_map_lat  = size_map_lat +2*cut_map;
size_map_long = size_map_long+2*cut_map;
% Расчет границы выбранного участка карты
%для записи в файл данных
lat_r  = [start_lat_grad  (start_lat_grad +size_map_lat_r *step_map_sec/3600)];
long_r = [start_long_grad (start_long_grad+size_map_long_r*step_map_sec/3600)];
%для алгоритма
lat  = [(start_lat_grad -cut_map*step_map_sec/3600)  (start_lat_grad +(size_map_lat -cut_map)*step_map_sec/3600)];
long = [(start_long_grad-cut_map*step_map_sec/3600)  (start_long_grad+(size_map_long-cut_map)*step_map_sec/3600)];
% Проверка на выход выбранного участка рельефа за границы загруженной карты
if ((lat(1)<b0) || (lat(2)>b1) || (long(1)<l0) || (long(2)>l1))
  error ('Wrong coords!');
end

%% Загрузка выбранного участка карты
Z = load_map(start_lat_grad -cut_map*step_map_sec/3600, start_long_grad-cut_map*step_map_sec/3600, size_map_lat, size_map_long, step_map_sec);

figure
mesh(Z); 
axis([cut_map size_map_lat-cut_map cut_map size_map_long-cut_map]);
title('Карта рельефа местности');

fprintf ('Size gradient MAP:\n');
fprintf ('   lat: [%g %g]\n', lat_r);
fprintf ('  long: [%g %g]\n', long_r);

alt_km= floor(alt/1000);
alt = mean2(Z)+alt; % Height above average terrain calc mean2-среднее
alt_EGM96 =alt+ mean2(Z);
%перевод высоты из численного формата в символьный (необходимо для файлов данных)
if (alt_km<10)
  alt_km_char = ['0' num2str(alt_km)];
else
  alt_km_char = num2str(alt_km);
end

% Запись в файл информации о используемом участке карты рельефа
file_data_of_map=[M 'DATA' alt_km_char];
fid = fopen(['..\..\MAPS\GradientMaps\' file_data_of_map],'w'); % Opening a file for the EGM96 Gradients
fprintf (fid ,'%d %d %d %10.6f %10.6f %10.6f %10.6f \n',size_map_lat_r, size_map_long_r, step_map_sec, lat_r, long_r);
fclose (fid);

%% Computer Earth 's Radius at central grid point
a = 6378137.0; % earth semimajor axis in meters
f = 1/298.257223563; % reciprocal flattening - сжатие сферы вдоль диаметра с образованием эллипсоида

e2 = 2*f-f^2; % eccentricity squared - эксцентриситет в квадрате
lat_middle =(min(lat)+max(lat))/2;
lat_middle2 = atand ((1 - e2)* tand ( lat_middle )); % convert geodetic lat to geocentric lat
long_middle =( min( long )+ max( long )) /2;
N=a/ sqrt (1- e2*sin( lat_middle2 )^2);

X_ECEF =N* cosd ( lat_middle2 )*cos ( long_middle );
Y_ECEF =N*cos ( lat_middle2 )*sin ( long_middle );
Z_ECEF =(N*(1 - e2))*sin ( lat_middle2 );

R1= sqrt ( X_ECEF ^2+ Y_ECEF ^2+ Z_ECEF ^2) ; % Earth 's radius at grid midpoint
dlong =max(abs(long )) - min(abs(long)); % distance longitude
dlat =max(lat) - min(lat); % distance latitude

x_dist =R1 *( pi /180) * dlong * cosd ( lat_middle );
y_dist =R1 *( pi /180) * dlat ;
y_int =R1 *( pi /180) *(step_map_sec/3600) ; %3 arc sec spacing ( DTED Level 1)3/3600, Level 0 - 30/3600
x_int = cosd ( lat_middle )* y_int;

Y1 = 0: y_int : y_dist;
X1 = 0: x_int : x_dist;
% Проверка! почему-то периодически выходит за границы 
if (length(Y1)>size_map_lat) || (length(X1)>size_map_long)  
  Y1 = y_int: y_int : y_dist;
  X1 = x_int: x_int : x_dist;
end 
[X2 ,Y2 ]= meshgrid (X1 ,Y1);

%% Gradients by Parker 's Methods ( Terrain )
del_x1 = x_int ; % x interval
del_x2 = y_int ; % y interval
m1 = length (X1);
m2 = length (Y1);
p1 = [0:1: m1 -1];
p2 = [0:1: m2 -1];
for ctr =1: m1
  if p1( ctr) <= (m1 /2)-1
    f1_p(ctr)=p1( ctr)/( del_x1 *m1); % spatial frequency x
    f2_p(ctr)=p2( ctr)/( del_x2 *m2); % spatial frequency y
  else
    f1_p(ctr) =(p1( ctr)-m1)/(del_x1 *m1); % spatial frequency x
    f2_p(ctr)=(p2( ctr)-m2)/(del_x2 *m2); % spatial frequency y
  end
end

[f1m_p, f2m_p]= meshgrid (f1_p, f2_p); % gridded spatial freq data
f_p = sqrt(f1m_p^2+ f2m_p^2);
sig =0;
for ctr =1:20
%  sig = sig+((1./ factorial(ctr)).*(2*pi.*f_p).^(ctr-2).*fft2((Z.^ctr))); % perform FFT (see Jekeli & Zhu)
  sig1 = (1./ factorial(ctr)).*(2*pi.*f_p);
  sig2 = sig1.^(ctr-2);
  sig3 = fft2((Z.^ctr));
  sig = sig + sig2.*sig3;
end

sig (1) =1E15; % prevent infinite value
mu_xx = -((2* pi)^2) .* f1m_p .^2;
mu_xy = -((2* pi)^2) .* f1m_p .* f2m_p ;
mu_xz =i *((2* pi)^2) .* f1m_p .* f_p;
mu_yy = -((2* pi)^2) .* f2m_p .^2;
mu_yz =i *((2* pi)^2) .* f2m_p .* f_p;
mu_zz =((2* pi)^2) .* f_p .^2;
Txx_parker =(2* pi*p*G.* ifft2 ( mu_xx .* exp ( -2* pi*alt .* f_p).* sig));
Txy_parker =(2* pi*p*G.* ifft2 ( mu_xy .* exp ( -2* pi*alt .* f_p).* sig));
Txz_parker =(2* pi*p*G.* ifft2 ( mu_xz .* exp ( -2* pi*alt .* f_p).* sig));
Tyy_parker =(2* pi*p*G.* ifft2 ( mu_yy .* exp ( -2* pi*alt .* f_p).* sig));
Tyz_parker =(2* pi*p*G.* ifft2 ( mu_yz .* exp ( -2* pi*alt .* f_p).* sig));
Tzz_parker =(2* pi*p*G.* ifft2 ( mu_zz .* exp ( -2* pi*alt .* f_p).* sig));
Txx_parker = real ( Txx_parker )./ Eotvos ;
Txy_parker = real ( Txy_parker )./ Eotvos ;
Txz_parker = real ( Txz_parker )./ Eotvos ;
Tyy_parker = real ( Tyy_parker )./ Eotvos ;
Tyz_parker = real ( Tyz_parker )./ Eotvos ;
Tzz_parker = real ( Tzz_parker )./ Eotvos ;

%% Gradients from EGM96 ( Long wavelength subterranean effects )
% Taken from Kiamehr & Eshagh and Modified
% Could also use geopot97 .v0 .4e.f code
phi_south =min(lat);
phi_north =max(lat);
lambda_west = min(long);
lambda_east = max(long);
phi_step =30/60; % EGM96 provides a 30 arcmin resolution 
lambda_step =30/60; % EGM96 provides a 30 arcmin resolution 
X_EGM96 = lambda_west+step_map_sec/3600 :step_map_sec/3600: lambda_east ; % setup array for griddata function у нас 30
Y_EGM96 = phi_south+step_map_sec/3600 :step_map_sec/3600: phi_north ; % setup array for griddata function у нас 30
[X2_EGM96 , Y2_EGM96]= meshgrid (X_EGM96, Y_EGM96);

[Nmax ,Ae ,GM ,C,S,dC ,dS] = Modelread('..\..\MAPS\geocoeff.txt'); % Read spherical harmonic model  %
[a,b,c,d,g,h,beta ,psi ,mu ,eta] = coefficients(Nmax+3); % calculateLegendre coeffs.            
CN= Normal(GM,Ae,Nmax);
for i = 3:1:11
  C(i,1) =C(i,1) - CN(i); % Generation of the Potetial Anomaly                              
end

filename ='EGM96Gradients';
fid = fopen(filename,'w'); % Opening a file for the EGM96 Gradients

for phi = phi_south : phi_step : phi_north
  phigeodetic = phi;
  phi =phi*pi/180;

  % Compute the Geocentric latitude via geodetic latitude
  e2 =.00669437999013; %1st eccentricity squared
  phi = atan((1-e2)*tan(phi));

  % Compute the Associated Legendre functions 
  [pnm ,dP ]= Pnm(phi *180/ pi , Nmax +3, Nmax +3);                                          
  
  for lambda = lambda_west : lambda_step : lambda_east
    lambda = lambda *pi /180;
    sum    =0; % Initialize summations
    sum1   =0;
    sum2   =0;
    sum3   =0;
    sum4   =0;
    sum5   =0;
    sumN   =0;
    sumdg  =0;
    sumeta =0;
    sumpsi =0;

    % Computation of geocenric distance
    N=Ae/ sqrt (1- e2*sin(phi)^2);                                                              
    X_ECEF =(N+ alt_EGM96 )*cos (phi )*cos ( lambda );
    Y_ECEF =(N+ alt_EGM96 )*cos (phi )*sin ( lambda );
    Z_ECEF =(( N+ alt_EGM96 )*(1 - e2))*sin (phi );
    r= sqrt ( X_ECEF ^2+ Y_ECEF ^2+ Z_ECEF ^2);

    for n=3: Nmax +1
      for m=1:n
        CS =(C(n,m)*cos ((m -1) * lambda )+S(n,m)*sin ((m -1) *lambda ));
        AA =( Ae/r)^(n+2) ;
        AA1 =( Ae/r)^n;
        CS1 =(-S(n,m)*cos ((m -1) * lambda )+C(n,m)*sin ((m -1) * lambda ));
        PNM =pnm (n,m);
        if (abs (m) -2) <= 0
          PP =( -1) ^( abs(m -4) -1)*pnm (n,abs ((m) -4));
          PP1 =( -1) ^( abs(m -4) -1)*pnm(n -1, abs ((m) -4));
        else
          PP= pnm(n, abs(m) -2);
          PP1 =pnm (n -1, abs(m) -2);
        end
    
        if  (abs(m) -1) <= 0
          QQ =( -1) ^( abs(m -3) -1)*pnm (n,abs (abs (m) -3));
          QQ1 =( -1) ^( abs(m -3) -1)*pnm(n -1, abs(abs(m) -3));
        else
          QQ= pnm(n, abs(m) -1);
          QQ1 =pnm (n -1, abs(m) -1);
        end
        % Computing the Txx summation
        sum1 = sum1 +AA*CS *(a(n,abs (m))*PP+b(n,abs (m))*pnm(n,abs(m)) + c(n,abs(m))*pnm(n,abs(m)+2));
        % Computing the Txy summation
        sum3 = sum3 +AA*CS1 *(d(n,m)*PP1 +g(n,m)*pnm (n -1 ,(m))+h(n,m)*pnm(n -1 ,(m)+2));
        % Computing the Txz summation
        sum4 = sum4 +AA*CS *( beta (n,m)*QQ+psi (n,m)*pnm (n ,(m)+1));
        % Computing the Tyz summation
        sum5 = sum5 +AA*CS1 *( mu(n,m)*QQ1 +eta (n,m)*pnm (n-1 ,(m) +1));
        % Computing the Tzz summation
        sum2 = sum2 +(n*(n+1) )*AA*CS*PNM ;
      end % of m
    end % of n

    %The gravity gradient tensor components
    Txx =-GM/Ae ^3* sum1 / Eotvos ;
    Tzz = GM/Ae ^3* sum2 / Eotvos ;
    Txy = (-GM/Ae ^3* sum3 / Eotvos ) /10;
    Tyz = -GM/Ae ^3* sum4 / Eotvos ;
    Txz = GM/Ae ^3* sum5 / Eotvos ;
    fprintf (fid ,'%g %g %e %e %e %e %e %e \n',phigeodetic ,lambda*180/pi ,Txx ,-(Txx+Tzz), Tzz, Txy,Txz,Tyz );
  end % of lambda
end % of phi
fclose (fid);
U= load (filename);

% Интерполяция составляющих тензора гравитации, полученных из гармоник (EGM96)
% т.к. точность расчитанных данных порядка 30мин, на выходе необходимо получить 3 сек.
Txx_EGM96 = griddata (U(: ,2) ,U(: ,1) ,U(: ,3) ,X2_EGM96 , Y2_EGM96 ,'v4');
Tyy_EGM96 = griddata (U(: ,2) ,U(: ,1) ,U(: ,4) ,X2_EGM96 , Y2_EGM96 ,'v4');
Tzz_EGM96 = griddata (U(: ,2) ,U(: ,1) ,U(: ,5) ,X2_EGM96 , Y2_EGM96 ,'v4');
Txy_EGM96 = griddata (U(: ,2) ,U(: ,1) ,U(: ,6) ,X2_EGM96 , Y2_EGM96 ,'v4');
Txz_EGM96 = griddata (U(: ,2) ,U(: ,1) ,U(: ,7) ,X2_EGM96 , Y2_EGM96 ,'v4');
Tyz_EGM96 = griddata (U(: ,2) ,U(: ,1) ,U(: ,8) ,X2_EGM96 , Y2_EGM96 ,'v4');

% %%Получение градиентов
Txx = Txx_parker + Txx_EGM96;
Txy = Txy_parker + Txy_EGM96;
Txz = Txz_parker + Txz_EGM96;
Tyy = Tyy_parker + Tyy_EGM96;
Tyz = Tyz_parker + Tyz_EGM96;
Tzz = Tzz_parker + Tzz_EGM96;

%%Получение градиентов
% Txx = Txx_EGM96;
% Txy = Txy_EGM96;
% Txz = Txz_EGM96;
% Tyy = Tyy_EGM96;
% Tyz = Tyz_EGM96;
% Tzz = Tzz_EGM96;

% ЗАПЛАТКА: почему-то на краях карты наблюдаются ошибки вычисления
% градиента, потому их обрезаем
Txx = abs(Txx(cut_map+1:size_map_lat-cut_map,cut_map+1:size_map_long-cut_map));
Txy = abs(Txy(cut_map+1:size_map_lat-cut_map,cut_map+1:size_map_long-cut_map));
Txz = abs(Txz(cut_map+1:size_map_lat-cut_map,cut_map+1:size_map_long-cut_map));
Tyy = abs(Tyy(cut_map+1:size_map_lat-cut_map,cut_map+1:size_map_long-cut_map));
Tyz = abs(Tyz(cut_map+1:size_map_lat-cut_map,cut_map+1:size_map_long-cut_map));
Tzz = abs(Tzz(cut_map+1:size_map_lat-cut_map,cut_map+1:size_map_long-cut_map));

%%Запись составляющих тензора гравитации
map_file_Txx=[M 'TxxFL' alt_km_char];
dlmwrite(['..\..\MAPS\GradientMaps\' map_file_Txx], Txx, 'delimiter', '\t', 'precision', 8);

map_file_Txy=[M 'TxyFL' alt_km_char];
dlmwrite(['..\..\MAPS\GradientMaps\' map_file_Txy], Txy, 'delimiter', '\t', 'precision', 8);

map_file_Txz=[M 'TxzFL' alt_km_char];
dlmwrite(['..\..\MAPS\GradientMaps\' map_file_Txz], Txz, 'delimiter', '\t', 'precision', 8);

map_file_Tyy=[M 'TyyFL' alt_km_char];
dlmwrite(['..\..\MAPS\GradientMaps\' map_file_Tyy], Tyy, 'delimiter', '\t', 'precision', 8);

map_file_Tyz=[M 'TyzFL' alt_km_char];
dlmwrite(['..\..\MAPS\GradientMaps\' map_file_Tyz], Tyz, 'delimiter', '\t', 'precision', 8);

map_file_Tzz=[M 'TzzFL' alt_km_char];
dlmwrite(['..\..\MAPS\GradientMaps\' map_file_Tzz], Tzz, 'delimiter', '\t', 'precision', 8);

%%Построение 6 составляющих гравитационного тензора на одной вкладке
% figure
% subplot (3,2,1);
% Z = Z(cut_map+1:size_map_lat-cut_map,cut_map+1:size_map_long-cut_map);
% mesh(Z); 
% % axis([cut_map size_map_lat-cut_map cut_map size_map_long-cut_map]);
% axis([0 size_map_lat-2*cut_map 0 size_map_long-2*cut_map]);
% xlabel('Долгота');
% ylabel('Широта');
% title('Карта рельефа местности');
% set(axis, 'XTick', [0 size_map_lat-2*cut_map]);
% set(axis, 'YTick', [0 size_map_long-2*cut_map]);
% set(axis, 'XTickLabel', { start_lat_grad;   (start_lat_grad+(size_map_lat -2*cut_map)*3/3600) });
% set(axis, 'YTickLabel', { start_long_grad; (start_long_grad+(size_map_long-2*cut_map)*3/3600) });

% subplot (3,2,2);
% mesh(Txx);
% axis([0 size_map_lat-2*cut_map 0 size_map_long-2*cut_map]);
% %axis([ start_lat_grad (start_lat_grad+(size_map_lat-2*cut_map)*3/3600) start_long_grad (start_long_grad+(size_map_long-2*cut_map)*3/3600) ]);
% xlabel('Долгота');
% ylabel('Широта');
% title('Составляющая Txx');
% 
% subplot (3,2,3);
% mesh(Txy);
% axis([0 size_map_lat-2*cut_map 0 size_map_long-2*cut_map]);
% xlabel('Долгота');
% ylabel('Широта');
% title('Составляющая Txy');
% 
% subplot (3,2,4);
% mesh(Txz);
% axis([0 size_map_lat-2*cut_map 0 size_map_long-2*cut_map]);
% xlabel('Долгота');
% ylabel('Широта');
% title('Составляющая Txz');
% 
% subplot (3,2,5);
% mesh(Tyz);
% axis([0 size_map_lat-2*cut_map 0 size_map_long-2*cut_map]);
% xlabel('Долгота');
% ylabel('Широта');
% title('Составляющая Tyz');
% 
% subplot (3,2,6);
% mesh(Tzz);
% axis([0 size_map_lat-2*cut_map 0 size_map_long-2*cut_map]);
% xlabel('Долгота');
% ylabel('Широта');
% title('Составляющая Tzz');

fprintf ('Maps Successfully Created!\n');