%% generating a circle
format long
%% generating random cylinder cuts
n=50;
p=10;
dat=1000;
%% generating random cylinder cuts
x_values =[];
y_values=[];
p_values=[];
for i = 1 : 1 : dat

x=0;% centered at the origin
y=0;
th= 0:pi/n:(2*pi-0.00000001); % pi 'ın paydasını değiştirince MoMa yaklaşıyoruz
a=rand(1)/15; % to randomize the shape of our structure
b=rand(1)/15;
be=rand(1)/15;
ber=rand(1)/15;
berk=rand(1)/15;
r=1-(a*(cos(8*th)).^2)+(b*sin(4*th)+((be*cos(4*th)).^3+(ber*(sin(3*th)).^5+(berk*cos(8*th)))));%the radius function that has a function in it to generate a random shape
circle(x,y,r,n); 
x_value =(r .* cos(th) + x);
y_value =(r .* sin(th) + y);
p_value = sqrt((x_value.^2)+(y_value.^2));
x_values = [x_values; x_value];
y_values = [y_values; y_value];
p_values = [p_values; p_value]; 

end


%% MoM direct problem

%stable variables
cd= 2*n;                          %for döngülerinde kolaylık için olan bi variable
c=299792458;                        % ışık hızı
f=300000000;                        %frekans
Lambda= c/f;                        % Lambda
eps0 =(8.8541878178)*(10^(-12));    %epsilon zero
mu0=4*pi*(10^(-7));                   %mu zero   
wn = (2*pi*r)/(2*n) ;               %hücre boyutu
k= 2*pi*f*sqrt(eps0*mu0);          %k aynı k 
nu = sqrt(mu0/eps0);                %nu 0
Y=1.781072418;                      % yot sabiti kitaptaki değer bu idi
x_value = r .* cos(th) + x;
y_value = r .* sin(th) + y;
p_value = sqrt((x_value.^2)+(y_value.^2)); %cylindrical coordinatestaki yarıçap değeri
wn1=sqrt(((x_value(1,2)-x_value(1,1))^2)+((y_value(1,2)-y_value(1,1))^2)); % hücre boyutu alternatif hesap



%üstteki kodla birleştirecek componentler               


Hmnl=[];
Hmnc=[];
Zmnl=[];
Zmnll=[];
Z_mn=[];
R_mnc=[];
R_mnl= [];
ESCATTERED=[];
for c =1:1:dat
    
  for a= 1:1:cd
    
      for b= 1:1:cd
      R_mn= sqrt( (( x_values(c,a) - x_values(c,b) )^2) + ( ( y_values(c,a) - y_values(c,b) )^2)); 
      arg= k*R_mn;                 % hassel inner argumanı
     
      R_mnl=[R_mnl,R_mn];
                % H 'a array ekle
          
          % Direct Funct => Zmn & Hmnl
          if a-b == 0
          Zmnll= ((k*nu*wn1)/4)*(1- (1i*(2/pi)) * (log(((Y*k*wn1)/4))-1) );
          Hmnl1= wn1 - (1i*(2/pi)*wn1)*(log((Y*k*wn1)/4)-1);
          Hmnl=[Hmnl,Hmnl1];
          end

          if not(a-b ==0)
          H_0 = besselh(0,2,arg);      %hassel  
          Zmnll=(k*nu/4)*wn1 *H_0;  % Z'ye array ekle  
          Hmnl=[Hmnl,H_0];
          end
          
      Zmnl=[Zmnl,Zmnll];
      end
  Hmnc=[Hmnc; Hmnl];               % H'a satır ekle
  Z_mn=[Z_mn; Zmnl];               % Z'ye satır ekle
  R_mnc=[R_mnc;R_mnl];
  R_mnl=[];
  Zmnl=[];
  Hmnl=[];

  end 
  
  
%Direct Funcitons => Einc

E_i=transpose(E_inc(x_values(c,:),y_values(c,:),k,0));
E_iabs=abs(E_i);


%Direct Functions => J

J = inv(Z_mn)* E_i;
J_abs= abs(J);


x=0;                % centered at the origin
y=0;
n=50;
th= 0:pi/n:2*pi-0.00001;    % pi 'ın paydasını değiştirince MoMa yaklaşıyoruz
th1= th * 57.2957795131;
ro=5;                %the radius function that has a function in it to generate a random shape
circle(x,y,ro,n);
x_test_value = ro .* cos(th) + x;
y_test_value = ro .* sin(th) + y;
Hmnltest=[];
Hmnctest=[];
R_mnc_test=[];
R_mnl_test= [];

for a= 1:1:cd
    
      for b= 1:1:cd
      R_mn= sqrt( (( x_test_value(1,a) - x_values(c,b) )^2) + ( ( y_test_value(1,a) - y_values(c,b) )^2)); 
      arg= k*R_mn;                 % hassel inner argumanı
     
      R_mnl_test=[R_mnl_test,R_mn];
                % H 'a array ekle
          
          % Direct Funct => Zmn & Hmnl
          if a-b == 0
          Hmnl1= wn1 - (1i*(2/pi)*wn1)*(log((Y*k*wn1)/4)-1);
          Hmnltest=[Hmnltest,Hmnl1];
          end

          if not(a-b ==0)
         
          H_0 = besselh(0,2,arg);      %hassel   
          Hmnltest=[Hmnltest,H_0];
          end
          
      
      end
  Hmnctest=[Hmnctest; Hmnltest];               % H'a satır ekle
  R_mnc_test=[R_mnc_test;R_mnl_test];
  R_mnl_test=[];
  Hmnltest=[];

end 
  

% Scattered Wave Function
Escat=zeros(1,cd);
% Hmnc observation pozisyonuna göre
for ae= 1:1:cd    
Escattered= -(((k*nu)/4)* (J(ae,1)*Hmnctest(ae,:))* wn1);
Escat=  Escat+ Escattered;
end


Escat_abs_MOM=abs(Escat);

ESCATTERED=[ESCATTERED;Escat];

ESCATTERED_REAL=real(ESCATTERED);
ESCATTERED_IMAG=imag(ESCATTERED);
% figure (2)
% plot (Escat_abs_MOM);
Z_mn=[];
end

% %% NOT MOM
% a = 1;
% ro = 5;
% arg1 = 2*pi*a;
% arg2 = 2*pi*ro;
% 
% n = 0;
% Escat = 0;
% Es = [];
% K = 60;
% nm=50; % nokta sayıs
% 
% 
% for fi = 0:(pi/nm):2*pi-0.1
%     in=1;
%     for n = 0:1:K
%             if n == 0
%                 Escat1 = -E_i(in,1)*((-1i)^n) * (besselj(n, arg1) / besselh(n, 2, arg1)) * besselh(n, 2, arg2) * cos(n*fi);
%                 Escat = Escat + Escat1;
%                 
%             end
%             if not(n==0)
%                 Escat2 = -E_i(in,1)*(2*(-1i)^n) * (besselj(n, arg1) / besselh(n, 2, arg1)) * besselh(n, 2, arg2) * cos(n*fi);
%                 Escat = Escat + Escat2;
%                 
%             end 
%   
%     end
%       Es = [Es, Escat];    
%       Escat = 0;  
%       in = in +1;
% end
%     
% Es_abs= abs(Es);
% figure (3)
% 
% plot (Es_abs);
% hold on
% plot(Escat_abs_MOM);
% 
% Diff=[];
%     for b= 1:1:50
%  Difference= (abs((Es_abs(1,b)-Escat_abs_MOM(1,b)))/abs(Es_abs(1,b)));
%  Diff=[Diff,Difference];
%  
%     end
% 
% Error=norm(mean(Diff));
% 
% 
% 
% 
% 


%% functions

function E = E_inc(x,y,k,fi)
E = exp(-1i*k*((x*cos(fi))+(y*sin(fi))));
end
function h = circle(x,y,r,n) %the function to create the circle, x,y = centre coordinates r=radius
hold on
th = 0:pi/n:2*pi-0.000000001;
x_value = r .* cos(th) + x;
y_value = r .* sin(th) + y;
p_value = sqrt((x_value.^2)+(y_value.^2)); %cylindrical coordinatestaki yarıçap değeri
% figure (1)
% h = scatter(x_value, y_value);
% h1= plot(x_value,y_value);
% hold on
end