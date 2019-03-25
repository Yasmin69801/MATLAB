% This code was written to solve for the six coupled ODEs that describe the problem of self-association.
% Ref: Yogurtcu, Osman N., and Margaret E. Johnson. "Cytosolic proteins can exploit membrane localization to trigger functional assembly." PLoS computational biology 14.3 (2018): e1006031
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%function to solve for the self-association problem, and to output the K
% effective in order to measure equilibrium and enhancement.
%kmain=((1/kapp)+(1/(4*s*12000)))^-1

function []=odesolverdimer()
 clear;  
 %close all; 
 clc;
 %species we will be solving for:
% p
% m
% pm
% ppm
% mppm
% pp

%initial conditions:
tspan=[0,1];      %in s
P0= 16.87;     %uM  
M0= 136.33;    %uM  
%needed parameters
kapp=0.214;       %for D=12            %(uM.s)^-1                   %0.3556 nm^3/us                
kfpp= 1.0;        %s^-1                %1.074 ; %s^-1
kapm=40.3;        %40.3;               %(uM.s)^-1                   %67.0    %nm^3/us    
kfpm= 1.0;        %s^-1                %1.226 ; %s^-1
s=1.00E-03 ;      %sigma (um)
V=0.9993;         %in um^3             %9.9933e-16 in Liters        % V = (1811.0*(10^-9))*(1811.0*(10^-9))*(304.7*(10^-9))
A=3.2797;         %in um^2                                          % A= (1811.0*(10^-9))*(1811.0*(10^-9))
g=152.5;          %V/(2*A*s) unitless


%coefficents matrix:
theta=[kapp,kfpp,kapm,kfpm,s,g];

Y0 = zeros(4,1);         %designing an empty list of 6 items,solving for 6 species
Y0(1)=P0;                %sets the first item as the intital amount of Protein
Y0(2)=M0;               %sets the second item as the intital amount of Membrane lipid
Y0(3)=0;
Y0(4)=0;
% Y0(5)=0;
% Y0(6)=0;

disp(Y0)
%Y0: (P0,M0,0,0,0,0)

%%%%%%%%%Solving the ODE %%%%%%%%%
%arguments for ode23s : function, time span,initial condition, tolerence, coeeficients,options
%opt=odeset('RelTol',1E-4,'AbsTol',1E-7);
[time, y]  = ode23s(@ode_sys, tspan, Y0, odeset('MaxStep',1e-5), theta) ;


%Plotting the 6 species 

Plotting the 6 species 

f1=figure;
figure(f1);
semilogx(time,y(:,1));
title('P');
xlabel('Time in seconds')
ylabel('Concentration in uM')
legend('P');
hold on

f2=figure;
figure(f2);
semilogx(time,y(:,2));
hold on;
title('M');
xlabel('Time in seconds')
ylabel('Concentration in uM')
legend('M');
hold on

f3=figure;
figure(f3);
semilogx(time,y(:,3));
hold on;
title('PM');
xlabel('Time in seconds')
ylabel('Concentration in uM')
legend('PM');

f4=figure;
figure(f4);
semilogx(time,y(:,4));
hold on;
title('PP');
xlabel('Time in seconds')
ylabel('Concentration in uM')
legend('PP');


f5=figure;
figure(f5);
semilogx(time,y(:,5));
hold on;
title('PPM');
xlabel('Time in seconds')
ylabel('Concentration in uM')
legend('PPM');

f6=figure;
figure(f6);
semilogx(time,y(:,6));
hold on;
title('MPPM');
xlabel('Time in seconds')
ylabel('Concentration in uM')
legend('MPPM');

%storing the info in each specific species vector 
p=y(end,1);
m=y(end,2);
pm=y(end,3);
pp=y(end,4);
ppm=y(end,5);
mppm=y(end,6);

Ktheory = kapp*((g*kapm*kapm*M0*M0+2*kapm*M0+1))/;

Species={'P';'M';'PM';'PP';'PPM';'MPPM';'Time'};
Equilibrium_Value = [p;m;pm;pp;ppm;mppm;time(end)];
Initial_Values = [16.87;136.33;0;0;0;0;0];
T=table(Species,Equilibrium_Value,Initial_Values)


%the equilibrium we were looking for and where trying to identify
keff=(PPsol+PPmembrane)/((Psol+Pmem)^2)

Kaeff= (mppm+ppm+pp)/(p+pm)^2;
calcEnh = Kaeff/kapp;
enhancement = (g*kapm*kapm*m*m + (kapm+kapm)*m +1) /((1+kapm*m)^2);
fprintf('The effective equilibrium constant is %4.2f and the enhancement is %4.2f \n',Kaeff,enhancement)

end

%%%%%%%%% the ODE system function %%%%%%%%%

function dydt= ode_sys(t,y,theta)
%theta=[kapp,kfpp,kapm,kfpm,V,A,s,g];

kapp=theta(1);
kfpp=theta(2);
kapm=theta(3);
kfpm=theta(4);
s=theta(5);
g=theta(6);


p=y(1);
m=y(2);
pm=y(3);
pp=y(4);
%  ppm=y(5);
% mppm=y(6);



dydt(1,1) = -kapp*p^2+kfpp*pp-kapm*p*m+kfpm*pm;      %P
dydt(2,1) = -kapm*p*m+kfpm*pm;  %M
dydt(3,1) = 2*kapm*p*m-2*kfpm*pm; %PM
dydt(4,1) = kapp*p^2-kfpp*pp;                            %PP
% dydt(5,1) = 2*kapp*pm*p-kfpp*ppm+2*kapm*pp*m-kfpm*ppm-g*kapm*ppm*m+2*kfpm*mppm;%PPM
% dydt(6,1) = g*kapp*pm^2-kfpp*mppm+kapm*g*ppm*m-2*kfpm*mppm;                   %MPPM



return;

end
