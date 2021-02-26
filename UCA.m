clear all
close all
clc

N = 10;% so phan tu 
kr =1.3*pi;
n = 1:N;
%%%Thong so ban dau
% Teta = 90;
Omega_n = (2*pi*(n-1/2))/N;%
Omega =[0:360]*pi/180;% goc quay phuong vi
% phuong trinh AF
u = kr*cos(Omega - Omega_n.');
a = exp(j*u);
%uniform
w_u = ones(1,N);
AF_u =w_u*a;

%PHASE ARRAY
Omega_choose =108*2*pi/360;
beta = kr*cos(Omega_choose-Omega_n);
w_p = exp(-j*beta);
AF_p = w_p*a;

%pha2
Omega_ch2 = 0*2*pi/360;
beta_2 = kr*cos(Omega_ch2-Omega_n);
w_p2 = exp(-j*beta_2);
AF_p2=w_p2*a;
%Taylor
w_t = [ 0.76955  0.97473  1.845486  1.845486 0.97473 0.76955 0.02409 0.37537 0.37537 0.02409 ];
AF_t = w_t.*w_p*a;
%ga
%w1 = [0.460	1.945	0.145	1.201	0.779 1.913	1.034	0.230	1.551	1.256]
w2 = [1.603	0.934	0.995	1.103	1.686 0.615	0.635	1.062	0.245	0.763]
%AFu_ga = w1.*w_p*a;
AFt_ga = w2.*w_p*a;
%pso
wt =[0.882253396252479,1.351137556741155,2.616779112858612,1.440375969169201,1.582463526357811,1.474906907948593,0.944182888551725,0.032167582109034,0.518751630939649,-0.752364945640637]
%wu = [0.480771261416160,1.287046025001676,1.225007806305886,1.876854659518131,0.018098237320116,1.050280851205014,1.485809611646276,1.017435476635856,0,0.552334388183211]

AFt_pso = wt.*w_p*a;
%AFu_pso = wu.*w_p*a;
%plot
figure(1)
plot(Omega*180/pi,(abs(AF_p/max(AF_p))),'--');
hold on;
plot(Omega*180/pi,(abs(AF_t/max(AF_t))));
hold on;
legend('Uniform','Taylor')

figure(2)
plot(Omega*180/pi,(abs(AFt_ga/max(AFt_ga))));
hold on;
plot(Omega*180/pi,(abs(AFt_pso/max(AFt_pso))),'--');
hold on;
plot(Omega*180/pi,(abs(AF_p/max(AF_p))),'-.');
hold on;
legend('GA Sym','Pso Symmetry','Uniform');

