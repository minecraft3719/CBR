clear all
close all
clc

N = 10;% so phan tu la 10
k_nor = 2*pi%/lamda=1 ;% tinh k=2*pi/lamda
n = 1:N;
Omega_n = (2*pi*(n-1))/N;%
Omega =[0:360]*2*pi/360;% goc quay phuong vi
% phuong trinh AF
u = k_nor*cos(Omega - Omega_n.');%sin(Teta_0) = 1
a = exp(j*u);

%PHASE ARRAY
Omega_pack = 120*2*pi/360;
Teta_pack = 90*2*pi/360;
beta = k_nor*sin(Teta_pack)*cos(Omega_pack - Omega_n);
w_p = exp(-j*beta);
AF_p = w_p*a;
%phase null
Omega = 68*2*pi/360;
Teta_n = 90*2*pi/360;
beta_n = k_nor*sin(Teta_n)*cos(Omega-Omega_n);
w_n = exp(j*beta_n);
w_p_null = w_p.*w_n;
%AF_n = w_n*a;
%Ga--zero
Delta_x = 1;
x_c = ones(1,N);
%fun =@(w) (abs(w(1)*w_p(1)+w(2)*w_p(2)+w(3)*w_p(3)+w(4)*w_p(4)+w(5)*w_p(5)+w(6)*w_p(6)+w(7)*w_p(7)+w(8)*w_p(8)+w(9)*w_p(9)+w(10)*w_p(10)));
fun =@(w) (abs(w(1)*w_p_null(1)+w(2)*w_p_null(2)+w(3)*w_p_null(3)+w(4)*w_p_null(4)+w(5)*w_p_null(5)+w(6)*w_p_null(6)+w(7)*w_p_null(7)+w(8)*w_p_null(8)+w(9)*w_p_null(9)+w(10)*w_p_null(10)));
x = ga(fun,N,[],[],[],[],x_c-Delta_x,x_c+Delta_x);
AF_g = x.*w_p*a;
%PSO

x2 = particleswarm(fun,N,x_c-Delta_x,x_c+Delta_x);
AF_pso = x2.*w_p*a;
%plot
figure(1)
plot(Omega*180/pi,abs((AF_p)/max(AF_p)));
hold on
figure(2)
plot(Omega*180/pi,abs((AF_g)/max(AF_g)));
plot(Omega*180/pi,abs((AF_pso)/max(AF_pso)));
legend('GA','PSO');
hold on



