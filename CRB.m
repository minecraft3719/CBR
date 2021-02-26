close all;
clear all;
clc;

M = 4;%Number of path
Nr = 2;% receive antenna
DOA = [pi/2 pi/4 pi/6 pi/8; pi/3 pi/4 pi/7 pi/8]; %DOA(each path have different DOA);
lambda = 0.42; %wavelength
d = 1e3;%length between two antenna
a = zeros(Nr,M);
for i = 1:Nr
    for l = 1:M
        a(i,l) =-j*2*pi*d*(i-1).*cos(DOA(i,l))/lambda;
        a_steer = exp(a);
    end
end

t_delay = [ 0.4,0.6,0.1,0.01;0.3,0.9,0.5,0.3];

h_fading = [0.4 0.6 0.1 0.01;0.3 0.9 0.5 0.3]; %for l path

t = 10; %assume time transmit = 1

h = zeros(Nr);
h_t = zeros(Nr,t); %impulse response of each antenna by time t
for i = 1 : Nr
   for l = 1:M
      h(i,l) = h_fading(i,l) .* a_steer(i,l) .* sinc(t-t_delay(i,l)); 
   end
end

h_t = sum(h,2);

%%after calculate the fading, we have channel coefficient

K = 64; %subcarrier

F = dftmtx(K) * 1/sqrt(K); %DFT matrix

Np = 1; %pilot simbol
 
Nd = 40; %data signal simbol
CP = K - (Np + Nd + 1) ;

W = F(:,1:CP);


%%calculate theta
for i = 1:Nr
    LAMBDA(i,:,:) = (h_t(i) .* W).';
    R_LAMBDA(i,:,:) = ((h_t(i) .* W)').';
end

 theta = LAMBDA * R_LAMBDA;



