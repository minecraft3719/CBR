clc;
clear all;
close all;

Nt = 2; %transmitter antenna
Nr = 2; %receiver antenna

K = 64; %subcarrier

F = dftmtx(K) * 1/sqrt(K); %DFT matrix

Np = 1; %pilot simbol
 
Nd = 40; %data signal simbol
CP = K - (Np + Nd + 1) ;

W = F(1,1:CP);


