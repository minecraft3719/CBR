clc;
clear all;
close all;

Nt = 2; %number of transmitting antenna
Nr = 2; %number of receiving antenna
N = 4; %number of channel tap

M = 4; %number of channel path

K = 64; %number of OFDM subcarrier
Np = 1; %number of OFDM pilot signal
Nd = 40;%number of OFDM data signal

Pxp = [23 13]; %pilot signal power dBm
Pxd = [20 18.8]; %data signal power dBm

CP = K - (Nt + Nd + 1);
F = dftmtx(K);

%W = F(: , 1:CP ); if N = CP
W = F(:,1:N); %if N = channel tap;
%----------------------------------------------------------------------------------

fading = [0.4 0.6 0.1 0.01 ; 0.3 0.9 0.5 0.3].';
delay = [ 0.4 0.6 0.1 0.4; 0.3 0.9 0.5 0.1].';
DOA = pi./[2 4 6 8;3 4 7 8].';

% calculate H----------------------------------------------------------------------------------
H = zeros(Nr,N,Nt);
for i = 1:Nt
    for j = 1 :Nr
        for k = 1:N
            H(j,k,i) = fading(k,i) *sinc(delay(k,i)) *exp(-1i*pi*(j-1)*sin(DOA(k,i)));
        end
    end
end


%%lambda code 
Cyy_l = zeros(K*Nr);
Cyy_N = zeros(K*Nr);
Cyy_lambda1 = zeros(Nr);

h=[];h_ifft=[];
lambda3=[];
lambda4=[];

for ii = 1:Nt
    h =[h reshape(transpose(H(:,:,ii)),1,Nr*N)];
    lambda(:,:,ii) = transpose(fft(transpose(H(:,:,ii)),K));
    h_ifft = [h_ifft reshape(fft(transpose(H(:,:,ii))),1,Nr*N)];
    lambda_i = [];
    lambda_2 = [];

    for jj = 1:Nr
        lambda_i = [lambda_i; (W'./sqrt(K)) * diag(F*transpose([H(ii,:,jj) zeros(1,K-N)]))];
        lambda_2 = [lambda_2;F(20,:)*transpose(H(jj,:,ii))];
    end

    lambda3 = [lambda3 lambda_i];
    lambda4 = [lambda4 lambda_2];
    lambda2(:,:,ii) = lambda_i;
    temp_lambd(:,:,ii)=lambda_i*lambda_i;
    temp_lambd2(:,:,ii)=lambda_2*lambda_2;
    temp_lambd_t(:,:,ii)=lambda_i*transpose(lambda_i);
    Cyy_l = Cyy_1+temp_lambd(:,:,ii)*sigmas2(ii);
    Cyy_lambda1=Cyy_lambda1+temp_lambd2(:,:,ii)*sigmas2(ii);
    Cyy_Nc= Cyy_Nc+temp_lambd_t(:,:,ii)*Nc_P(ii);
end