clear all;
close all;
clc; 


%%
tic
Nt = 4;    % number of transmit antennas
Nr = 4;    % number of receive antennas
L   = 5;    % channel order
M   = L;    % Number of multipaths (assumption: M  = L)   
Pxp = 50;
Pxd = [90 92 88 86];
Pxt = 100;
K  = 64;        % OFDM subcarriers
Np      = 4;      % Pilot symbols
Ns_data = 48;      % Data symbols
F  = dftmtx(K);
FL  = F(:,1:L);
d_nor=1/2;
 
sigma_x = 0.01;
%% Channel generation 
% Fading, delay, DOA matrix of size(M,Nt), M - the number of multipath
%fading = rand(M,Nt)+1i*rand(M,Nt);
fading = 0.8 + (0.9-0.8).*rand(M,Nt);
% fading=[0.8,0.6,0.4,0.2;0.9,0.7,0.5,0.3];
delay  = 0.00001 + (0.0002-0.00001).*rand(M,Nt);
% delay=[0.1,0.2,0.3,0.4;0.2,0.3,0.4,0.5];
% DOA = pi * rand(M,Nt);
DOA = pi * zeros(M,Nt);
% DOA=[pi/2,pi/4,pi/6,pi/8;pi/3,pi/5,pi/7,pi/9];
% AOA = pi*rand(M,Nr);
AOA = pi * zeros(M,Nr);
% AOA =[pi/2,pi/4,pi/6,pi/8;pi/3,pi/5,pi/7,pi/9];

H = spec_chan(fading,delay,DOA,AOA,Nr,L,M,Nt);

% fading = ones(M,Nt);
% delay  = zeros(M,Nt);
% DOA=ones(M,Nt);
% AOA = DOA;

%% Derivative
% w.r.t. fading
%Br_fading = zeros(Nt,M,L);
dev_h_fading_tmp=[];
dev_h_delay = [];
dev_h_DOA = [];
%dev_h_fading_tmp=cell(Nr,1);

dev_h_fading=[];
for Nr_index=1:Nr
D_fading = spec_chan_derive_fading(fading,delay,DOA,AOA,d_nor,Nr_index,L,M,Nt);
dev_h_fading=[dev_h_fading; transpose(D_fading)];
D_delay =  spec_chan_derive_delay(fading,delay,DOA,AOA,d_nor,Nr_index,L,M,Nt);
dev_h_delay = [dev_h_delay; transpose(D_delay)];
D_DOA =  spec_chan_derive_DOA(fading,delay,DOA,AOA,d_nor,Nr_index,L,M,Nt);
dev_h_DOA = [dev_h_DOA; transpose(D_DOA)];
end

D_AOA = spec_chan_derive_AOA(fading,delay, DOA,AOA, Nr, L,M, Nt);

%% Derivation of $h$ w.r.t. (bar{h},tau,a4lpha) %% channel specular parameters

G = [dev_h_fading dev_h_delay dev_h_DOA D_AOA]; 

% G = dev_h_fading+ dev_h_delay+ dev_h_DOA+ D_AOA;
%% channel matrix
LAMBDA  = [];
for jj = 1 : Nt
    lambda_j =[];
    for r = 1 : Nr
        lambda_rj  = diag(FL*H(r,:,jj)');
        lambda_j   = [lambda_j; lambda_rj];
    end
    LAMBDA = [LAMBDA lambda_j];
end

 partial_LAMBDA  = cell(1,Nr*Nt*L);
for idx = 1:Nr*Nt*L
    for ll = 1 : L
        partial_LAMBDA_ll = [];
        for jj = 1 : Nt
            lambda_jj =[];
            for r = 1 : Nr
                lambda_rj_ll = diag(FL(:,ll));
                lambda_jj    = [lambda_jj; lambda_rj_ll];
            end
            partial_LAMBDA_ll = [partial_LAMBDA_ll lambda_jj];
        end
        partial_LAMBDA{1,idx} = partial_LAMBDA_ll;
    end
end
%% Signal Generation
% we use the Zadoff-Chu sequences
U = 1:2:7;
ZC_p = [];
ZC_d = [];
for u = 1 : Nt
    for k = 1 : K
        ZC(k,u) = sqrt(Pxp) * exp( ( -1i * pi * (u*2-1) * (k-1)^2 ) / K );
        ZCd(k,u) = sqrt(Pxd(u))*exp((-1i*pi*(u*2-1)*(k-1)^2)/K);
    end
    ZC_p = [ZC_p; ZC(:,u)];
    ZC_d = [ZC_d;ZCd(:,u)];
end
%% 

X = [];
for ii = 1 : Nt
    X        = [X diag(ZC(:,ii))*FL];
end

%% CRB 

D = G*G';

 %% Loop SNR
SNR = -15:10:15;%%dBm
CRB_op = zeros(length(SNR),1);
CRB_op_spec = zeros(length(SNR),1);
CRB_OD = zeros(length(SNR),1);
CRB_OD_spec = zeros(length(SNR),1);
SB = zeros(length(SNR),1);
CRB_SB_spec = zeros(length(SNR),1);
for snr_i = 1 : length(SNR)
    
    sigmav2 = 10^(-SNR(snr_i)/10);

    Cyy      = sigma_x * LAMBDA * LAMBDA'  + sigmav2 * eye(Nr*K);
    Cyy_inv  = pinv(Cyy);
    
   I_D = zeros(Nr*Nt*L);
    for ii = 1 : Nr*Nt*L
        partial_Cyy_hii = sigma_x * LAMBDA * partial_LAMBDA{1,ii}';
        for jj = ii : Nr*Nt*L
            partial_Cyy_hjj = sigma_x * LAMBDA * partial_LAMBDA{1,jj}';
            % Slepian-Bangs Formula
            I_D(ii,jj) = trace(Cyy_inv * partial_Cyy_hii * Cyy_inv * partial_Cyy_hjj);
            I_D(ii,jj) = I_D(ii,jj)';
        end
    end
    
    X_nga=kron(X,eye(Nr));
    Iop      = Np*X_nga'*X_nga / sigmav2;

    Iop_spec = Np*(-1/sigmav2)*D*( X_nga'*eye(Nr*K)*X_nga)*D';
    
    Iop_SB_usual =Ns_data * I_D+Iop;
    Iop_SB_spec =Ns_data*D*I_D*D' + Iop_spec;
   
    CRB_op(snr_i) = abs(trace(pinv(Iop)));
    CRB_op_spec(snr_i) = abs(trace(pinv(Iop_spec)));
    CRB_OD(snr_i) = abs(trace(pinv(I_D)));
    CRB_OD_spec(snr_i) = abs(trace(pinv(D*I_D*D')));
%     SB(snr_i) = abs(trace(inv(Iop_SB_usual)));
%     CRB_SB_spec(snr_i) = abs(trace(pinv(Iop_SB_spec)));
    
end


    figure
    semilogy(SNR,CRB_op,'-b>');
    hold on; semilogy(SNR,CRB_op_spec,'-r>');
     hold on; semilogy(SNR, CRB_OD,'-g');
     hold on; semilogy(SNR, CRB_OD_spec,'-m');
%      hold on; semilogy(SNR, SB,'--bx');
%      hold on; semilogy(SNR,CRB_SB_spec,'--ro'); 
    
    grid on
    ylabel('Normalized CRB')
    xlabel('SNR(dB)')
    legend('CRB_{pilot usual}','CRB_{pilot specular}','CRB_{data CRB}','CRB_{data but specular}','CRB_{semiblind usual}','CRB_{semiblind specular}')
    title(' ')
    %axis([-10 20 1e-4 1e2])
% % end
% end