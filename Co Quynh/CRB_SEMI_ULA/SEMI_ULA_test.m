clear all;
close all;
clc; 


%%
tic
Nt = 2;    % number of transmit antennas
Nr = 16;    % number of receive antennas
L   = 4;    % channel order                     Chua ro
M   = 2;    % Number of multipaths 
Pxp = 10;
Pxd = [90 92 88 86];
Pxt = 100;
K  = 64;        % OFDM subcarriers
F  = dftmtx(K);
FL  = F(:,1:L);
sigmax2=1;      % Phuong sai cua du lieu

%% Signal Generation
% we use the Zadoff-Chu sequences
U = 1:2:7;
ZC_p = [];
for u = 1 : Nt
    for k = 1 : K
        ZC(k,u) = sqrt(Pxp) * exp( ( -1i * pi * U(u) * (k-1)^2 ) / K );
    end
    ZC_p = [ZC_p; ZC(:,u)];
end

%% Channel generation 
% Fading, delay, DOA matrix of size(M,Nt), M - the number of multipath
%fading = rand(M,Nt)+1i*rand(M,Nt);
%fading = 0.1+(0.6-0.1)*rand(M,Nt);
fading=[0.8,0.6,0.4,0.2;0.9,0.7,0.5,0.3];
%delay  = 0.1+(0.15-0.1)*rand(M,Nt)*0.01;
delay=[0.1,0.2,0.3,0.4;0.2,0.3,0.4,0.5]*0.001;
%DOA    = (-0.5+(0.5-(-0.5))*rand(M,Nt))*pi;
DOA=[pi/2,pi/4,pi/6,pi/8;pi/3,pi/5,pi/7,pi/9];

AOA = pi./[8,5,8,4,4,3,2,4,9,5,4,4,8,3,1,9; 6,2,6,5,1,2,3,1,9,5,9,2,4,4,2,9];
% AOA = zeros(M,16);
d_nor=1/2;

%H      = spec_chan(fading,delay,DOA,Nr,L,Nt);  % H(Nr,L,Nt)

%h=[];
%for ii = 1 : Nt
%    h              = [h reshape(transpose(H(:,:,ii)),1,Nr*L)];
%end

%hmod2      = (h*h')^2; 
%% Derivative
% w.r.t. fading
%Br_fading = zeros(Nt,M,L);
dev_h_fading_tmp=[];
dev_h_delay_tmp=[];
dev_h_angle_tmp=[];
%dev_h_fading_tmp=cell(Nr,1);

for test = 1:3
    
dev_h_fading=[];
dev_h_delay=[];
dev_h_angle=[];
dev_h_AOA = [];
for Nr_index=1:Nr
Br_fading = SEMI_spec_chan_derive_fading_ULA(fading,delay,DOA,AOA,d_nor,Nr_index,L,M,Nt);
dev_h_fading=[dev_h_fading; transpose(Br_fading)];

Br_delay = SEMI_spec_chan_derive_delay_ULA(fading,delay,DOA,AOA,d_nor,Nr_index,L,M,Nt);
dev_h_delay=[dev_h_delay; transpose(Br_delay)];

Br_angle = SEMI_spec_chan_derive_angle_ULA(fading,delay,DOA,AOA,d_nor,Nr_index,L,M,Nt);
dev_h_angle=[dev_h_angle; transpose(Br_angle)];

Br_arive_angle = SEMI_spec_chan_derive_arive_angle_ULA(fading,delay,DOA,AOA,d_nor,Nr_index,L,M,Nt);
dev_h_AOA=[dev_h_AOA; transpose(Br_arive_angle)];
%dev_h_fading=[dev_h_fading;Br_fading];
end

%% Derivation of $h$ w.r.t. (bar{h},tau,alpha) %% channel specular parameters

G = [dev_h_fading,dev_h_delay,dev_h_angle, dev_h_AOA]; 
%% ------------------------------------------------------------------------
 
X = [];
for ii = 1 : Nt
    X        = [X diag(ZC(:,ii))*FL];
end


[H, h_true]        = gen_chan_specular(fading,delay,DOA,AOA,Nr,L,Nt);


%% CRB 
% Loop SNR


%============================================
%LAMBDA                            Chua hieu lam
LAMBDA  = [];
for jj = 1 : Nt
    lambda_j =[];
    for r = 1 : Nr
        h_rj       = transpose(H(r,:,jj));
        lambda_rj  = diag(FL*h_rj);
        lambda_j   = [lambda_j; lambda_rj];
    end
    LAMBDA = [LAMBDA lambda_j];
end


% Partial derivative of LAMBDA w.r.t. h_i                  
% Dao ham Cyy theo h => Lambda theo h truoc
partial_LAMBDA  = cell(1,Nr*Nt*L);
idx = 1;
while  idx <= Nr*Nt*L
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
        idx = idx + 1;
    end
end


% for test = 0:3
N_total=4;
N_pilot=N_total - test;
N_data=test; % Dung 1 phan data
%============================================
SNR = -10:5:30;
for snr_i = 1 : length(SNR)
    sigmav2 = 10^(-SNR(snr_i)/10);
%============================================
%Only Pilot    
    X_nga=kron(eye(Nr),X);
%============================================
%Only Pilot Normal
    Iop      = N_pilot * X_nga'*X_nga / sigmav2;            % FIM
    Iop_wp      = N_pilot *Iop;
    CRB_op(snr_i) = abs(trace(pinv(Iop_wp)));
%============================================
%Only Pilot Specular   
    %Iop_spec = ((-1)/(sigmav2)^2)*G*G'* X_nga'*sigmav2*eye(Nr*K)*X_nga*G*G';
    Iop_spec = G*G'*Iop_wp*G*G';
    CRB_op_spec(snr_i) = abs(trace(pinv(Iop_spec)));
%============================================
%SemiBlind
    Cyy      = sigmax2 * LAMBDA * LAMBDA'  + sigmav2 * eye(K*Nr);
    Cyy_inv  = pinv(Cyy);   
       
    for ii = 1 : Nr*Nt*L
        partial_Cyy_hii = sigmax2 * LAMBDA * partial_LAMBDA{1,ii}';
        for jj = ii : Nr*Nt*L
            partial_Cyy_hjj = sigmax2 * LAMBDA * partial_LAMBDA{1,jj}';
            % Slepian-Bangs Formula
            I_D(ii,jj) = trace(Cyy_inv * partial_Cyy_hii * Cyy_inv * partial_Cyy_hjj);
            I_D(ii,jj) = I_D(ii,jj)';
        end
        disp(ii);
    end
%============================================
%Semiblind Normal
    I_SB       = N_data*I_D+N_pilot*Iop;
    CRB_SB_i           = pinv(I_SB);
    CRB_SB(snr_i)      = abs(trace(CRB_SB_i));
    
%============================================
%Semiblind Specular
   I_SB_spec=G*G'*I_SB*G*G';
   CRB_SB_spec_i           = pinv(I_SB_spec);
   CRB_SB_spec(snr_i)      = abs(trace(CRB_SB_spec_i));  
end
% 
% figure
C = {'b','r',[0, 0.5, 0]};
Color = {'k',[0.4940, 0.1840, 0.5560],'m'};
semilogy(SNR,CRB_op,'color',C{test},'marker','>');
hold on; semilogy(SNR,CRB_op_spec,'color',Color{test},'marker','+');

semilogy(SNR,CRB_SB,'color',C{test},'marker','x');
hold on; semilogy(SNR,CRB_SB_spec,'color',Color{test},'marker','o')
end
grid on
ylabel('Normalized CRB')
xlabel('SNR(dB)')
legend('normal OP 3 pilot','spec OP 3 pilot','normal SB 3 pilot 1 data','spec SB 3 pilot 1 data','normal OP 2 pilot','spec OP 2 pilot','normal SB 2 pilot 2 data','spec SB 2 pilot 2 data','normal OP 1 pilot','spec OP 1 pilot','normal SB 1 pilot 3 data','spec SB 1 pilot 3 data')
%legend('normal_OP','spec_OP','normal_SB','spec_SB')
title(' ')
%axis([-10 20 1e-4 1e2])
hold on;
