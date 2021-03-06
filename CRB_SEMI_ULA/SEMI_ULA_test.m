clear all;
close all;
clc; 


%%
tic
Nt = 2;    % number of transmit antennas
Nr =8;    % number of receive antennas
L   = 4;    % channel order                     Chua ro
M   = 2;    % Number of multipaths 
Pxp = 10;
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
end

%% Channel generation 

fading=[0.8,0.6,0.4,0.2;0.9,0.7,0.5,0.3];
%delay  = 0.1+(0.15-0.1)*rand(M,Nt)*0.01;
delay=[0.1,0.2,0.3,0.4;0.2,0.3,0.4,0.5];
%DOA    = (-0.5+(0.5-(-0.5))*rand(M,Nt))*pi;
DOA=[pi/2,pi/4,pi/6,pi/8;pi/3,pi/5,pi/7,pi/9];

AOA = pi./[8,5,8,4,4,3,2,4,9,5,4,4,8,3,1,9; 6,2,6,5,1,2,3,1,9,5,9,2,4,4,2,9];
% 
% AOA = zeros(M,16);
% fading = [0.82, 0.77; 0.43, 0.4; 0.89, 0.81; 0.39, 0.76];
% % delay  = rand(M,Nt);
% delay = [0.11,0.19;0.14,0.5;0.68,0.15;0.5,0.05];
% % DOA =pi* rand(M,Nt);
% DOA =[0.26, 2.8; 1.97, 3.09; 2.08, 2.42; 2.29, 1.83];
% % AOA =pi*rand(M,Nr);
% AOA = [2.18 ,0.39, 0.85, 1.31, 0.33, 1.8, 2.32, 3.09, 0.56, 2.95, 1.47, 1.76, 0.17, 2.82, 2.22, 1.46; 1.57, 1.54, 0.650000000000000,0.650000000000000,0.450000000000000,0.160000000000000,0.200000000000000,2.70000000000000,1.25000000000000,0.950000000000000,2.04000000000000,2.68000000000000,0.560000000000000,0.370000000000000,3.14000000000000,2.40000000000000;1.68000000000000,2.68000000000000,1.77000000000000,2.98000000000000,0.520000000000000,2.93000000000000,2.70000000000000,2.47000000000000,0.420000000000000,0.930000000000000,0.0800000000000000,1.09000000000000,2.08000000000000,3.11000000000000,0.900000000000000,2.57000000000000;1.40000000000000,2.75000000000000,2.01000000000000,0.260000000000000,1.95000000000000,2.29000000000000,2.94000000000000,1.61000000000000,0.100000000000000,1.05000000000000,2.65000000000000,1.40000000000000,1.04000000000000,1.70000000000000,1.30000000000000,0.310000000000000];
% % 
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
    X        = [X 1/K*F*diag(ZC(:,ii))*FL];
end


[H, h_true]        = gen_chan_specular(fading,delay,DOA,AOA,Nr,L,Nt);


%% CRB 
% Loop SNR


%============================================
%LAMBDA
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

N_total=4;
N_pilot=N_total - test;
N_data=test;
%============================================
SNR = -15:5:15;
for snr_i = 1 : length(SNR)
    tic
    sigmav2 = 10^(-SNR(snr_i)/10);
%============================================
%Only Pilot    
    X_nga=kron(eye(Nr),X);
%============================================
%Only Pilot Normal
    Iop      = X_nga'*X_nga / sigmav2;
    Iop_4p   = N_pilot * Iop;
    CRB_op(snr_i) = abs(trace(pinv(Iop_4p)));
%============================================
%Only Pilot Specular   
    %Iop_spec = ((-1)/(sigmav2)^2)*G*G'* X_nga'*sigmav2*eye(Nr*K)*X_nga*G*G';
    Iop_spec = G*G'*Iop_4p*G*G';
    CRB_op_spec(snr_i) = abs(trace(pinv(Iop_spec)));
%============================================
%SemiBlind
    Cyy      = sigmax2 * LAMBDA * LAMBDA'  + sigmav2 * eye(K*Nr);
    Cyy_inv  = pinv(Cyy);
    paths = Nr*Nt*L;
%     partial_Cyy_hii = zeros(Nr*K , Nr*K, paths);
%     for ii = 1 : paths
%         partial_Cyy_hii(:,:,ii) = sigmax2 * LAMBDA * partial_LAMBDA{1,ii}';
%     end
    
    for ii =1: paths
          partial_Cyy_hii = sigmax2 * LAMBDA * partial_LAMBDA{1,ii}';
        for jj = ii : paths
            partial_Cyy_hjj =sigmax2 * LAMBDA * partial_LAMBDA{1,jj}';
            % Slepian-Bangs Formula
            I_D(ii,jj) = trace(Cyy_inv * partial_Cyy_hii* Cyy_inv * partial_Cyy_hjj);
            I_D(ii,jj) = I_D(ii,jj)';
        end
    end
%============================================
%Semiblind Normal
    I_SB       = N_data* I_D+N_pilot*Iop;
    CRB_SB_i           = pinv(I_SB);
    CRB_SB(snr_i)      = abs(trace(CRB_SB_i));
    
%============================================
%Semiblind Specular
   I_SB_spec=G*G'*I_SB*G*G';
   CRB_SB_spec_i           = pinv(I_SB_spec);
   CRB_SB_spec(snr_i)      = abs(trace(CRB_SB_spec_i));
   log = [num2str(snr_i), ' times.'];
   disp(log);
   toc
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