clear all 
clc
close all
tic
N_t = 4; %s? antenna phát
N_r =2; %s? antenna thu
L = 5; %channel tap ?  hay là maximum channel delay
M=L;
Pxp = 10; 
Pxd = [90 92 88 86];
%Pxt = 100;
Ns_data = 20;
Np = 1;
N_OFDM = 64;

W = dftmtx(N_OFDM);

F = W(:,1:L);
sigmas2 = Pxd;
sigmas = sqrt(sigmas2);

U = [1 3 5 7];
ZC_p = []; ZC_d = [];
for u = 1:N_t
    for k =1 :N_OFDM
        ZC(k,u) = sqrt(Pxp)*exp((-1i*pi*U(u)*(k-1)^2)/N_OFDM);
        ZCd(k,u) = sqrt(Pxd(u))*exp((-1i*pi*U(u)*(k-1)^2)/N_OFDM);
    end
    ZC_p = [ZC_p;ZC(:,u)];
    ZC_d = [ZC_d;ZCd(:,u)];
end

%test_listing = zeros(1,100);

%for test = 1:100

fading = rand(L,N_t);
delay = rand(L,N_t);
DOA = pi*rand(L,N_r)/6; %theo ma tr?n thu
AOA = pi*rand(L,N_t)/6; %theo ma tr?n phát 

H = spec_chan(fading, delay, DOA,AOA,N_r,L,N_t);

Cyy_H = zeros(N_OFDM*N_r);
Cyy_T = zeros(N_OFDM*N_r);
Cyy_lambda1 = zeros(N_r);

lambda3 =[];
lambda4 = [];
h = []; h_ifft = [];
for ii = 1:N_t
    h = [h reshape(transpose(H(:,:,ii)),1,N_r*L)];
    h_ifft = [h_ifft reshape(fft(transpose(H(:,:,ii))),1,N_r*L)];
end

hmod2 = (h*h')^2;
hmod2_ifft = (h_ifft*h_ifft')^2;

M1 = cell(1,N_t);
for ii = 1:N_t
    M_l = diag(1./fading(:,ii));
    MM_l = repmat(M_l,N_r);
    M1{1,ii} = MM_l;
end
Mat1 = blkdiag(M1{:});

D1_fading = repmat(h,size(h,2),1);
D_fading = Mat1 .* D1_fading;

Hderiv = spec_chan_derive_sinc(fading,delay,DOA,AOA,N_r,L,N_t);
hnew = [];
 for ii = 1:N_t
     Hderiv_i = Hderiv(:,:,ii);
     hnew = [hnew Hderiv_i(:)'];
 end
 
  M2 = cell(1,N_t);
  for ii = 1:N_t
      M2{1,ii} = repmat(eye(L),N_r);
  end
  Mat2 = blkdiag(M2{:});
  D_delay = repmat(hnew,size(hnew,2),1);
  D_delay = Mat2.*D_delay;
  
    M3 = cell(1,N_t);
  for ii = 1:N_t
      Mtetal = diag( - 1i*pi.*cos(AOA(:,ii)));
      Mtetall = repmat(Mtetal,[N_r 1]);
      M3_l = zeros(L*N_r);
      M3_l(:,end-L+1:end) = Mtetall;
      M3{1,ii} = M3_l;
  end
  Mat3 = blkdiag(M3{:});
  
  D_AOA = repmat(h,size(h,2),1);
  D_AOA = Mat3 .* D_AOA;

  M4 = cell(1,N_r);
  for ii = 1:N_r
      Mtetal = diag( - 1i*pi.*cos(DOA(:,ii)));
      Mtetall = repmat(Mtetal,[N_t 1]);
      M4_l = zeros(L*N_t);
      M4_l(:,end-L+1:end) = Mtetall;
      M4{1,ii} = M4_l;
  end
  Mat4 = blkdiag(M4{:});
  
  D_DOA = repmat(h,size(h,2),1);
  D_DOA = Mat4 .* D_DOA;
  
  DD = [D_fading D_delay D_AOA D_DOA];
  
      
  
  Xt = [];
  Xt_lambda = [];
  for ii = 1:N_t
      Xt = [Xt (W'/sqrt(N_OFDM))*diag(ZC(:,ii))*F];
      Xt_lambda = [Xt_lambda sqrt(N_OFDM) * diag(ZC(:,ii))];
  end
  XX_h = Xt' *Xt;

  G = DD*DD';
  X_i = kron(eye(N_r),Xt);
  
  SNR = -10:5:20;
  for snr_i = 1:length(SNR)
      sigmav2 = 10^(-SNR(snr_i)/10);
      Ip1 = Np * X_i'*X_i/sigmav2;
      Iop_spec = - Np * G' *X_i'*eye(N_r*N_OFDM)*X_i *G/sigmav2;
      CRB_op1(snr_i) = abs(trace(inv(Ip1)));
      CRB_op_spec(snr_i) = abs(trace(pinv(Iop_spec)));
  end

  
%end

test_listing = CRB_op1(5)-CRB_op_spec(5);
  figure
  semilogy(SNR,CRB_op1,'-b>')
  hold on; semilogy(SNR,CRB_op_spec,'-r+')
  grid on
  ylabel('Normalized CRB')
  xlabel('SNR(dB)')
  legend('ususal','specular')
  axis([-10 20 1e-4 1e-2])
  title('Not apple to apple')
    