
function Br_DOA = spec_chan_derive_DOA(fading,delay,DOA,AOA,d_nor,Nr,L,M,Nt)

Br_DOA_tmp = zeros(M,L,Nt);
for jj = 1 : Nt
    for mm = 1 : M
        for l = 1 : L-1
                Br_DOA_tmp(mm,l,jj) =fading(mm,jj)*sinc(l - delay(mm,jj))*(-1i*2*pi*d_nor*(Nr-1)*cos(DOA(mm,jj))*sin(AOA(mm,Nr)))*exp(-1i*2*pi*d_nor*(Nr-1)*sin(DOA(mm,jj))*sin(AOA(mm,Nr)));
        end
    end
end
Br_DOA_tmp1=cell(1,Nt);
for jj = 1 : Nt
Br_DOA_tmp1{1,jj}=Br_DOA_tmp(:,:,jj);
end
Br_DOA=blkdiag(Br_DOA_tmp1{:});
end
