function Br_fading= spec_chan_derive_fading(fading,delay,DOA,AOA,d_nor,Nr_index,L,M,Nt)

Br_fading_tmp = zeros(M,L,Nt);
for jj = 1 : Nt
    for mm = 1 : M
        for l = 1 : L-1
            Br_fading_tmp(mm,l,jj)=sinc(l-delay(mm,jj))*exp(-1i*2*pi*d_nor*(Nr_index-1)*sin(DOA(mm,jj))*sin(AOA(mm,Nr_index))); %deriveation function false, fixed
        end
    end
end
Br_fading_tmp1=cell(1,Nt);
for jj = 1 : Nt
Br_fading_tmp1{1,jj}=Br_fading_tmp(:,:,jj);
end
Br_fading=blkdiag(Br_fading_tmp1{:});
end

