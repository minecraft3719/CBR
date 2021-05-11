function Br_delay = spec_chan_derive_delay(fading,delay,DOA,AOA,d_nor,Nr,L,M,Nt)

Br_delay_tmp = zeros(M,L,Nt);
for jj = 1 : Nt
    for mm = 1 : M
        for l = 1 : L-1
                Br_delay_tmp(mm,l,jj) =fading(mm,jj) * ((sin(l-delay(mm,jj))-((l-delay(mm,jj))*cos(l-delay(mm,jj))))/((l-delay(mm,jj))^2))*exp(-1i*2*pi*d_nor*(Nr-1)*sin(DOA(mm,jj))*sin(AOA(mm,Nr)));
        end
    end
end
Br_delay_tmp1=cell(1,Nt);
for jj = 1 : Nt
Br_delay_tmp1{1,jj}=Br_delay_tmp(:,:,jj);
end
Br_delay=blkdiag(Br_delay_tmp1{:});
end
