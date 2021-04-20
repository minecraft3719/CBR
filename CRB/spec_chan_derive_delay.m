function H = spec_chan_derive_delay(fading,delay, DOA,AOA, Nr, L,M, Nt)
h_derive = zeros(M,L);

H = [];
for r = 1:Nr
    B_derive = [];
    for j = 1:Nt
        for m = 1:M
            for l = 1:L-1   
                  derive_sinc = (sin(l-delay(m,j))-cos(l-delay(m,j)))/(l-delay(m,j))^2;
                h_derive(m,l) =fading(m,j) * derive_sinc*exp(-1i*pi/2*(r-1)*sin(DOA(m,r)+AOA(m,j)));
            end
        end
        B_derive = blkdiag(B_derive,h_derive);
    end
    H = [H B_derive'];
end
