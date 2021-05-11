function H = spec_chan_derive_AOA(fading,delay, DOA,AOA, Nr, L,M, Nt)
h_derive = zeros(M,L);

H = [];
for r = 1:Nr
    B_derive = [];
    for j = 1:Nt
        for m = 1:M
            for l = 1:L-1   
                    h_derive(m,l) =fading(m,j)* sinc(l - delay(m,j))*(-1i*pi/2*sin(DOA(m,j))*cos(AOA(m,r)))*exp(-1i*pi/2*(r-1)*sin(DOA(m,j))*sin(AOA(m,r)));
            end
        end
        B_derive = blkdiag(B_derive,h_derive);
    end
    H = [H; B_derive'];
end
