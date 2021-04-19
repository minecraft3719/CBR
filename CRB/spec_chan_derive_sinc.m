function H = spec_chan_derive_sinc(fading, delay, DOA,AOA, Nr, N,Nt)
H = zeros(Nr,N,Nt);
for i = 1:Nt
    for j = 1:Nr
        for k = 1:N
            dev = (delay(k,i)*cos(delay(k,i))-sin(delay(k,i)))/delay(k,i)^2;
             H(j,k,i) = fading(k,i)*dev*exp(-1i*pi*(j-1)*(sin(DOA(k,j))+sin(AOA(k,i))));
        end
    end
end