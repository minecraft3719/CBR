function H = spec_chan(fading,delay,DOA,AOA,Nr,N,Nt)

for i = 1:Nt
    for j = 1:Nr
        for k = 1:N
            H(j,k,i) = fading(k,i)*sinc(delay(k,i))*exp(-1i*pi*(j-1)*(sin(DOA(k,j))+sin(AOA(k,i))));
        end
    end
end
