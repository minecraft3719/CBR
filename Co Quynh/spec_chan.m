function H = spec_chan(fading,delay,DOA,AOA,Nr,L,M,Nt)
H = zeros(Nr,L,Nt);
for j = 1:Nt
    for r = 1:Nr
	for l = 1:L
	    h = 0;
       	    for m = 1:M
                h = h + (fading(m,j)*sinc(l - delay(m,j))*exp(-1i*pi*(r-1)*(sin(DOA(m,j))+sin(AOA(m,r)))));
 	    end 
	    H(r,l,j) = h;      
	end
    end
end
