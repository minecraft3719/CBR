function AF = circular(a, N, theta, phi, In)
%a: radius (m), %N: number of antenna, %theta & phi: angle (degree),
%In: magnitude coefficient = 1 for UCA
k = 8.98*10e9;
alpha = zeros(N,1);
phin = zeros(N,1);
AF1 = zeros(N,1);

for n = 1:N
    phin(n,1) = 2*pi*(n-1)/N;
    alpha(n,1) = -k*a*sin(theta)*cos(phi - phin(n,1));
    AF1(n,1)= k*a*(sin(theta)*cos(phi-phin(n,1)))+alpha(n,1);
end
AF1 = In*exp(AF1*1i);
plot(AF1);
AF = sum(AF1);
end




