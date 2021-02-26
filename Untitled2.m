clc;
clear all;

% linear 
theta = [0:360] * pi/180; 
phi = [0:360] * pi/180;

d = 1;
f = 32e+8;
lambda = (3e+8)/f;

angle = (pi*d/lambda)*sin(theta).*sin(phi).' ;

N = 10;
n = [1:N];

for i = 1: N
    u(:,:,i) = (2*i-1) .*angle;
end
AF = cos(u);
AF_u = sum(AF,3);

[PHI,THETA] = meshgrid(phi,theta);
[X,Y,Z] = sph2cart(THETA,PHI,abs(AF_u));

figure(1)
surf(X,Y,Z);
legend('Uniform');
