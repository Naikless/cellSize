function [ lambda ] = cellSizeCorr( epsilon, delta, sigma, ucj )

A0 = 30.465860763763;
a1 = 89.55438805808153;
a2 = -130.792822369483;
a3 = 42.02450507117405;
b1 = -0.02929128383850;
b2 = 1.026325073064710e-5;
b3 = -1.031921244571857e-9;

zeta = epsilon .* delta .* sigma ./ ucj;

A = A0 + a3./zeta.^3 + a2./zeta.^2 + a1./zeta + b1.*zeta + b2.*zeta.^2 + b3.*zeta.^3; 
lambda = A .* delta;

end

