function logL_CIR = CIRlog(data,dt,alpha,mi,sigma)
Nobs = length(data);
DataL = data(1:end-1); %empirical data
DataF = data(2:end);

c = (2*alpha)/((sigma^2)*(1-exp(-alpha*dt)));
q = ((2*alpha*mi)/(sigma^2))-1;
u = c*exp(-alpha*dt)*DataL;
v = c*DataF; 
z = 2*sqrt(u.*v);
bf = besseli(q,z,1);
logL_CIR = -(Nobs - 1)*log(c) + sum(u + v - 0.5*q*log(v./u) - log(bf) - z);
% the transition density function 
%c*exp(-u-v).*(v./u).^(q/2).*bf.*exp(z)
end
