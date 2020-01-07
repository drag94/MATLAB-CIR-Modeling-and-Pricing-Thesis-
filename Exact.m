 % The exact algorithm for CIR model
function [t, xExact] = Exact(alpha,mi,sigma,T,x0,Nsteps,nsim)
 timestep = T/(Nsteps-1);
 n = (4*alpha*exp(-alpha*timestep))/((sigma^2)*(1-exp(-alpha*timestep)));
 % x(T) or n(t+timestep)
 df = (4*alpha*mi)/(sigma^2); % degrees of freedom
 x = zeros(Nsteps,nsim);
 x(1, :) = x0*ones(1,nsim);
 for t = 1:Nsteps-1;
 ncp = x(t,:)*n; % non-centrality parameter x(t)x(T) or x(t)n(t+timestep)
 % ncsrd = ncx2rnd(d, ncp, 1, 1); % non-central chi-square random var
 Po = poissrnd(ncp/2);
 v = df+2*Po;
 % chi2prime = (randn+sqrt(ncp)ˆ2)+ncsrd when d ...greater than 1
 % x(t+1) = chi2prime*
 chi2v = chi2rnd(v);
 x(t+1,:)= chi2v*(exp(-alpha*timestep)/n);
 end
 t = linspace (0,T,Nsteps);
 xExact = x;
end


