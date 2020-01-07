 % The Euler approximation method for CIR model
 function [t, yEuler] = Euler(alpha,mi,sigma,T,x0,Nsteps,nsim)
timestep = T/(Nsteps-1); % arbitrary increment _ t
 y = zeros(Nsteps,nsim); % create a Nsteps-by-nsim array of zeros
 y(1, :) = x0*ones(1,nsim); % nsim is the number of simulations
 for i = 1:(Nsteps-1);
 y(i+1,:)=y(i,:)+alpha*(mi-y(i,:))*timestep+sigma*sqrt(abs(y(i,:))).*randn(1,nsim)*sqrt(timestep);
 end
 t = linspace (0,T,Nsteps); % generating Nsteps points and the space between the points is _ t
 yEuler = y;
end

