function PiecewiseCapacityConstraint

mu = 15;
sigma = 2;
Qup = 25;

% the expected inventory
syms x t;
EI=matlabFunction(int((x-t)*exp(-0.5*(t-mu)^2/sigma^2)/(sigma*sqrt(2*pi)),t,-inf,x)); 
x = 0:0.1:Qup;
plot(x, EI(x));

% given a maximum approximation error
error = 0.1;

% the tangent line through the capacity point
Derivative = @(q) normcdf(q, mu, sigma);



fprintf('derivative at capacity point: %4.2f\n',  Derivative(Qup));


end