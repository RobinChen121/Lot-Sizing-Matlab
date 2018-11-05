function Ly = ComputeLy
% ************************************************************************
% Description: compute L(y) for the cash constrainted stochastic lot
% sizing problem.
%Matlab can't get integral for poisson pdf, because it is discrete
%
%
% author: Zhen Chen
% time: 2018-09-13, 17:16
% ************************************************************************

price = 4;
holdCost = 1;	
meand = 4.4;
variCost = 1;

I = 0;
y =  poissinv((price - variCost) / (price+ holdCost), meand);
y = 5;
for i = 0 : y
    I = (y - i ) * poisspdf(i, meand) + I;
end

Ly = price * y - variCost * y - (price + holdCost) * I; 
fprintf('final value is %.2f \n', Ly);

S = poissinv((price - variCost) / (price + holdCost), meand);
fprintf('Order up to level S is %.2f \n', S);

end