function ComputeIy(y)
% ************************************************************************
% Description: compute I(y) for the cash constrainted stochastic lot
% sizing problem.
%Matlab can't get integral for poisson pdf, because it is discrete
%
%
% author: Zhen Chen
% time: 2018-09-11, 17:16
% ************************************************************************

price = 5;
holdCost = 1;	
meand1 = 13;
meand2 = meand1 + 20;
meand3 = meand2 + 16;
I1 = 0;
I2 = 0;
I3 = 0;
for i = 0 : y
    I1 = (y - i ) * poisspdf(i, meand1) + I1;
    I2 = (y - i ) * poisspdf(i, meand2) + I2;
    I3 = (y - i ) * poisspdf(i, meand3) + I3;
end

Ly = price * y -(price + holdCost) * I3 - holdCost * I2 - holdCost * I1;

fprintf('final value is %.2f \n', Ly);
end