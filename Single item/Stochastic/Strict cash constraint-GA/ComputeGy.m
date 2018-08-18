function ComputeGy(y, meand)

% ************************************************************************
% Description: compute G(y) for the cash constrainted stochastic lot
% sizing problem.
%Matlab can't get integral for poisson pdf, because it is discrete
%
%
% author: Zhen Chen
% time: 2018-07-31, 13:16
% ************************************************************************



variCost = 1;
price = 8;
holdCost = 1;	
penalCost = 0;

holdCosts = 0;
penalCosts = 0;
for i = 0 : y
    holdCosts = holdCost * (y - i ) * poisspdf(i, meand) + holdCosts;
end
for i = y : 500
    penalCosts = penalCost * (i - y) * poisspdf(i, meand) + penalCosts;
end

Gy = price * y - variCost * y - holdCosts - penalCosts - price/holdCost * holdCosts;

fprintf('final value is %.2f \n', Gy);


end