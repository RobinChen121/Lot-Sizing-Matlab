function Ly = ComputeLy(y)
% ************************************************************************
% Description: compute L(y) for the cash constrainted stochastic lot
% sizing problem.
%Matlab can't get integral for poisson pdf, because it is discrete
%
%
% author: Zhen Chen
% time: 2018-09-13, 17:16
% ************************************************************************

price = 8;
holdCost = 3;	
meand = 22.2;
variCost = 1;

I = 0;
%y =  poissinv((price - variCost) / (price+ holdCost), meand);
for i = 0 : y
    I = (y - i ) * poisspdf(i, meand) + I;
end

Ly = price * y - variCost * y - (price + holdCost) * I; 
fprintf('final value is %.2f \n', Ly);

% Ly = zeros(1, 50);
% for y = 1 : 40
% for i = 0 : y
%     I = (y - i ) * poisspdf(i, meand) + I;
% end
% 
% if y <= 10
%     Ly1(y) = price * y -(price + holdCost) * I + 2*(y-10)*0.2;
% else
%     Ly2(y) = price * y -(price + holdCost) * I + 2*(y-10)*0.1;
% end
% end
% 
% plot( Ly1,'blue');
% figure;
% plot(Ly2,'red');

end