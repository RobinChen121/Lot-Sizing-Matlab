function result = GasCS(input)
% ************************************************************************
% Description: single item stochastic lot sizing problem with strict cash constraint;
% GA function
%
% author: Zhen Chen
% time: 2018-07-25, 13:16
% ************************************************************************

global iniCash price fixCost holdCost variCost T demandSample;
s = input(1 : T);
C = input(T+1 : 2*T);
S = input(2*T + 1 : 3*T);

N = length(demandSample);
results = zeros(N,1);
for i = 1 : N
    B=zeros(1, T); I = zeros(1, T); deltaX = zeros(1, T);
    orderQ = max(0, min((iniCash - fixCost)/variCost, S(1)));
    realD = demandSample(i, :);
    I(1) = orderQ - realD(1);
    I(1)  = max(I(1), 0);
    B(1) = iniCash+price*min(orderQ, realD(1)) - fixCost*deltaX(1) - variCost*orderQ - holdCost*I(1);

    for t = 2 : T
        if B(t - 1) > C(t) && I(t - 1) < s(t)
            orderQ = max(0, min((B(t - 1) - fixCost)/variCost, S(t) - I(t - 1)));
        else
            orderQ = 0;
        end
        if orderQ > 1e-1
            deltaX(t) = 1;
        end
        I(t) = I(t - 1) + orderQ - realD(t);
        B(t) = B(t-1) + price*min(I(t) + orderQ, realD(t)) - fixCost*deltaX(t) - variCost*orderQ - holdCost*I(t);
    end
    results(i)=B(T);
end
result = -mean(results);
end