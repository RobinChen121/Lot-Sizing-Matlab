function StrictCashConstraint

% ************************************************************************
% Description: single item stochastic lot sizing problem with strict cash constraint;
% ordering policy (s, C, S);
% adopt simulation-based ga to obtain (s, C, S) values ;
%
% author: Zhen Chen
% email: 15011074486@163.com
% time: 2018-07-24, 13:16
% ************************************************************************


%% parameters
demand = [41.8,6.6,0.4,21.8,44.8,9.6,2.6,17];
global iniCash price fixCost holdCost variCost T demandSample;
iniCash = 20;
price = 4;
fixCost = 10;
holdCost = 1;
variCost = 1;
T = length(demand);


%% generate latin hypercube samples
sampleNum = 1000;
demandSample = poissinv(lhsdesign(sampleNum, T), repmat(demand, sampleNum, 1));
save('demandSample.mat', 'demandSample');


%% invoke ga
nvars = 3*T;
intcon = [];
lb = zeros(1, 3* T);
meanD = mean(demand);
ub = [2*meanD, iniCash*50, 5*meanD];
options = optimoptions('ga','Display', 'off');
tic
[x,fval] = ga(@GasCS, nvars, [], [], [], [], lb, ub, [], intcon, options);
toc
fprintf('iniGaResult = %.2f\n', -fval);


%% simulation (s, C, S)
sampleNum = 10000;
demandSamples = poissinv(lhsdesign(sampleNum, T), repmat(demand, sampleNum, 1));
s = x(1 : T);
C = x(T+1 : 2*T);
S = x(2*T + 1 : 3*T);

N = length(demandSamples);
results = zeros(N,1);
for i = 1 : N
    B=zeros(1, T); I = zeros(1, T); deltaX = zeros(1, T);
    orderQ = max(0, min((iniCash - fixCost)/variCost, S(1)));
    realD = demandSamples(i, :);
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
simValue = mean(results);
fprintf('simGAResult = %.2f\n', simValue);
end


%% fuction for (s, C, S) policy
function result = GasCS(input)
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
