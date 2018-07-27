function StrictCashConstraint

% ************************************************************************
% Description: single item stochastic lot sizing problem with strict cash constraint;
% ordering policy (s, C, S);
% adopt simulation-based ga to obtain (s, C, S) values ;
% invoke GasCs.m and simulatesCS.m
%
% author: Zhen Chen
% time: 2018-07-24, 13:16
% ************************************************************************


%% parameters
global iniCash price fixCost holdCost variCost T demandSample demand;
demand = [41.8,6.6,0.4,21.8,44.8,9.6,2.6,17];
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
fprintf('simGAResult = %.2f\n', simulatesCS);

end


