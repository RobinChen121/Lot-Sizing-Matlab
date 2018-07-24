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

%%

demand = [15, 15, 15, 15, 15, 15, 15, 15];
iniCash = 20;
price = 4;
fixCost = 10;
holdCost = 1;
variCost = 2;

T = length(demand);

%% generate latin hypercube samples
sampleNum = 1000;
demandSample = poissinv(lhsdesign(sampleNum, T), repmat(demand, sampleNum, 1));


%% invoke ga





end