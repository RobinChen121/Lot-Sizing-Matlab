function StrictCashConstraintTesting

% ************************************************************************
% Description: single item stochastic lot sizing problem with strict cash constraint;
% ordering policy (s, C, S);
% adopt simulation-based ga to obtain (s, C, S) values ;
%
% author: Zhen Chen
% time: 2018-07-24, 13:16
% ************************************************************************

%% parameters
demands = [15, 15, 15, 15, 15, 15, 15, 15; 
				   21.15, 18.9, 17.7, 16.5, 15.15, 13.95, 12.75, 11.55;
				   6.6, 9.3, 11.1, 12.9, 16.8, 21.6, 24, 26.4; 
				   9, 13, 20, 16, 10, 16, 22, 15;
				   22, 18, 11, 16, 22, 12, 14, 21;
				   41.8, 6.6, 0.4, 21.8, 44.8, 9.6, 2.6, 17;
				   4.08, 12.16, 37.36, 21.44, 39.12, 35.68, 19.84, 22.48;
				   4.7, 8.1, 23.6, 39.4, 16.4, 28.7, 50.8, 39.1;
				   4.4, 11.6, 26.4, 14.4, 14.6, 19.8, 7.4, 18.3;
				   4.9, 18.8, 6.4, 27.9, 45.3, 22.4, 22.3, 51.7];
               
fixCosts = [10, 20]; 
variCosts = [1, 2];
iniCashs = [20, 30];   % can change to fixCost + variCost * d
prices = [4, 8];

global iniCash price fixCost holdCost variCost T demandSample demand;
holdCost = 1;
[N, T] = size(demands);
nvars = 3*T;
intcon = [];
options = optimoptions('ga','Display', 'off');
sampleNum = 1000;

headString = {'K', 'v', 'h', 'I0', 'price', 'B0', 'DemandPatt', 'Time(s)', 'gaIniValue', 'gaSimValue'};
%% circulation
recordResult = zeros(160, 10);
index = 1;
for iDemand = 1 : N    
    demand = demands(iDemand, :);
    demandSample = poissinv(lhsdesign(sampleNum, T), repmat(demand, sampleNum, 1));
    save('demandSample.mat', 'demandSample');
    for iFixCost = 1 : length(fixCosts)
        for iVariCost = 1 : length(variCosts)
            for iPrice = 1 : length(prices)
                for iIniCash = 1 : length(iniCashs)             
                    fixCost = fixCosts(iFixCost);
                    variCost = variCosts(iVariCost);
                    price = prices(iPrice);
                    iniCash = iIniCash * 5 * variCost + fixCost;                    
                    lb = zeros(1, 3* T);
                    meanD = mean(demand);
                    ub = [2*meanD, iniCash*50, 5*meanD];
                    tic
                    [x,fval] = ga(@GasCS, nvars, [], [], [], [], lb, ub, [], intcon, options);
                    time = toc;
                    fprintf('iniGaResult = %.2f\n', -fval);
                    simValue = simulatesCS(x);
                    fprintf('simGAResult = %.2f\n', simValue);
                    recordResult(index, :) = [fixCost, variCost, holdCost, 0, price, iniCash, iDemand, time, -fval, simValue];
                    index = index + 1;
                end
            end
        end
    end   
end
xlswrite('GATestResults.xls', headString);
xlswrite('GATestResults.xls', recordResult, 1, 'A2');



end