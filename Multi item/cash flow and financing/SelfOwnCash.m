function SelfOwnCash
% multi item lot sizing problem with onle self owned cash
% consider 2 items, 8 periods first
% initial cash: B0, set up cost: s, variable production cost: v
% inventory flow: I, cash flow: B
% unit holding cost: h, selling price: p
% payment delay length: L
% discount rate: r
% overhead costs: K
% minimum cash balance: Bmin
% whether or not set up x; how much to produce: y; lost sales quantity: w

%% parameter values
B0 = 70000;       M = 1000000;
p1 = [77,79,74,77,72,70,68,67,64,61,59,57,53,55,55,56,60,58,59,58,64,59,64,64,66,64,64,71,67,68,72,71,70,75,76,78,75,74,73,71,73,73,70,72,65,63,58,59,53,55,57,61];    
p2 = [56,56,60,57,56,60,58,60,63,63,64,61,62,60,58,52,52,48,50,48,45,50,46,50,53,54,59,59,62,60,64,66,71,73,72,74,78,81,80,86,82,82,80,75,74,73,66,67,65,64,58,58];
d1= [275,284,303,304,318,324,344,357,364,368,383,392,393,412,421,400,395,383,381,374,355,358,344,324,315,307,297,297,289,263,261,279,288,284,292,301,327,334,346,347,364,369,371,392,384,379,365,359,330,333,324,314];    
%d2 = [243,260,261,280,271,294,296,305,316,335,335,354,350,377,382,384,375,353,352,330,326,322,300,299,292,279,264,262,286,291,297,318,311,321,343,344,352,363,373,384,392,395,389,378,366,342,350,333,317,318,299,289];
T=length(p1);
d2=zeros(1,T);
s1 = 5000*ones(1,T);  s2 = 5000*ones(1,T);        
v1 = 31*ones(1,T);    v2 = 30*ones(1,T);  
h1= 10*ones(1,T);     h2= 10*ones(1,T);  
I10=0; I20=0;

L = 3;  alpha = 0.25;
K = 6000*ones(1,T);     Bmin = 0*ones(1,T);


%% decision variables
x1 = binvar(T,1);   x2 = binvar(T,1);
y1 = sdpvar(T ,1);      y2 = sdpvar(T,1);
w1 = sdpvar(T ,1);      w2 = sdpvar(T,1);


%% objective function
for t = 1 : T
    if t == 1
        I1 = I10 + y1(t) - d1(t) + w1(t);       I2 = I20 + y2(t) - d2(t) + w2(t);
        if t <= L
            B = B0 - s1(t)*x1(t) - s2(t)*x2(t) - v1(t)*y1(t) - v2(t)*y2(t) - h1(t)*I1(t) - h2(t)*I2(t)-K(t);
        else
            B = B0 + p1(t-L)*(d1(t-L) - w1(t-L)) + p2(t-L)*(d2(t-L) - w2(t-L)) - s1(t)*x1(t) - s2(t)*x2(t) - v1(t)*y1(t) - v2(t)*y2(t) - h1(t)*I1(t) - h2(t)*I2(t)-K(t);
        end
    end
    
    if t > 1
        I1 = [I1; I1(t-1) + y1(t) - d1(t) + w1(t)];       
        I2 = [I2; I2(t-1) + y2(t) - d2(t) + w2(t)];
        if t <= L
            B = [B; B(t-1)- s1(t)*x1(t) - s2(t)*x2(t) - v1(t)*y1(t) - v2(t)*y2(t) - h1(t)*I1(t) - h2(t)*I2(t)-K(t)];
        else
            B = [B; B(t-1) + p1(t-L)*(d1(t-L) - w1(t-L)) + p2(t-L)*(d2(t-L) - w2(t-L)) - s1(t)*x1(t) - s2(t)*x2(t) - v1(t)*y1(t) - v2(t)*y2(t) - h1(t)*I1(t) - h2(t)*I2(t)-K(t)];
        end
    end
    
    if t == T
        for j = T-L+1 : T
            tempB = p1(j)*(d1(j) - w1(j)) / (1 + alpha)^(L+ j - T) + p2(j)*(d2(j) - w2(j)) / (1 + alpha)^(L+ j - T) ;
            B(t) = B(t) + tempB;
        end
    end
end


%% constraints
xyC = [y1(1:T) <= M*x1(1:T);     y2(1:T) <= M*x2(1:T)];
wdC = [w1(1:T) <= d1(1:T)' ;       w2(1:T) <= d2(1:T)' ];
BC = [B0 - s1(1)*x1(1)-s2(1)*x2(1)-v1(1)*y1(1)-v2(1)*y2(1) >= Bmin(1)];
for t = 2:T
    BC = [BC; B(t-1)-s1(t)*x1(t)-s2(t)*x2(t)-v1(t)*y1(t)-v2(t)*y2(t) >= Bmin(t)];
end
IC = [I1 >= 0; I2 >= 0];
wC = [w1 >= 0; w2 >= 0];
yC = [y1 >= 0; y2 >= 0];
constraints = [xyC; wdC; BC; IC; wC; yC];

%% solve
tic
options = sdpsettings('solver','cplex');
sol = optimize(constraints, -B(T), options);
toc

%% output results
if sol.problem == 0
    fprintf('\n\n');
    fprintf('payment delay length: L = %d\n', L);
    fprintf('inventory level at each periods: I1 = [');
    for i = 1:T
        fprintf('  %.2f ',value(I1(i)));
    end
    fprintf(' ]\n');
    fprintf('inventory level at each periods: I2 = [');
    for i = 1:T
        fprintf('  %.2f ',value(I2(i)));
    end
    fprintf(' ]\n');
    
    fprintf('cash balance at each periods: B = [');
    for i = 1:T
        fprintf('  %.2f ',value(B(i)));
    end
    fprintf(' ]\n');
    
        fprintf('whether set up production: x1 = [');
    for i = 1:T
        fprintf('  %.2f ',value(x1(i)));
    end
    fprintf(' ]\n');
    fprintf('whether set up production: x2 = [');
    for i = 1:T
        fprintf('  %.2f ',value(x2(i)));
    end
    fprintf(' ]\n');
    
    fprintf('how much to produce: y1 = [');
    for i = 1:T
        fprintf('  %.2f ',value(y1(i)));
    end
    fprintf(' ]\n');
    fprintf('how much to produce: y2 = [');
    for i = 1:T
        fprintf('  %.2f ',value(y2(i)));
    end
    fprintf(' ]\n');
    
    fprintf('lost sales: w1 = [');
    for i = 1:T
        fprintf('  %.2f ',value(w1(i)));
    end
    fprintf(' ]\n');
    fprintf('lost sales: w2 = [');
    for i = 1:T
        fprintf('  %.2f ',value(w2(i)));
    end
    fprintf(' ]\n');
    
    fprintf('\n');
    fprintf('optimal final cash increment is %4.2f\n', value(B(T))-B0);
    fprintf('\n\n');
else
    fprintf('no feasible result\n');
end

end