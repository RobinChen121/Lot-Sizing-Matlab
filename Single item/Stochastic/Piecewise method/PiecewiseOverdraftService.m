function PiecewiseOverdraftService

% parameters:
a=50;h=1;pi=0;
meand=[5,15,26,44];%,24,15,22,10];
sigma=[1,1,1,1];%[1.5,4.5,7.8,13.2];%,7.2,4.5,6.6,3];

I0=0;v=1;B0=50;rate=0.1;price=10;alpha=0.95;
M=5000;T=length(meand);

% piecewise values:
nbpartitions = 5;
prob = [0.1324110437406592, 0.23491250409192982, 0.26535290433482195, 0.23491250409192987, 0.13241104374065915];
means = [-1.6180463502161044, -0.6914240068499904, 0, 0.6914240068499903, 1.6180463502161053];
error = 0.022270929512393414;

% standard variance values
sigma_sum=zeros(T,T);
for t=1:T
    for j=1:t
        sigma_sum(j,t)=sqrt(sum(sigma(j:t).^2));
    end
end

% decision values
delta=binvar(T,1);
I=sdpvar(T,1);
Iplus=sdpvar(T,1);
Iminus=sdpvar(T,1);
P=binvar(T,T);
P=triu(P);

% objective function
%B=B0+price*(sum(meand)-Iminus(T))-(a*sum(delta)+h*sum(Iplus)+pi*sum(Iminus)+v*I(T)-v*I0+v*sum(meand));
B=B0+price*(meand(1)-Iminus(1))-(a*delta(1)+h*Iplus(1)+pi*Iminus(1)+v*(I(1)+meand(1)-I0));
for t=2:T
    B=[B;B(t-1)+rate*min(B(t-1),0)+price*(meand(t)-Iminus(t)+Iminus(t-1))-(a*delta(t)+h*Iplus(t)+pi*Iminus(t)+v*(I(t)+meand(t)-I(t-1)))];
end
B=B+rate*min(B(t),0);
B=-B;

% constraints
deltaC=[I(1)-I0+meand(1)<=M*delta(1);I(2:T)+meand(2:T)'-I(1:T-1)<=M*delta(2:T)];
QC0=[I(1)+meand(1)>=I0;I(2:T)+meand(2:T)'>=I(1:T-1)];
PC1=(sum(P)==1);
PC2=[];
for t=1:T
    for j=1:t
        PC2=[PC2;P(j,t)>=delta(j)-sum(delta(j+1:t))];
    end
end
PieceC1=[];IC=[];
for t=1:T
    for i=1:nbpartitions
        %PieceC1=[PieceC1;Iplus(t)>=sum(prob(1:i))*I(t)-prob(1:i)*means(1:i)'*P(1:t,t)'*sigma_sum(1:t,t)];
        PieceC1=[PieceC1;Iplus(t)>=sum(prob(1:i))*I(t)-prob(1:i)*means(1:i)'*P(1:t,t)'*sigma_sum(1:t,t)+error*P(1:t,t)'*sigma_sum(1:t,t)];
    end
    %PieceC1=[PieceC1;Iplus(t)>=0];
    PieceC1=[PieceC1;Iplus(t)>=error*P(1:t,t)'*sigma_sum(1:t,t)];
    IC=[IC;I(t)>=norminv(alpha)*sigma_sum(1:t,t)'*P(1:t,t)];
end

PieceC2=[];
for t=1:T
    for k=1:nbpartitions
        %PieceC2=[PieceC2;Iminus(t)>=-I(t)+sum(prob(1:k))*I(t)-prob(1:k)*means(1:k)'*P(1:t,t)'*sigma_sum(1:t,t)];
      PieceC2=[PieceC2;Iminus(t)>=-I(t)+sum(prob(1:k))*I(t)-prob(1:k)*means(1:k)'*P(1:t,t)'*sigma_sum(1:t,t)+error*P(1:t,t)'*sigma_sum(1:t,t)];
    end
    %PieceC2=[PieceC2;Iminus(t)>=-I(t)];
    PieceC2=[PieceC2;Iminus(t)>=-I(t)+error*P(1:t,t)'*sigma_sum(1:t,t)];
end

constraints=[deltaC;QC0;PC1;PC2;PieceC1;PieceC2,IC];
options=sdpsettings('solver','cplex');
optimize(constraints,B(end),options);

fprintf('各阶段库存水平 I= \n');
for i=1:T
    fprintf('  %.2f',value(I(i)));
end
fprintf('\n');
fprintf('各阶段库存量 I+ = \n');
for i=1:T
    fprintf('  %.4f',value(Iplus(i)));
end
fprintf('\n');
fprintf('各阶段缺货量 I- =\n');
for i=1:T
    fprintf('  %.4f',value(Iminus(i)));
end
fprintf('\n');
fprintf('各阶段是否订货\n');
for i=1:T
    fprintf('  %.2f',value(delta(i)));
end
fprintf('\n');
fprintf('各阶段补货上限 R=\n');
for i=1:T
    fprintf('  %.4f',value(I(i)+meand(i)));
end
fprintf('\n');
fprintf('各阶段订货量 Q=\n');
fprintf('  %.2f',value(I(1)+meand(1)-I0));
for i=2:T
    fprintf('  %.2f',value(I(i)+meand(i)-I(i-1)));
end
fprintf('\n');
fprintf('各阶段是否生产：\n');
for i=1:T
    fprintf('  %.2f',value(delta(i)));
end
fprintf('\n');
fprintf('初始资金：%.2f    ',B0);
fprintf('销售收益：%.2f    ',value(price*(sum(meand)-Iminus(T))));
fprintf('生产总成本：');
fprintf('%.2f    ',value(a*sum(delta)+h*sum(Iplus)+pi*sum(Iminus)+v*I(T)-v*I0+v*sum(meand)));
fprintf('\n');
fprintf('最终资金：');
interest=0;
for t=1:length(B)
    fprintf('%.4f  ',value(-B(t)));
    interest=-rate*min(value(-B(t)),0)+interest;
end
fprintf('\n');
fprintf('支付的利息 %.2f',value(interest));
fprintf('\n');
end