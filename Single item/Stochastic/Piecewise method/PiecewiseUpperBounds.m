function PiecewiseUpperBounds

% coded by Chen Zhen, 2016-11-5
% this programm for getting the upper bounds of a normal distribution
% output the coordinates of the breakpoints and a figure
% num is the number of the upper bounds parts
% 上界部分比下界部分多一个

mu = 10;
sigma = 2;
num = 3;


% the complementary loss function
syms x t;
CL=matlabFunction(int((x-t)*exp(-0.5*(t-mu)^2/sigma^2)/(sigma*sqrt(2*pi)),t,-inf,x)); 

startp=norminv(1e-2,mu,sigma);
endp=norminv(1-1e-2,mu,sigma);

% figure(4);
% fplot(CL,[startp endp],'--b','Linewidth',1.5);
% axis([startp endp -0.5 CL(endp)]);
% hold on;
% fplot(CL(mu),[startp,mu],'r','Linewidth',1.5);
% fplot(@(x)x-mu+CL(mu),[mu,endp],'r','Linewidth',1.5);
% legend({'\int_{-\infty}^{x}(x-t)f(t)dt','lowwer bounds'},'FontSize',12,'Location','Northwest');
% title('example for one lower bounds');
% hold off;
% 
% figure(2);
% fplot(CL,[startp endp],'--b','Linewidth',1.5);
% axis([startp endp -0.5 CL(endp)]);
% hold on;
% fplot(CL(-1.598),[startp,-1.598],'r','Linewidth',1.5);
% plot([-1.598,1.598],[CL(-1.598),CL(1.598)],'r','Linewidth',1.5);
% fplot(@(x)x+CL(1.598)-1.598,[1.598,endp],'r','Linewidth',1.5);
% legend({'\int_{-\infty}^{x}(x-t)f(t)dt','lowwer bounds'},'FontSize',12,'Location','Northwest');
% title('example for one lower bounds');
% hold off;

%误差最大点的横坐标 输入两个两个相邻的横坐标a, b， 论文中的公式
maxe_point=@(x1,x2)norminv((CL(mu+x2)-CL(mu+x1))/(x2-x1),mu,sigma);  

%下界值，其实是纵坐标的下界,输入两个两个相邻的横坐标a, b
CL_ub=@(x1,x2)CL(mu+x1)*(mu+x2-maxe_point(x1,x2))/(x2-x1)+CL(mu+x2)*(maxe_point(x1,x2)-x1-mu)/(x2-x1);

e=@(x1,x2)CL_ub(x1,x2)-CL(maxe_point(x1,x2));
delta_e0=@(x1,x2)CL(mu+x1)-e(x1,x2);
delta_efinal=@(x1,x2)e(x1,x2)+x2-CL(mu+x2);
delta_e=@(x1,x2,x3)e(x1,x2)-e(x2,x3);
switch num  
    case 3   
        fe=@(x)delta_e0(-x,x)^2+delta_efinal(-x,x)^2; 
    case 4
        fe=@(x)delta_e0(-x,0)^2+delta_e(-x,0,x)^2+delta_efinal(0,x)^2;
    case 5
        fe=@(x)delta_e0(-x(2),-x(1))^2+delta_e(-x(2),-x(1),x(1))^2+delta_e(-x(1),x(1),x(2))^2+delta_efinal(x(1),x(2))^2;
    case 6
        fe=@(x)delta_e0(-x(2),-x(1))^2+delta_e(-x(2),-x(1),0)^2+delta_e(-x(1),0,x(1))^2+delta_e(0,x(1),x(2))^2+delta_efinal(x(1),x(2))^2;
    case 7
        fe=@(x)delta_e0(-x(3),-x(2))^2+delta_e(-x(3),-x(2),-x(1))^2+delta_e(-x(2),-x(1),x(1))^2+delta_e(-x(1),x(1),x(2))^2+delta_e(x(1),x(2),x(3))^2+delta_efinal(x(2),x(3))^2;
    case 8
        fe=@(x)delta_e0(-x(3),-x(2))^2+delta_e(-x(3),-x(2),-x(1))^2+delta_e(-x(2),-x(1),0)^2+delta_e(-x(1),0,x(1))^2+delta_e(0,x(1),x(2))^2+delta_e(x(1),x(2),x(3))^2+delta_efinal(x(2),x(3))^2;
    case 9
        fe=@(x)delta_e0(-x(4),-x(3))^2+delta_e(-x(4),-x(3),-x(2))^2+delta_e(-x(3),-x(2),-x(1))^2+delta_e(-x(2),-x(1),x(1))^2+delta_e(-x(1),x(1),x(2))^2+delta_e(x(1),x(2),x(3))^2+...
            +delta_e(x(2),x(3),x(4))^2+delta_efinal(x(3),x(4))^2;
    case 10
        fe=@(x)delta_e0(-x(4),-x(3))^2+delta_e(-x(4),-x(3),-x(2))^2+delta_e(-x(3),-x(2),-x(1))^2+delta_e(-x(2),-x(1),0)^2+delta_e(-x(1),0,x(1))^2+delta_e(0,x(1),x(2))^2+delta_e(x(1),x(2),x(3))^2+...
            +delta_e(x(2),x(3),x(4))^2+delta_efinal(x(3),x(4))^2;
    case 11
        fe=@(x)delta_e0(-x(5),-x(4))^2+delta_e(-x(5),-x(4),-x(3))^2+delta_e(-x(4),-x(3),-x(2))^2+delta_e(-x(3),-x(2),-x(1))^2+delta_e(-x(2),-x(1),x(1))^2+delta_e(-x(1),x(1),x(2))^2+delta_e(x(1),x(2),x(3))^2+...
            +delta_e(x(2),x(3),x(4))^2+delta_e(x(3),x(4),x(5))^2+delta_efinal(x(4),x(5))^2;
end

%options = optimoptions('fmincon','algorithm','sqp','PlotFcns',@optimplotfval); % the solution is related with the algorithm selected
var_num=fix((num-1)/2);
A=eye(var_num)+diag(-ones(var_num-1,1),1);A(var_num,:)=[];
b=zeros(var_num-1,1);
start_point=linspace(sigma/10,norminv(1-1e-2,0,sigma),var_num);
[x,fval]=fmincon(fe,start_point,A,b,[],[],zeros(1,var_num));%,[],[],options

if mod(num-1,2)~=0
    x=[rot90(mu-x,2),mu,mu+x];
else
    x=[rot90(mu-x,2),mu+x];
end
error=CL(x(1));
fprintf('the maximum approximation error for %d lower bounds parts is %.4f\n',num,error);
fprintf('the upper bounds are:\n');
fprintf('[%.4f   %.4f]\n',-inf,x(1));
for i=1:length(x)-1
    fprintf('[%.4f   %.4f]\n',x(i),x(i+1));
end
fprintf('[%.4f   %.4f]\n',x(end),inf);
fprintf('the slops of each tangent lines are:\n');
slops=zeros(num,1);
slops(1)=0;slops(end)=1;
for i=2:num-1
    slops(i)=(CL(x(i))-CL(x(i-1)))/(x(i)-x(i-1));
end
for i=1:num
    fprintf('%.4f    ',slops(i));
end
fprintf('\n');
bp=zeros(num,1);bp_value=zeros(num,1);
bp(1)=0.5*(x(1)+startp);bp_value(1)=CL(x(1));
bp(num)=0.5*(endp+x(end));bp_value(num)=bp(num)+CL(x(end))-x(end);
for i=2:num-1
    bp(i)=maxe_point(x(i-1)-mu,x(i)-mu);bp_value(i)=bp(i)*slops(i)+(x(i)*CL(x(i-1))-x(i-1)*CL(x(i)))/(x(i)-x(i-1));
end
fprintf('the vertical interceps are:\n');  
intercepts=zeros(1,num);   
intercepts(1)=CL(x(1));
for i=2:num-1  
     intercepts(i)=slops(i)*(-x(i-1))+CL(x(i-1));  
end 
intercepts(end)=slops(end)*(-x(end))+CL(x(end));
for i=1:num  
    fprintf('%.4f    ',intercepts(i));  
end  
fprintf('\n'); 

%-----------------------------------------------
% draw picture for the results
figure(3);
fplot(CL,[startp endp],':b','Linewidth',1.5);
axis([startp endp -0.5 CL(endp)]);
hold on;
plot([startp,x(1)],[CL(x(1)),CL(x(1))],'r','Linewidth',1.5);
for i=2:num-1
    fplot(@(t)slops(i)*(t-x(i-1))+CL(x(i-1)),[x(i-1),x(i)],'r','Linewidth',1.5);
end
fplot(@(t)slops(end)*(t-x(end))+CL(x(end)),[x(end),endp],'r','Linewidth',1.5);
%title(['the results of the piecewise approximation I^{+} having ',num2str(num),' upbounds']);
for i=1:num
    text(bp(i),bp_value(i),['\leftarrow',num2str(i)]);
end
legend('期望库存量曲线', '拟合线段', 'location', 'Northwest');
hold off;

end