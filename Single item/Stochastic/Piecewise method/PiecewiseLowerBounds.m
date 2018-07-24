function PiecewiseLowerBounds  
  
% coded by Chen Zhen, 2016-11-5  
% this programm for getting the lower bounds of a normal distribution  
% output the coordinates of the breakpoints and a figure  
% num is the number of the lower bounds parts  
 
mu = 0;
sigma = 1;
num = 5;

% the complementary loss function  
syms x t;  
CL=matlabFunction(int((x-t)*exp(-0.5*(t-mu)^2/sigma^2)/(sigma*sqrt(2*pi)),t,-inf,x));  
CL1=matlabFunction(int((x-t)*exp(-0.5*t^2)/sqrt(2*pi),t,-inf,x));  
startp=norminv(1e-2,mu,sigma);  
endp=norminv(1-1e-2,mu,sigma);  
  
% fplot(CL,[startp endp],'--b','Linewidth',1.5);  
% axis([startp endp -0.5 CL(endp)]);  
% hold on;  
% fplot(CL(mu),[startp,mu],'r','Linewidth',1.5);  
% fplot(@(x)x-mu+CL(mu),[mu,endp],'r','Linewidth',1.5);  
% legend({'\int_{-\infty}^{x}(x-t)f(t)dt','lowwer bounds'},'FontSize',12,'Location','Northwest');  
% title('example for one lower bounds');  
% hold off;  
  
  
  
%折点的横坐标 输入两个两个相邻的横坐标a, b；下界与上界的折点横坐标不一样
breakp_x=@(x1,x2)(normpdf(x1)-normpdf(x2))/(normcdf(x2)-normcdf(x1));  
  
%下界值，输入两个两个相邻的横坐标a, b， 根据直线公式算出
CL_lb=@(x1,x2)normcdf(x2)*(normpdf(x1)-normpdf(x2))/(normcdf(x2)-normcdf(x1))+normpdf(x2);  
  
e=@(x1,x2)CL1(breakp_x(x1,x2))-CL_lb(x1,x2);  
delta_e=@(x1,x2,x3)e(x1,x2)-e(x2,x3);  
switch num  %论文中的num-1  
    case 3     
        fe=@(x)delta_e(-inf,-x,x)^2+delta_e(-x,x,inf)^2;  
    case 4  
        fe=@(x)delta_e(-inf,-x,0)^2+delta_e(-x,0,x)^2+delta_e(0,x,inf)^2;  
    case 5  
        fe=@(x)delta_e(-inf,-x(2),-x(1))^2+delta_e(-x(2),-x(1),x(1))^2+delta_e(-x(1),x(1),x(2))^2+delta_e(x(1),x(2),inf)^2;  
    case 6  
        fe=@(x)delta_e(-inf,-x(2),-x(1))^2+delta_e(-x(2),-x(1),mu)^2+delta_e(-x(1),mu,x(1))^2+delta_e(mu,x(1),x(2))^2+delta_e(x(1),x(2),inf)^2;  
    case 7  
        fe=@(x)delta_e(-inf,-x(3),-x(2))^2+delta_e(-x(3),-x(2),-x(1))^2+delta_e(-x(2),-x(1),x(1))^2+delta_e(-x(1),x(1),x(2))^2+delta_e(x(1),x(2),x(3))^2+delta_e(x(2),x(3),inf)^2;  
    case 8  
        fe=@(x)delta_e(-inf,-x(3),-x(2))^2+delta_e(-x(3),-x(2),-x(1))^2+delta_e(-x(2),-x(1),mu)^2+delta_e(-x(1),mu,x(1))^2+delta_e(mu,x(1),x(2))^2+delta_e(x(1),x(2),x(3))^2+delta_e(x(2),x(3),inf)^2;  
    case 9  
        fe=@(x)delta_e(-inf,-x(4),-x(3))^2+delta_e(-x(4),-x(3),-x(2))^2+delta_e(-x(3),-x(2),-x(1))^2+delta_e(-x(2),-x(1),x(1))^2+delta_e(-x(1),x(1),x(2))^2+delta_e(x(1),x(2),x(3))^2+delta_e(x(2),x(3),x(4))^2+...  
            delta_e(x(3),x(4),inf)^2;  
    case 10  
        fe=@(x)delta_e(-inf,-x(4),-x(3))^2+delta_e(-x(4),-x(3),-x(2))^2+delta_e(-x(3),-x(2),-x(1))^2+delta_e(-x(2),-x(1),mu)^2+delta_e(-x(1),mu,x(1))^2++delta_e(0,x(1),x(2))^2+delta_e(x(1),x(2),x(3))^2+delta_e(x(2),x(3),x(4))^2+...  
            delta_e(x(3),x(4),inf)^2;  
end  
  
options = optimoptions('fmincon','algorithm','sqp'); % the solution is related with the algorithm selected  
var_num=fix((num-1)/2);  
A=eye(var_num)+diag(-ones(var_num-1,1),1);A(var_num,:)=[];  
b=zeros(var_num-1,1);  
start_point=linspace(0.1,norminv(1-1e-2),var_num);  
[x,fval]=fmincon(fe,start_point,A,b,[],[],zeros(1,var_num),[],[],options);  
if num==3  
    error=sigma*abs(e(x,-x));  
else if num==4  
        error=sigma*abs(e(0,x));  
    else  
        error=sigma*abs(e(x(1),x(2)));  
    end  
end  
if mod(num-1,2)~=0  
    x=[-rot90(x,2),0,x];  
else  
    x=[-rot90(x,2),x];  
end  
  
  
bp=zeros(1,num);  
bp(1)=breakp_x(-inf,x(1));  
bp(end)=breakp_x(x(end),inf);  
for i=2:length(x)  
    bp(i)=breakp_x(x(i-1),x(i));  
end  
  
bp_value=zeros(1,num);  
bp_value(1)=CL_lb(-inf,x(1));  
bp_value(end)=CL_lb(x(end),inf);  
for i=2:length(x)  
    bp_value(i)=CL_lb(x(i-1),x(i));  
end  
bp_value=bp_value*sigma;  

x=sigma*x+mu;  
fprintf('the maximum approximation error for %d lower bounds parts is %.4f\n',num,error);  
fprintf('the lower bounds are:\n');  
fprintf('[%.4f   %.4f]\n',-inf,x(1));  
for i=1:length(x)-1  
    fprintf('[%.4f   %.4f]\n',x(i),x(i+1));  
end  
fprintf('[%.4f   %.4f]\n',x(end),inf);  
fprintf('the x coordinates of breaking points are:\n');  
bp=sigma*bp+mu;  
for i=1:num  
    fprintf('%.4f    ',bp(i));  
end  
fprintf('\n');  
fprintf('the slops of each tangent lines are:\n');  
slops=zeros(1,num+1);  
slops(1)=normcdf(-inf,mu,sigma);  
for i=1:num-1  
    slops(i+1)=normcdf(x(i),mu,sigma);  
end  
slops(end)=normcdf(inf,mu,sigma);  
for i=1:num+1  
    fprintf('%.4f    ',slops(i));  
end  
fprintf('\n');  
fprintf('the vertical interceps are:\n');  
intercepts=zeros(1,num+1);   
intercepts(1)=slops(1)*(-startp)+CL(startp);
for i=2:num  
     intercepts(i)=slops(i)*(-x(i-1))+CL(x(i-1));  
end 
intercepts(end)=slops(end)*(-endp)+CL(endp);
for i=1:num+1  
    fprintf('%.4f    ',intercepts(i));  
end  
fprintf('\n'); 

%-----------------------------------------------  
% draw picture for the results  
figure(2);  
fplot(CL,[startp endp],':b','Linewidth',1.5);  
axis([startp endp -0.5 CL(endp)]);  
hold on;  
fplot(@(t)slops(1)*(t-startp)+CL(startp),[startp,bp(1)],'r','Linewidth',1.5);  
for i=2:num  
    fplot(@(t)slops(i)*(t-x(i-1))+CL(x(i-1)),[bp(i-1),bp(i)],'r','Linewidth',1.5);  
end  
fplot(@(t)slops(end)*(t-endp)+CL(endp),[bp(end),endp],'r','Linewidth',1.5);  
title(['the results of the piecewise approximation I^{+} having ',num2str(num+1),' lowbounds']);  
% for i=1:num  
%     text(bp(i),bp_value(i),['\leftarrow',num2str(i)]);  
% end  
text(0.5*bp(1)+0.5*startp,0.5*bp_value(1)+0.5*CL(startp),['\leftarrow',num2str(1)]); 
for i=1:num-1  
    text(x(i),CL(x(i)),['\leftarrow',num2str(i+1)]);  
end  
text(0.5*bp(end)+0.5*endp,0.5*bp_value(end)+0.5*CL(endp),['\leftarrow',num2str(num+1)]); 
hold off;  
  
end  