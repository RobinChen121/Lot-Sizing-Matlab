function piecewise_lower_bounds  

% this programm for getting the lower bounds of a normal distribution��
% ��׼��̬�ֲ���ֱ���㣬���任������ȣ��������
% output the coordinates of the breakpoints and a figure

% the complementary loss function
syms x t;
%CL=matlabFunction(int((x-t)*(1/sqrt(2*pi))*exp(-t^2/2),t,-inf,x));
CL=matlabFunction(2*int((x-t)*(1/sqrt(2*pi))*exp(-t^2/2),t,-inf,x)+10*int((t-x)*(1/sqrt(2*pi))*exp(-t^2/2),t,x,inf));
fplot(CL,[-20 20],'b','Linewidth',1.5);
axis([-2 2 -0.5 2]);
hold on;
fplot(0,[-2,0],'--r','Linewidth',1.5);
fplot(@(x)x,[0,2],'--r','Linewidth',1.5);
legend({'\int_{-\infty}^{x}(x-t)f(t)dt','lowwer bounds'},'FontSize',12,'Location','Northwest');

hold off;
%�½磬�۵�������꣬���������������ڵĺ�����
comp_loss_lb=@(x)normcdf(x(2))*(normpdf(x(1))-normpdf(x(2)))/(normcdf(x(2))-normcdf(x(1)))+...
                  normpdf(x(2));

%��ʵֵ������һ��������
comp_loss=matlabFunction(int((t-x)*(1/sqrt(2*pi))*exp(-x^2/2),-inf,t));
%comp_loss1=matlabFunction(2*(-Q+2*Q*int((1/sqrt(2*3.14159))*exp(-x^2/2),-inf,Q)+2*(1/sqrt(2*3.14159))*exp(-Q^2/2)));
%comp_loss1=matlabFunction(-Q+2*Q*(int((1/sqrt(2*3.14159))*exp(-x^2/2),-inf,Q))+2*(1/sqrt(2*3.14159))*exp(-Q^2/2));
%comp_loss1=matlabFunction(-Q+2*Q*(int((1/sqrt(2*3.14159))*exp(-x^2/2),-inf,Q))+2*(1/sqrt(2*3.14159))*exp(-Q^2/2));
%comp_loss2=matlabFunction(-Q+1+2*int((Q-x)*(1/sqrt(2*3.14159*2))*exp(-(x-1)^2/2),-inf,Q));
%comp_loss2=matlabFunction(int((Q-x)*(1/sqrt(2*3.14159))*exp(-(x-1)^2/2),-inf,Q)-int((Q-x)*(1/sqrt(2*3.14159))*exp(-(x-1)^2/2),Q,+inf));
%ezplot(comp_loss1);
%hold on;
%ezplot(comp_loss2);

%legend('E\{C\}','\sigma f(y)');

%�۵�ĺ����꣬ ����ֱ�߹�ʽ�����������
breakpoint_x=@(x)(comp_loss(x(2))-comp_loss(x(1))+x(1)*normcdf(x(1))-...
                  x(2)*normcdf(x(2)))/(normcdf(x(1))-normcdf(x(2)));


%���ܾ���׾���������x������������
%3���䣬һ��δ֪��������ͨ�� 
% 4�����䣬һ��δ֪��������ͨ��
%error=@(x)(comp_loss(breakpoint_x([-10,-x]))-comp_loss_lb([-10,-x]))^2 ...
          %+(comp_loss(breakpoint_x([-x,0]))-comp_loss_lb([-x,0]))^2 ...
          %+(comp_loss(breakpoint_x([0,x]))-comp_loss_lb([0,x]))^2 ...
          %+(comp_loss(breakpoint_x([x,10]))-comp_loss_lb([x,10]))^2;

%5�����䣬����δ֪��
error=@(x)(comp_loss(breakpoint_x([-10,-x(1)]))-comp_loss_lb([-10,-x(1)]))^2 ...
          +(comp_loss(breakpoint_x([-x(1),-x(2)]))-comp_loss_lb([-x(1),-x(2)]))^2 ...
          +(comp_loss(breakpoint_x([-x(2),x(2)]))-comp_loss_lb([-x(2),x(2)]))^2 ...
          +(comp_loss(breakpoint_x([x(2),x(1)]))-comp_loss_lb([x(2),x(1)]))^2 ...
          +(comp_loss(breakpoint_x([x(1),10]))-comp_loss_lb([x(1),10]))^2;
      
      
%options=optimoptions(@fmincon,'Algorithm','sqp');
[xx,fval]=fmincon(error,[0.5,0.8],[],[],[],[],[-10,-10],[10,10]); %��Ҫ�޸�,��Լ��������Ҫ��fmincon Լ������,@mycon
error_value=comp_loss(breakpoint_x([-10,-xx(1)]))-comp_loss_lb([-10,-xx(1)]);
breakpoint_x([-10,-xx(1)])
end