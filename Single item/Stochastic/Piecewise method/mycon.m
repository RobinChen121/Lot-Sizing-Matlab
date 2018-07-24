function [c,ceq]=mycon(x)
    syms Q y;
    comp_loss=matlabFunction(int((Q-y)*(1/sqrt(2*3.14159*2))*exp(-(y-1)^2/2),-inf,Q));
    comp_loss_lb=@(x)normcdf(x(2))*(normpdf(x(1))-normpdf(x(2)))/(normcdf(x(2))-normcdf(x(1)))+...
                  normpdf(x(2));
    breakpoint_x=@(x)(comp_loss(x(2))-comp_loss(x(1))+x(1)*normcdf(x(1))-...
                  x(2)*normcdf(x(2)))/(normcdf(x(1))-normcdf(x(2)));
    ceq=[];
    
    %4个区间，一个未知数
    %c=[comp_loss(breakpoint_x([-10,-x]))-comp_loss_lb([-10,-x])-comp_loss(breakpoint_x([-x,0]))+comp_loss_lb([-x,0])-1e-4;
    %comp_loss(breakpoint_x([-x,0]))-comp_loss_lb([-x,0])-comp_loss(breakpoint_x([0,x]))+comp_loss_lb([0,x])-1e-4;    
    %comp_loss(breakpoint_x([0,x]))-comp_loss_lb([0,x])-comp_loss(breakpoint_x([x,10]))+comp_loss_lb([x,10])-1e-4];
    
    %5个区间，两个未知数
    c=[comp_loss(breakpoint_x([-10,-x(1)]))-comp_loss_lb([-10,-x(1)])-comp_loss(breakpoint_x([-x(1),-x(2)]))+comp_loss_lb([-x(1),-x(2)])-1e-4;
    comp_loss(breakpoint_x([-x(1),-x(2)]))-comp_loss_lb([-x(1),-x(2)])-comp_loss(breakpoint_x([-x(2),x(2)]))+comp_loss_lb([-x(2),x(2)])-1e-4;    
    comp_loss(breakpoint_x([-x(2),x(2)]))-comp_loss_lb([-x(2),x(2)])-comp_loss(breakpoint_x([x(2),x(1)]))+comp_loss_lb([x(2),x(1)])-1e-4;
    comp_loss(breakpoint_x([x(2),x(1)]))-comp_loss_lb([x(2),x(1)])-comp_loss(breakpoint_x([x(1),10]))+comp_loss_lb([x(1),10])-1e-4];
end
