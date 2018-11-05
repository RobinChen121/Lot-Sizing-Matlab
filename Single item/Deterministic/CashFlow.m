function CashFlow

% ************************************************************************
% Description: compute cash flow for single item deterministic lot sizing
% problem
% d: demands
% p: price
% s: fixed ordering cost
% c: unit vari ordering cost
% h: unit holding cost
%
% author: Zhen Chen
% time: 2018-09-27, 17:16
% ************************************************************************


d = [6, 18.8, 6.4];  
N = length(d);
p = 4 * ones(1, N);
s = 20 * ones(1, N);
c = 1 * ones(1, N);
h = 1 * ones(1, N);
pai = 0 * ones(1, N);
B0 = 30;
BL = 0; L = 0; r = 0;

BB=zeros(N,N);
B=zeros(1,N); %各阶段期末资金
x=zeros(1,N);
y=zeros(1,N);
w=zeros(1,N);
mark=zeros(1,N);
initial_B=zeros(1,N+1);

for i=1:N
    for j=i:N
        temp_h=zeros(j-i+1,j-i+1);
        temp_h1=zeros(j-i+1,j-i+1);
        for m=1:j-i+1
            for n=m:j-i+1
                if n==m
                    temp_h1(m,n)=0;
                else
                    temp_h1(m,n)=h(m+i-1);%错误在这里
                end
            end
        end
        if j==i
            temp_h=0;
        else
            temp_h(1,:)=temp_h1(1,:);
            for k=2:j-i+1
                temp_h(k,:)=temp_h(k-1,:)+temp_h1(k,:);
            end
        end
        f=c(i)*ones(1,j-i+1)-p(i:j)-pai(i:j)+temp_h(j-i+1,:);
        f=f';
        if i==1
            initial_B(i)=B0+BL;
        end
        lb=zeros(j-i+1,1); ub=d(i:j)';
        A=zeros(j-i+1,j-i+1);b=zeros(j-i+1,1);
        b(j-i+1)=initial_B(i)-s(i);
        for k=1:j-i
            if j<L
                b(k)=initial_B(i)-s(i)-pai(i:k+i-1)*d(i:k+i-1)';
            else
                b(k)=initial_B(i)-s(i)-pai(i:k+i-1)*d(i:k+i-1)'-BL*(1+r)^L;
            end
        end
        A(j-i+1,:)=c(i)*ones(1,j-i+1);
        temp_p=tril(repmat(p(i:j)+pai(i:j),n,1));
        for k=1:j-i
            A(k,:)=c(i)*ones(1,j-i+1)-temp_p(k,:)+temp_h(k,:);  %
        end
        
        %options.algorithm='simplex';   
        [~,fval,exitflag]=linprog(f,A,b,[],[],lb,ub); %,zeros(j-i+1,1),options);%%默认采用内点法
        if exitflag==1
            BB(i,j)=-fval-s(i)-pai(i:j)*d(i:j)'+initial_B(i);
            if BB(i,j)<initial_B(i)-pai(i:j)*d(i:j)' %原因在这里
                BB(i,j)=initial_B(i)-pai(i:j)*d(i:j)';
            end
        else
            BB(i,j)=initial_B(i)-pai(i:j)*d(i:j)';
        end
    end
    [initial_B(i+1),mark(i)]=max(BB(1:i,i));
    if abs(initial_B(i+1)+pai(mark(i):i)*d(mark(i):i)'-initial_B(mark(i)))<1e-4
        mark(i)=-mark(i);
    end
end

%%根据现金流矩阵反推
%back track
j=N;
temp_flag=zeros(1,N);%用来记录生不如死的点
while j>=1
    if mark(j)>0
        x(mark(j))=1;
        j=mark(j);
    else
        temp_flag(-mark(j))=1;
        j=-mark(j);
    end
    j=j-1;
end

sa_demand=zeros(1,N);
i=1;
while i<=N
    if x(i)==1&& i~=N
        for j=i+1:N
            if x(j)~=1&&j<N&&temp_flag(j)~=1
                continue;
            else
                if x(j)==1
                    index=j-1;
                end
                if temp_flag(j)==1
                    index=j-1;
                end
                if j==N&&x(j)~=1&&temp_flag(j)~=1
                    index=N;
                end
                
                temp_h=zeros(index-i+1,index-i+1);
                temp_h1=zeros(index-i+1,index-i+1);
                for m=1:index-i+1
                    for n=m:index-i+1
                        if n==m
                            temp_h1(m,n)=0;
                        else
                            temp_h1(m,n)=h(m+i-1);
                        end
                    end
                end
                if index==i
                    temp_h=0;
                else
                    temp_h(1,:)=temp_h1(1,:);
                    for k=2:index-i+1
                        temp_h(k,:)=temp_h(k-1,:)+temp_h1(k,:);
                    end
                end
                f=c(i)*ones(1,index-i+1)-p(i:index)-pai(i:index)+temp_h(index-i+1,:);
                f=f';
                %if i==1
                    %initial_B=B0;
                %else
                    %initial_B=max(nonzeros(BB(:,i-1)));
                %end
                lb=zeros(index-i+1,1); ub=d(i:index)';
                A=zeros(index-i+2,index-i+1);b=zeros(index-i+2,1);
                b(index-i+2)=initial_B(i)-s(i);
                for k=1:index-i+1
                    if index<L
                        b(k)=initial_B(i)-s(i)-pai(i:k+i-1)*d(i:k+i-1)';
                    else
                        b(k)=initial_B(i)-s(i)-pai(i:k+i-1)*d(i:k+i-1)'-BL*(1+r)^L;
                    end
                end
                A(index-i+2,:)=c(i)*ones(1,index-i+1);
                temp_p=tril(repmat(p(i:index)+pai(i:index),n,1));
                for k=1:index-i+1
                    A(k,:)=c(i)*ones(1,index-i+1)-temp_p(k,:)+temp_h(k,:);
                end
                [sa_demand(i:index),~,~]=linprog(f,A,b,[],[],lb,ub);%%默认采用内点法
                %if exitflag~=1
                    %sa_demand(i:index)=zeros(1,index-i+1);
                    %x(i)=0;
                    %%%出现了一个小问题，为啥exitflag没有返回1
                %end
                i=j;
                break;
            end
        end
    end
    if i==N
        if x(i)==1
            sa_demand(N)=min((initial_B(i)-s(i))/c(N),d(N)); %最后一个点直接算了，不用调算法
        end
        break;
    end
    if x(i)==0
        i=i+1;
    end
end

for i=1:N
    w(i)=d(i)-sa_demand(i);
end

i=1;
while i<=N
    if x(i)==1&& i~=N
        for j=i+1:N
            if x(j)~=1&&j<N
                continue
            else
                if x(j)==1
                    index=j-1;
                end
                if j==N&&x(j)~=1
                    index=j;
                end
                y(i)=sum(sa_demand(i:index));
                i=j;
                break;
            end
        end
    end
    if i==N
        if x(i)==1
            y(N)=sa_demand(N);
        end
        break;
    end
    if x(i)==0
        i=i+1;
    end
end

%启发式调整1
temp_x1=zeros(1,N);%用来记录起始生产的阶段
temp_x2=zeros(1,N);%用来记录生产终止的阶段
temp_i=0;
for i=1:N
    if x(i)==1
        temp_i=i;
    end
    temp_x1(i)=temp_i;
end
temp_i=0;
for i=1:N
    if x(i)==1
        for j=i+1:N
            if x(j)==1
                temp_i=j-1;
                break;
            end
            if j==N
                temp_i=j;
            end
        end
        if i==N
            temp_i=N;
        end
    end
    temp_x2(i)=temp_i;
end

y_increase=zeros(1,N);
for i=1:N-1
    if x(i)==1 
        if initial_B(i)-s(i)-c(i)*y(i)>1e-4 && temp_x2(i)<N && initial_B(temp_x2(i)+1)-s(temp_x2(i)+1)-c(temp_x2(i)+1)*y(temp_x2(i)+1)>1e-4
            if c(i)+sum(h(i:temp_x2(i)))<c(temp_x2(i)+1)
                y_increase(i)=(initial_B(i)-s(i))/c(i)-y(i);%%防止生产太多破产
                y(i)=y(i)+y_increase(i);
                y(temp_x2(i)+1)=y(temp_x2(i)+1)-y_increase(i);
                initial_B(temp_x2(i)+1)=initial_B(temp_x2(i)+1)-y_increase(i)*(c(i)+sum(h(i:temp_x2(i))));
                for j=temp_x2(temp_x2(i)+1)+1:N+1
                    initial_B(j)=initial_B(j)+y_increase(i)*(c(temp_x2(i)+1)-c(i)-sum(h(i:temp_x2(i))));
                end
            end
        end
    end
end

%计算期末资金
I=zeros(1,N);
I(1)=y(1)-(d(1)-w(1));
if I(1)<1e-3
    I(1)=0;
end
for t=2:N
    I(t)=I(t-1)+y(t)-(d(t)-w(t));
    if I(t)<1e-3
        I(t)=0;
    end
end
B(1)=B0+p(1)*(d(1)-w(1))-h(1)*I(1)-y(1)*c(1)-x(1)*s(1)-pai(1)*w(1);
for t=2:N
    if t==L
       B(t)=B(t-1)+p(t)*(d(t)-w(t))-h(t)*I(t)-y(t)*c(t)-x(t)*s(t)-pai(t)*w(t)-BL*(1+r)^L; 
    else
        B(t)=B(t-1)+p(t)*(d(t)-w(t))-h(t)*I(t)-y(t)*c(t)-x(t)*s(t)-pai(t)*w(t);
    end
end 

fprintf('各阶段是否生产:\n');
fprintf('x= ');
for i=1:N
    fprintf('%d\t',x(i));
    if mod(i,20)==0
        fprintf('\n');
    end
end
fprintf('\n');
fprintf('各阶段生产量:\n');
fprintf('y= ');
for i=1:N
    fprintf('%6.2f\t',y(i));
    if mod(i,20)==0
        fprintf('\n');
    end
end
fprintf('\n');
fprintf('各阶段缺货量:\n');
fprintf('w= ');
for i=1:N
    fprintf('%6.2f\t',w(i));
    if mod(i,20)==0
        fprintf('\n');
    end
end
fprintf('\n');
fprintf('各阶段库存水平:\n');
fprintf('I= ');
for i=1:N
    fprintf('%6.2f\t',I(i));
    if mod(i,20)==0
        fprintf('\n');
    end
end
fprintf('\n');
fprintf('各阶段资金量:\n');
fprintf('B= ');
for i=1:N
    fprintf('%6.2f\t',B(i));
    if mod(i,20)==0
        fprintf('\n');
    end
end
fprintf('\n');
end