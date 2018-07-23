function MyWW
% matlab code for WW algorithm
% x             binary variable, whether settting up production
% y             producing quantity
% I              inventory level
% I0            initial inventory level

I0=50;
n=8;
% d=randi([5 25],1,n); 
% S=randi([20 100],1,n); 
% h=randi([5 20],1,n); 
% c=randi([10 25],1,n); 

d=[9,12,9,25,9,20,20,25]; 
S=[100,100,100,100,100,100,100,100];
h=[5,5,5,5,5,5,5,5];
c=10*ones(1,n);

x=zeros(1,n);
y=zeros(1,n);
I=zeros(1,n);

C=1e4*ones(n,n);         % matrix for costs
opt_cost=zeros(n,1);
initial_I=zeros(n,1);
for i=1:n
    if i>1
        opt_cost=min(C(:,i-1));
    end
    for j=i:n
        if I0>=sum(d(1:j))
            if i==1
                C(i,j)=h(1:j)*(I0-cumsum(d(i:j)))';
            else
                C(i,j)=opt_cost+h(i:j)*(I0-cumsum(d(i:j)))';
            end
        else
            initial_I(i)=I0;
            if i>1
                if I0>sum(d(1:i-1))
                    initial_I(i)=I0-sum(d(1:i-1));
                else
                    initial_I(i)=0;
                end
            end
            prod_amount=sum(d(i:j))-initial_I(i);
            temp_I=initial_I(i)+prod_amount;
            h_sum=h(i:j)*(temp_I-cumsum(d(i:j)))';
            if i>1
                C(i,j)=opt_cost+S(i)+h_sum+c(i)*prod_amount;
            else
                C(i,j)=S(i)+h_sum+c(i)*prod_amount;
            end
        end
   end
end

%back track
j=n;
while j>=1
    [~,index]=min(C(:,j));
    if I0<sum(d(1:index))
        x(index)=1;
        if i>1
            y(index)=sum(d(index:j))-initial_I(index);            
        else
            y(index)=sum(d(index:j))-initial_I(index);
        end
    end
    j=index-1;
end
for i=1:n
    if i==1
        I(i)=I0+y(i)-d(i);
    else
        I(i)=I(i-1)+y(i)-d(i);
    end
    fprintf('  %d',I(i));
end
fprintf('\n');
fprintf('producing quantity for each period:\n');
for i=1:n
    fprintf('  %d',y(i));
end
fprintf('\n');
fprintf('whether producing in each period:\n');
for i=1:n
    fprintf('  %d',x(i));
end
fprintf('\n');

end