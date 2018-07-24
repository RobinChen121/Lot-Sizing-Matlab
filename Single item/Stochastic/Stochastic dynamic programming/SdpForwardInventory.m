function SdpForwardInventory
T=3;
state=4;
%a=10;v=1;pai=3;h=1;P=3;Q_max=10;
d=[1,2];
prod_cost1=0;
prod_cost2=@(x)3+2*x;
hold_cost=@(i,x)0.5*(i+x-d(1))+0.5*(i+x-d(2));
salv_cost1=0;
salv_cost2=@(i,x)0.5*2*(i+x-d(1))+0.5*2*(i+x-d(2));

xx1=cell(T,state);%record the optimal ordering quantity in each period
M=100;
CC1=cell(T,state);%record the optimal cost from t to T for a certain state in t

for t=1:T
    for temp_i=1:state
        i=temp_i-1;
        if t==1
            i=1;
        end
        %minc=M;
        CC1{t,i+1}=M*ones(3,1);xx1{t,i+1}=zeros(3,1);        
        for x=max(0,2-i):4-i
            if t==T
                salv_cost=salv_cost2(i,x);%CC_next=0;
            else
                salv_cost=salv_cost1;%CC_next=0.5*CC(t+1,i+x-d(1)+1)+0.5*CC(t+1,i+x-d(2)+1);
            end
            if x==0
                prod_cost=prod_cost1;
            else
                prod_cost=prod_cost2(x);
            end
            CC1{t,i+1}(x-max(0,2-i)+1)=prod_cost+hold_cost(i,x)-salv_cost; xx1{t,i+1}(x-max(0,2-i)+1)=x;
%             temp_c=prod_cost+hold_cost(i,x)-salv_cost;%+CC_next;
%             if temp_c<minc
%                 minc=temp_c;CC(t,i+1)=minc; xx(t,i+1)=x;
%             end
        end  
        if t==1
            break;
        end
    end   
end

%backward to get the solution

CC2=M*ones(T,state);xx2=zeros(T,state);

for temp_i=1:state
    [temp,index]=min(CC1{T,temp_i});
    CC2(T,temp_i)=temp;xx2(T,temp_i)=xx1{T,temp_i}(index);
end

for t=T-1:-1:1
    for temp_i=1:state
        i=temp_i-1;
        if t==1
            i=1;
        end
        for x=max(0,2-i):4-i
            CC1{t,i+1}(x-max(0,2-i)+1)=CC1{t,i+1}(x-max(0,2-i)+1)+0.5*CC2(t+1,i+x-d(1)+1)+0.5*CC2(t+1,i+x-d(2)+1);
        end
        [temp,index]=min(CC1{t,i+1});
        CC2(t,i+1)=temp;xx2(t,i+1)=xx1{t,i+1}(index);
        if t==1
            break;
        end
    end
end

end