function SdpBackwardInventory
a=10;v=1;pai=2;h=1;
I0=0;Qmax=25;Imin=-100;Imax=100;
d=ones(1,5)*10;
T=length(d);

D=cell(T,1);percent=cell(T,1);
for t=1:T
    D{t}=poissinv(0.025,d(t)):poissinv(0.975,d(t));
    percent{t}=poisspdf(D{t},d(t));
end

prod_cost1=0;
prod_cost2=@(x)a+v*x;

count1=0;count2=0;
costFull=cell(T,1);
costExpect=cell(T,1);
for t=1:T
    if t==1
        i=I0;
    else
        i=unique(initialI);
    end
    rowNum1=1;
    rowNum2=1;
    Ilength=length(i);
    xLength=Qmax+1;
    costFull{t}=zeros(xLength*length(D{t}),3);
    costExpect{t}=zeros(Ilength*xLength,3+length(D{t}));
    for k=1:Ilength
         I=i(k);      
        for x=0:Qmax
            if x==0
                prod_cost=prod_cost1;
            else
                prod_cost=prod_cost2(x);
            end
            if t==1
                prod_cost=v*I0;
            end
            hold_cost=h*max(I+x-D{t}',0);pena_cost=-pai*min(I+x-D{t}',0);
            costFull{t}(rowNum1:rowNum1+length(D{t})-1,1)=I;
            costFull{t}(rowNum1:rowNum1+length(D{t})-1,2)=x;
            costExpect{t}(rowNum2,1)=I;
            costExpect{t}(rowNum2,2)=x;
            costFull{t}(rowNum1:rowNum1+length(D{t})-1,3)=prod_cost+hold_cost+pena_cost;
            costExpect{t}(rowNum2,3)=percent{t}*(prod_cost+hold_cost+pena_cost);
            INext=I+x-D{t};
            costFull{t}(rowNum1:rowNum1+length(D{t})-1,4)=INext;
            INext(logical(INext<Imin))=Imin;
            count1=count1+sum(logical(INext>Imax+1e-1));
            INext(logical(INext>Imax))=Imax;
            count2=count2+sum(logical(INext<Imin+1e-1));
            costExpect{t}(rowNum2,4:3+length(D{t}))=INext; 
            rowNum2=rowNum2+1;rowNum1=rowNum1+length(D{t});
        end 
    end
    initialI=costFull{t}(:,4);
end

cost=10000*ones(Imax-Imin+1,t);
for t=T:-1:1
    rowLength=size(costExpect{t},1);
    rowIndex=1;
    xLength=Qmax+1;
    while rowIndex<=rowLength
        I=costExpect{t}(rowIndex,1);
        if t~=T
            for i=rowIndex:rowIndex+xLength-1
                INext=costExpect{t}(i,4:3+length(D{t}));
                tempCost=zeros(length(D{t}),1);
                for j=1:length(D{t})
                    tempCost(j)=cost(INext(j)-Imin+1,t+1);
                end
                costExpect{t}(i,3)=costExpect{t}(i,3)+percent{t}*tempCost;
            end
        end
        cost(I-Imin+1,t)=min(costExpect{t}(rowIndex:rowIndex+xLength-1,3));
        rowIndex=rowIndex+xLength;
    end
end

[maxExpecCost,index]=min(costExpect{1}(:,3));
optimalOrder=costExpect{1}(index,2);
fprintf('optimal expected total cost= %f\n',maxExpecCost);
fprintf('expected ordering quantity at the beginning = %f\n', optimalOrder);

end