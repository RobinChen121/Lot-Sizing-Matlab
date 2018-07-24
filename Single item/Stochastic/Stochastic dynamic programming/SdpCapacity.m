function SdpCapacity

a=10;v=1;pai=2;h=1;I0=0;
Imax=100;Imin=-100;Qmax=20;

meand=[3,4,3,5,4,3,5,4];
sigma=[1,1,1,1,1,1,1,1];

T=8;
meand=meand(1:T);

D=cell(T,1);percent=cell(T,1);
for t=1:T
    D{t}=round(norminv(0.025,meand(t),sigma(t))):round(norminv(0.975,meand(t),sigma(t)));
    percent{t}=normcdf(D{t}+0.5,meand(t),sigma(t))-normcdf(D{t}-0.5,meand(t),sigma(t));
    percent{t}=percent{t}/sum(percent{t});
end

DLength=zeros(T,1);
for t=1:T
    DLength(t)=length(D{t});
end

matrixRecurFull=cell(T,1);
matrixRecurExpection=cell(T,1);
recordOptimalDecision=cell(T,1);

count=0;% count the number of times that sdp hits the inventory bounds
xLow=0;
xUp=Qmax;
xLength=xUp-xLow+1;

tic
for t=1:T
    currentFullRowNum=1;currentExpeRowNum=1; 
    if t==1
        initialI=I0;
    end
    if t<T
        RecurExpeColumNum=3+DLength(t);
    else
        RecurExpeColumNum=3;
    end
    for i=1:length(initialI)
        I=initialI(i);  
        if t==1||(i==1)
            matrixRecurFull{t}=zeros(xLength*DLength(t),4); % x changes related with i, so i is not fixed
            matrixRecurExpection{t}=zeros(xLength,RecurExpeColumNum);
        else
            matrixRecurFull{t}(end+1:end+xLength*DLength(t),:)=zeros(xLength*DLength(t),4);
            matrixRecurExpection{t}(end+1:end+xLength,:)=zeros(xLength,RecurExpeColumNum);
        end
        
        for x=xLow:xUp
            if x>0
                prodCost=a+v*x;
            else
                prodCost=0;
            end
            holdCost=h*max(I+x-D{t}',0);penaCost=-pai*min(I+x-D{t}',0);
            cost=prodCost+holdCost+penaCost;
            INext=(I+x-D{t})';
            INext(logical(INext<Imin))=Imin;
            INext(logical(INext>Imax))=Imax;
            count=count+sum(logical(INext<Imin))+sum(logical(INext>Imax));
            
            matrixRecurFull{t}(currentFullRowNum:currentFullRowNum+DLength(t)-1,1)=initialI(i);
            matrixRecurExpection{t}(currentExpeRowNum,1)=initialI(i);
            matrixRecurFull{t}(currentFullRowNum:currentFullRowNum+DLength(t)-1,2)=x;
            matrixRecurExpection{t}(currentExpeRowNum,2)=x;
            matrixRecurFull{t}(currentFullRowNum:currentFullRowNum+DLength(t)-1,3)=INext;
            matrixRecurFull{t}(currentFullRowNum:currentFullRowNum+DLength(t)-1,4)=cost;
            matrixRecurExpection{t}(currentExpeRowNum,3)=percent{t}*cost;
            if t<T
                matrixRecurExpection{t}(currentExpeRowNum,4:RecurExpeColumNum)=INext; %next IStates
            end
            currentFullRowNum=currentFullRowNum+DLength(t);currentExpeRowNum=currentExpeRowNum+1;
        end
    end
    initialI=sort(unique(matrixRecurFull{t}(:,3)));
end

%backward
cost=zeros(Imax-Imin+1,T);
for t=T:-1:1
    rowLength=size(matrixRecurExpection{t},1); 
    rowIndex=1;
    while rowIndex<rowLength
        I=matrixRecurExpection{t}(rowIndex);
        if t~=T
            for i=rowIndex:rowIndex+xLength-1
                INext=matrixRecurExpection{t}(i,4:4+DLength(t)-1);
                INextIndex=INext-Imin+1;
                matrixRecurExpection{t}(i,3)=matrixRecurExpection{t}(i,3)+percent{t}*cost(INextIndex,t+1);
            end
        end
        [cost(I-Imin+1,t),orderIndex]=min(matrixRecurExpection{t}(rowIndex:rowIndex+xLength-1,3));
        if rowIndex==1
            recordOptimalDecision{t}=[I,matrixRecurExpection{t}(rowIndex+orderIndex-1,2),cost(I-Imin+1,t)];
        else
            recordOptimalDecision{t}(end+1,:)=[I,matrixRecurExpection{t}(rowIndex+orderIndex-1,2),cost(I-Imin+1,t)];
        end
        rowIndex=rowIndex+xLength;
    end   
end

[maxExpecCost,index]=min(matrixRecurExpection{1}(:,3));
optimalOrder=matrixRecurExpection{1}(index,2);
toc
fprintf('count=%d\n',count);
fprintf('optimal totoal expected cost= %f\n',maxExpecCost);
fprintf('expected ordering quantity at the beginning = %f\n', optimalOrder);
end