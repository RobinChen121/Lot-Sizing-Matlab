function SdpOverdraftPenalty

%parameters
a=10;v=1;pai=0;h=1;
I0=0;
B0=30;
price=5;rate=0;
Imax=100;Imin=0;
Bmax=200;Bmin=-80;
%Qmax=30;

%demands and their probabilities
 meand=[8, 10, 5, 10];

T=length(meand);
 D=cell(T,1);percent=cell(T,1);
 
% % for t=1:T
% %     D{t}=round(norminv(0.025,meand(t),sigma(t))):round(norminv(0.975,meand(t),sigma(t)));
% %     percent{t}=normcdf(D{t}+0.5,meand(t),sigma(t))-normcdf(D{t}-0.5,meand(t),sigma(t));
% %     percent{t}=percent{t}/sum(percent{t});
% % end
for t=1:T
    D{t}=poissinv(0,meand(t)):poissinv(0.99,meand(t));
    percent{t}=poisspdf(D{t},meand(t));
    percent{t}=percent{t}/sum(percent{t});
end


% states
Istep=1;Bstep=1;
IStates=(Imin:Istep:Imax); %inventory states
BStates=(Bmin:Bstep:Bmax); % capital states

DLength=zeros(T,1);
for t=1:T
    DLength(t)=length(D{t});
end

tic
matrixRecurFull=cell(T,1);
matrixRecurExpection=cell(T,1);
recordCapital=cell(T,1);
recordOptimalDecision=cell(T,1);

% first, forward recursion
for t=1:T     
    currentFullRowNum=1;currentExpeRowNum=1; 
    if t==1
        [~,B0Index]=min(abs(B0-BStates));
        [~,I0Index]=min(abs(I0-IStates));
        initialIBIndex=[I0Index,B0Index]; % initial inventory and capital marks
    end
    initialI=unique(initialIBIndex(:,1)); % duplicate checking, initial inventory capital marks
    for i=1:length(initialI)
        initialB=initialIBIndex(logical(initialIBIndex(:,1)==initialI(i)),2);% initial capital marks, one initialI may respond to several initialB
        initialB=sort(unique(initialB)); %duplicate checking and sorting
        for j=1:length(initialB)   
            I=IStates(initialI(i));B=BStates(initialB(j));    
            %xLow=max(0,max(D{t}-I+Imin));
            %xUp=Imax-I; % the upper bound for x  min(Q_max,I_max-I)
            xLow=0;
            xUp=max(0,round((B-a)/v));
            %xUp=min(xUp,Qmax);
            %xUp=Qmax;
            xLength=xUp-xLow+1; 
            if t<T
                RecurExpeColumNum=4+2*DLength(t);
            else
                RecurExpeColumNum=4;
            end
            if t==1||(i==1&&j==1)  
                matrixRecurFull{t}=zeros(xLength*DLength(t),6); % x changes related with i, so i is not fixed
                matrixRecurExpection{t}=zeros(xLength,RecurExpeColumNum);
            end
            if t~=1&&(i~=1||j~=1)
                matrixRecurFull{t}(end+1:end+xLength*DLength(t),:)=zeros(xLength*DLength(t),6);
                matrixRecurExpection{t}(end+1:end+xLength,:)=zeros(xLength,RecurExpeColumNum);
            end  
            
            for x=xLow:xUp
                if x>0
                    prodCost=a+v*x;
                else
                    prodCost=0;
                end
                revenue=price*min(D{t}'-min(I,0),x+max(I,0));
                holdCost=h*max(I+x-D{t}',0);penaCost=-pai*min(I+x-D{t}',0);
                capital=B+revenue-holdCost-penaCost-prodCost;
                capital=round(capital+rate*min(capital,0)); % round to the nearest integer
                capital(logical(capital<Bmin))=Bmin;%
                INext=max(0,(I+x-D{t})');
                INextIndex=zeros(DLength(t),1);
                BNextIndex=zeros(DLength(t),1);
                capital(logical(capital>Bmax))=Bmax;% bounds for B states
                for k=1:length(capital)
                    [~,BIndex]=min(abs(capital(k)-BStates));
                    [~,IIndex]=min(abs(INext(k)-IStates));
                    capital(k)=BStates(BIndex);INext(k)=IStates(IIndex);
                    INextIndex(k)=IIndex;BNextIndex(k)=BIndex;
                end
                
                matrixRecurFull{t}(currentFullRowNum:currentFullRowNum+DLength(t)-1,1)=initialI(i);
                matrixRecurExpection{t}(currentExpeRowNum,1)=initialI(i);
                matrixRecurFull{t}(currentFullRowNum:currentFullRowNum+DLength(t)-1,2)=initialB(j);
                matrixRecurExpection{t}(currentExpeRowNum,2)=initialB(j);
                matrixRecurFull{t}(currentFullRowNum:currentFullRowNum+DLength(t)-1,3)=x;
                matrixRecurExpection{t}(currentExpeRowNum,3)=x;
                matrixRecurFull{t}(currentFullRowNum:currentFullRowNum+DLength(t)-1,4)=INextIndex;
                matrixRecurFull{t}(currentFullRowNum:currentFullRowNum+DLength(t)-1,5)=BNextIndex;
                matrixRecurFull{t}(currentFullRowNum:currentFullRowNum+DLength(t)-1,6)=capital-B;%capital increase
                matrixRecurExpection{t}(currentExpeRowNum,4)=percent{t}*capital-B;
                if t<T
                    matrixRecurExpection{t}(currentExpeRowNum,5:2:RecurExpeColumNum-1)=INextIndex; %next IStates
                    matrixRecurExpection{t}(currentExpeRowNum,6:2:RecurExpeColumNum)=BNextIndex; %next BStates
                end
                currentFullRowNum=currentFullRowNum+DLength(t);currentExpeRowNum=currentExpeRowNum+1;
            end
        end
    end
    initialIBIndex=matrixRecurFull{t}(:,4:5);
end

                          
% backward to get the solution
for t=T:-1:1   
    StatesLength=size(matrixRecurExpection{t},1);
    rowNum=matrixRecurExpection{t}(end,1);
    columnNum=matrixRecurExpection{t}(end,2);
    capital=-10000*ones(rowNum,columnNum); % record capital
    rowIndex=1;
    while rowIndex<StatesLength
        IIndex=matrixRecurExpection{t}(rowIndex,1);BIndex=matrixRecurExpection{t}(rowIndex,2);
%         xLow=max(0,max(D{t})-IStates(IIndex)+Imin);
%         xUp=Imax-IStates(IIndex); %min(Q_max,I_max-I_states(IIndex));
        xLow=0;
        xUp=max(0,round((BStates(BIndex)-a)/v));
        %xUp=min(xUp,Qmax);
        %xUp=Qmax;
        xLength=xUp-xLow+1;
        if t~=T
            for i=rowIndex:rowIndex+xLength-1
                for j=1:DLength(t)
                    matrixRecurExpection{t}(i,4)=matrixRecurExpection{t}(i,4)+percent{t}(j)*recordCapital{t+1}(matrixRecurExpection{t}(i,5+(j-1)*2),matrixRecurExpection{t}(i,6+(j-1)*2));
                end
            end
        end
        [capital(IIndex,BIndex),orderIndex]=max(matrixRecurExpection{t}(rowIndex:rowIndex+xLength-1,4));
        if rowIndex==1
            recordOptimalDecision{t}=[IStates(IIndex),BStates(BIndex),matrixRecurExpection{t}(rowIndex+orderIndex-1,3),capital(IIndex,BIndex)];
        else
            recordOptimalDecision{t}(end+1,:)=[IStates(IIndex),BStates(BIndex),matrixRecurExpection{t}(rowIndex+orderIndex-1,3),capital(IIndex,BIndex)];
        end
        %recordOptimalDecision{t}=[recordOptimalDecision{t};IIndex,BIndex,matrixRecurExpection(rowIndex+orderIndex-1,3),capital(IIndex,BIndex)];
        rowIndex=rowIndex+xLength;
    end
    recordCapital{t}=capital;%recordOptimalDecision{t}=recordOptimalDecision{t};
end

[maxExpecProfit,index]=max(matrixRecurExpection{1}(:,4));
optimalOrder=matrixRecurExpection{1}(index,3);
toc
fprintf('optimal expected final capital= %f\n',maxExpecProfit+B0);
fprintf('expected ordering quantity at the beginning = %f\n', optimalOrder);


% realD=meand;
% % D=[2,1,2;2,1,1;2,2,2;1,1,2;1,2,1];
% % realD=D(5,:);
% Q=zeros(1,T);B=zeros(1,T);I=zeros(1,T);delta=zeros(1,T);
% Q(1)=optimalOrder;
% if Q(1)>0
%     delta(1)=1;
% end
% I(1)=I0+Q(1)-realD(1);
% B(1)=B0+price*min(I0+Q(1),realD(1))-a*delta(1)-v*Q(1)-h*max(I(1),0)-pai*max(-I(1),0);
% if B(1)<0
%     B(1)=round(B(1)+rate*B(1));
% end
% for t=2:T
%     areaIB=recordOptimalDecision{t}(logical(recordOptimalDecision{t}(:,1)==I(t-1)),1:3);
%     Q(t)=areaIB(areaIB(:,2)==B(t-1),3);
%     if Q(t)>0
%         delta(t)=1;
%     end
%     I(t)=I(t-1)+Q(t)-realD(t);
%     B(t)=B(t-1)+price*min(max(I(t-1),0)+Q(t),realD(t)+max(-I(t-1),0))-a*delta(t)-v*Q(t)-h*max(I(t),0)-pai*max(-I(t),0);
%     if B(t)<0
%         B(t)=round(B(t)+rate*B(t));
%     end
% end
% I(1)=I0+Q(1)-realD(1);B(1)=B0+price*min(I0+Q(1),realD(1))-a*delta(1)-v*Q(1)-h*max(I(1),0)-pai*max(-I(1),0);
% if B(1)<0
%     B(1)=B(1)+rate*B(1);
% end
% for t=2:T
%     I(t)=I(t-1)+Q(t)-realD(t);
%     B(t)=B(t-1)+price*min(max(I(t-1),0)+Q(t),realD(t)+max(-I(t-1),0))-a*delta(t)-v*Q(t)-h*max(I(t),0)-pai*max(-I(t),0);
%     if B(t)<0
%         B(t)=B(t)+rate*B(t);
%     end
% end
% display(price-v);
% display(Q);
end