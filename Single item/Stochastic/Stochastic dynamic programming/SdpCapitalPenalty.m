function SdpCapitalPenalty

%% mail setting
MailAddress = '15011074486@163.com';
password = 'chenzhen881267';
setpref('Internet','E_mail',MailAddress);
setpref('Internet','SMTP_Server','smtp.163.com');
setpref('Internet','SMTP_Username',MailAddress);
setpref('Internet','SMTP_Password',password);
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

%% demands and periods length
meands=[7	7	7	7	7	7
2	3	4	5	6	7
8	7	6	5	4	3
7	6	5	4	5	6
8	5	2	1	2	5
8	4	1	3	1	3					
1	3	8	4	8	7
1	4	7	3	5	8
3	8	4	4	6	2
3	1	5	8	4	4];
T=size(meands,2)-1;D=cell(T,1);percent=cell(T,1);

%% parameters
sampleValues=importdata('640samples.mat');sampleNum=length(sampleValues);
h=1;I0=0;Bmax=500;Bmin=-500;
solutions=zeros(sampleNum,12);

for kk=1:1 %sampleNum %
    a=sampleValues(kk,1);v=sampleValues(kk,2);price=sampleValues(kk,3);pai=sampleValues(kk,4);
    b=sampleValues(kk,5);dd=sampleValues(kk,6);B0=sampleValues(kk,7);
    
    %% bounds for Q and I
    if dd==1
        Imax=30*ones(1,6);
        Imin=-50*ones(1,6);
        Qmax=23*ones(1,6);
    end
    if dd==2
        Imax=30*ones(1,6);
        Imin=-50*ones(1,6);
        Qmax=23*ones(1,6);
    end
    if dd==3
        Imax=30*ones(1,6);
        Imin=-50*ones(1,6);
        Qmax=23*ones(1,6);
    end
    if dd==4
        Imax=50*ones(1,6);
        Imin=-50*ones(1,6);
        Qmax=23*ones(1,6);
    end
    if dd==5
        Imax=50*ones(1,6);
        Imin=-50*ones(1,6);
        Qmax=23*ones(1,6);
    end
    if dd==6
    Imax=[20,30,30,20,10,8];
    Imin=[-10,-15,-20,-23,-25,-30];
    Qmax=[20,15,10,10,8,6];
    %Qmax=30*ones(1,6);
    end
    if dd==7
    Imax=[20,40,60,60,30,10];
    Imin=[-3,-9,-25,-30,-40,-45];
    Qmax=[20,20,25,25,20,10];
    end
    if dd==8
    Imax=[15,35,55,50,30,10];
    Imin=[-3,-10,-15,-25,-40,-50];
    Qmax=[15,20,20,20,25,10];
    end
    if dd==9
    Imax=[18,35,55,40,25,12];
    Imin=[-6,-15,-20,-22,-30,-38];
    Qmax=[18,20,18,15,12,10];
    end
    if dd==10
    Imax=[13,33,55,45,25,10];
    Imin=[-6,-10,-18,-25,-30,-40];
    Qmax=[13,20,22,20,15,8];
    end
    meand=meands(dd,:);
    
%% demands possibilities
for t=1:T
    D{t}=poissinv(0,meand(t)):poissinv(0.99,meand(t));
    percent{t}=poisspdf(D{t},meand(t));
    percent{t}=percent{t}/sum(percent{t});
end

%% matrixes for recursion
DLength=zeros(T,1);
for t=1:T
    DLength(t)=length(D{t});
end
matrixRecurFull=cell(T,1);
matrixRecurExpection=cell(T,1);
recordCapital=cell(T,1);
recordOptimalDecision=cell(T,1);

%% forward recursion
tic
count1=0;count2=0;count3=0;
for t=1:T  
    currentFullRowNum=1;currentExpeRowNum=1; 
    if t==1
        initialIB=[I0,B0]; % initial inventory and capital marks
    end
    initialI=unique(initialIB(:,1)); % duplicate checking, initial inventory capital marks
    if t<T
        RecurExpeColumNum=4+2*DLength(t);
    else
        RecurExpeColumNum=4;
    end
    for i=1:length(initialI)
        initialB=initialIB(logical(initialIB(:,1)==initialI(i)),2);% initial capital marks, one initialI may respond to several initialB
        initialB=sort(unique(initialB)); %duplicate checking and sorting
        for j=1:length(initialB)   
            I=initialI(i);
            B=initialB(j);    
            xLow=0;
            if b>10
                xUp=max(0,round((B-a)/v));
                xUp=max(0,min(Qmax(t),xUp));
            else
                xUp=Qmax(t);
            end
            xLength=xUp-xLow+1; 
            if t==1||(i==1&&j==1)  
                matrixRecurFull{t}=zeros(xLength*DLength(t),6); % x changes related with i, so i is not fixed
                matrixRecurExpection{t}=zeros(xLength,RecurExpeColumNum);
            else
                matrixRecurFull{t}(end+1:end+xLength*DLength(t),:)=zeros(xLength*DLength(t),6);
                matrixRecurExpection{t}(end+1:end+xLength,:)=zeros(xLength,RecurExpeColumNum);
            end  
            
            for x=xLow:xUp
                if x>0
                    prodCost=a+v*x;
                else
                    prodCost=0;
                end
                revenue=price*min(D{t}'-min(I,0),x+max(I,0));holdCost=h*max(I+x-D{t}',0);penaCost=-pai*min(I+x-D{t}',0);
                capital=B+revenue-holdCost-penaCost-prodCost;
                if b<=10
                    if B+1e-1<prodCost
                        capital=capital-b*(prodCost-B);
                    end
                end
                capital=round(capital); % round to the nearest integer
                count1=count1+sum(logical(capital<Bmin+1e-1));  
                capital(logical(capital<Bmin))=Bmin;%
                %[~,index1]=min(abs(capital-BStates));   
                INext=(I+x-D{t})';
                count2=count2+sum(logical(INext<Imin(t)-1e-1));  
                INext(logical(INext<Imin(t)))=Imin(t);
                count3=count3+sum(logical(INext>Imax(t)+1e-1));  
                INext(logical(INext>Imax(t)))=Imax(t);
                                        
                matrixRecurFull{t}(currentFullRowNum:currentFullRowNum+DLength(t)-1,1)=I;
                matrixRecurExpection{t}(currentExpeRowNum,1)=I;
                matrixRecurFull{t}(currentFullRowNum:currentFullRowNum+DLength(t)-1,2)=B;
                matrixRecurExpection{t}(currentExpeRowNum,2)=B;
                matrixRecurFull{t}(currentFullRowNum:currentFullRowNum+DLength(t)-1,3)=x;
                matrixRecurExpection{t}(currentExpeRowNum,3)=x;
                matrixRecurFull{t}(currentFullRowNum:currentFullRowNum+DLength(t)-1,4)=INext;
                matrixRecurFull{t}(currentFullRowNum:currentFullRowNum+DLength(t)-1,5)=capital;
                matrixRecurFull{t}(currentFullRowNum:currentFullRowNum+DLength(t)-1,6)=capital-B;%capital increase
                matrixRecurExpection{t}(currentExpeRowNum,4)=percent{t}*capital-B;
                
                if t<T
                    matrixRecurExpection{t}(currentExpeRowNum,5:2:RecurExpeColumNum-1)=INext; %next IStates
                    matrixRecurExpection{t}(currentExpeRowNum,6:2:RecurExpeColumNum)=capital; %next BStates
                end
                currentFullRowNum=currentFullRowNum+DLength(t);currentExpeRowNum=currentExpeRowNum+1;
            end
        end
    end
    initialIB=matrixRecurFull{t}(:,4:5);
end
                          
%% backward to get the solution
for t=T:-1:1   
    rowLength=size(matrixRecurExpection{t},1);
    rowNum=Imax(t)-Imin(t)+1;
    columnNum=Bmax-Bmin+1;
    capital=zeros(rowNum,columnNum); % record capital
    rowIndex=1;
    while rowIndex<=rowLength
        I=matrixRecurExpection{t}(rowIndex,1);B=matrixRecurExpection{t}(rowIndex,2);
        xLow=0;
        if b>10
            xUp=max(0,round((B-a)/v));
            xUp=max(0,min(xUp,Qmax(t)));
        else
            xUp=Qmax(t);
        end
        xLength=xUp-xLow+1;
        if t~=T
            for i=rowIndex:rowIndex+xLength-1
                INext=matrixRecurExpection{t}(i,5:2:4+2*DLength(t));
                BNext=matrixRecurExpection{t}(i,6:2:4+2*DLength(t));
                INextIndex=INext-Imin(t+1)+1;
                BNextIndex=BNext-Bmin+1;
                tempCapitial=zeros(DLength(t),1);
                for j=1:DLength(t)
                    tempCapitial(j)=recordCapital{t+1}(INextIndex(j),BNextIndex(j));     
                end
                matrixRecurExpection{t}(i,4)=matrixRecurExpection{t}(i,4)+percent{t}*tempCapitial;
            end
        end
        [capital(I-Imin(t)+1,B-Bmin+1),orderIndex]=max(matrixRecurExpection{t}(rowIndex:rowIndex+xLength-1,4));
        if rowIndex==1
            recordOptimalDecision{t}=[I,B,matrixRecurExpection{t}(rowIndex+orderIndex-1,3),capital(I-Imin(t)+1,B-Bmin+1)];
        else
            recordOptimalDecision{t}(end+1,:)=[I,B,matrixRecurExpection{t}(rowIndex+orderIndex-1,3),capital(I-Imin(t)+1,B-Bmin+1)];
        end
        rowIndex=rowIndex+xLength;
    end
    recordCapital{t}=capital;
end
[maxExpecProfit,index]=max(matrixRecurExpection{1}(:,4));
optimalOrder=matrixRecurExpection{1}(index,3);
fprintf('count1=%d\n',count1);
fprintf('count2=%d\n',count2);
fprintf('count3=%d\n',count3);
fprintf('optimal expected final capital= %f\n',maxExpecProfit+B0);
%fprintf('expected ordering quantity at the beginning = %f\n', optimalOrder);
solution1=maxExpecProfit+B0;

%% simulation
N=100000;deviError=zeros(N,1);
string=['benchmark',num2str(dd),'.mat'];
D2=importdata(string);
for i=1:N
    Q=zeros(1,T);B=zeros(1,T);I=zeros(1,T);delta=zeros(1,T);
    Q(1)=optimalOrder;
    if Q(1)>1e-1
        delta(1)=1;
    end
    realD=D2(i,:);   
    I(1)=I0+Q(1)-realD(1);B(1)=B0+price*min(I0+Q(1),realD(1))-a*delta(1)-v*Q(1)-h*max(I(1),0)-pai*max(-I(1),0);
    if b<10
        if B0+1e-1<a*delta(1)+v*Q(1)
            B(1)=round(B(1)-(a*delta(1)+v*Q(1)-B0)*b);
        end
    end
    for t=2:T
        areaIB=recordOptimalDecision{t}(logical(recordOptimalDecision{t}(:,1)==I(t-1)),1:3);
        try
        Q(t)=areaIB(areaIB(:,2)==B(t-1),3);
        catch
            Q;
        end
        if b>10
            Q(t)=min(Q(t),round((B(t-1)-a)/v));
            Q(t)=max(Q(t),0);
        end
        if Q(t)>1e-1
            delta(t)=1;   
        end
        I(t)=I(t-1)+Q(t)-realD(t);
        B(t)=B(t-1)+price*min(max(I(t-1),0)+Q(t),realD(t)+max(-I(t-1),0))-a*delta(t)-v*Q(t)-h*max(I(t),0)-pai*max(-I(t),0);
        if b<10
            if B(t-1)+1e-1<a*delta(t)+v*Q(t)
                B(t)=round(B(t)-(a*delta(t)+v*Q(t)-B(t-1))*b);
            end
        end
    end
    deviError(i)=B(T);
end
fprintf('%.4f\n',sum(deviError)/N);
solution2=sum(deviError)/N;
toc

%% record solution
solutions(kk,1)=a;solutions(kk,2)=v;solutions(kk,3)=pai;solutions(kk,4)=h;
solutions(kk,5)=price;solutions(kk,6)=b;solutions(kk,7)=count1;
solutions(kk,8)=count2;solutions(kk,9)=count3;ttime=toc;solutions(kk,10)=ttime;
solutions(kk,11)=solution1;solutions(kk,12)=solution2;
if kk/8==0 %
    string1=['test-SDP-',num2str(kk),'.xlsx'];
    xlswrite(string1,solutions);
    sendmail('15011074486@163.com',['sdp-results',num2str(dd)],'',{string1});
end

end

string=['test-SDP','.xlsx'];
xlswrite(string,solutions);
sendmail('15011074486@163.com',['sdp-results',num2str(dd)],'',{string});

end