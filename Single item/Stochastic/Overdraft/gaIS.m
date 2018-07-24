function solutions= gaIS

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


%% Start with the default options
global a v pai h I0 B0 price b T;
a=10;v=1;pai=2;h=1;I0=0;B0=5;price=5;b=0.2;
T=3;

%% Default seeting of ga
nvars=2*T;
%intcon=[];
intcon=(1:2*T);
options = optimoptions('ga');
options = optimoptions(options,'Display', 'off');

%% Run ga
Num=1;
solutions=zeros(Num,6+2*T+2); 
for k=1:Num
%     ub=[meands(dNum,:),3*max(meands(dNum,:))*ones(1,T)];
%     if b<10&&B0>10
%         lb=zeros(1,2*T);
%     else
%         lb=[-2*max(meands(dNum,:))*ones(1,T),zeros(1,T)];
%     end
    ub=8*ones(1,2*T);lb=-8*ones(1,2*T);
    tic
    [x,fval] = ...
        ga(@IS,nvars,[],[],[],[],lb,ub,[],intcon,options);
    
    solutions(k,1:6)=[B0,a,v,pai,price,b];
    solutions(k,7:2*T+6)=x;
    solutions(k,2*T+7)=-fval;
    fprintf('iniResult=%.2f\n',-fval);
    
    
    %% Simultion
    x=[0,7,0,5,3,3];
    I1=x(1:T);
    S1=x(T+1:2*T);
    %string=['benchmark',num2str(dNum),'.mat'];
    string='benchmark-example.mat';
    D=importdata(string);
    D=[2,1,2;2,1,1;2,2,2;1,1,2;1,2,1];
    
    N=length(D);
    simuResult=zeros(N,1);
    for i=1:N
        B=zeros(1,T);I=zeros(1,T);delta=zeros(1,T);
        if I0>=I1(1)
            tempQ=0;
        else
            if b>10
                tempQ=max(0,min((B0-a)/v,S1(1)-I0));
            else
                tempQ=max(0,S1(1)-I0);
            end
        end
        if tempQ>1e-1
            delta(1)=1;
        end
        realD=D(i,:);
        I(1)=I0+tempQ-realD(1);
        B(1)=B0+price*min(I0+tempQ,realD(1))-a*delta(1)-v*tempQ-h*max(I(1),0)-pai*max(-I(1),0);
%         if b<10&&B0<a*delta(1)+v*tempQ
%             B(1)=B(1)-b*(B0-a*delta(1)-v*tempQ);
%         end
        if b<10&&B(1)<0
            B(1)=B(1)+b*B(1);
        end
        for t=2:T
            if I(t-1)>=I1(t)
                tempQ=0;
            else
                if b>10
                    tempQ=max(0,min((B(t-1)-a)/v,S1(t)-I(t-1)));
                else
                    tempQ=max(0,S1(t)-I(t-1));
                end
            end
            if tempQ>1e-1
                delta(t)=1;
            end
            I(t)=I(t-1)+tempQ-realD(t);
            B(t)=B(t-1)+price*min(max(I(t-1),0)+tempQ,realD(t)+max(-I(t-1),0))-a*delta(t)-v*tempQ-h*max(I(t),0)-pai*max(-I(t),0);
%             if b<10&&B(t-1)<a*delta(t)+v*tempQ
%                 B(t)=B(t)-b*(B(t-1)-a*delta(t)-v*tempQ);
%             end
            if b<10&&B(t)<0
                B(t)=B(t)+b*B(t);
            end
        end
        simuResult(i)=B(T);
    end
    toc %
    %sendmail('15011074486@163.com',['time',num2str(k)],num2str(ttime));
    solutions(k,2*T+8)=mean(simuResult);
    fprintf('simuResult=%.2f\n',mean(simuResult));
    fprintf('k=%d\n',k);
end
%     string='solutionIS.xlsx';
%     xlswrite(string,solutions);
%     sendmail('15011074486@163.com','solutionIS','',{string});
end



%% IS
function result=IS(input)
global a v pai h I0 B0 price T b;

Ib=input(1:T);
Sb=input(T+1:2*T);


%% load demand
string='example-1000.mat';
D=importdata(string);

%% Expectation capital
N=length(D);
results=zeros(N,1);
for i=1:N
    B=zeros(1,T);I=zeros(1,T);delta=zeros(1,T);
    if I0>=Ib(1)
        tempQ=0;
    else
        if b>10
            tempQ=max(0,min((B0-a)/v,Sb(1)-I0));
        else
            tempQ=max(0,Sb(1)-I0);
        end
    end
    if tempQ>1e-1
        delta(1)=1;
    end
    realD=D(i,:);
    I(1)=I0+tempQ-realD(1);
    B(1)=B0+price*min(I0+tempQ,realD(1))-a*delta(1)-v*tempQ-h*max(I(1),0)-pai*max(-I(1),0);
%     if b<10&&B0<a*delta(1)+v*tempQ
%         B(1)=B(1)-b*(B0-a*delta(1)-v*tempQ);
%     end
    if b<10&&B(1)<0
        B(1)=B(1)+b*B(1);
    end
    for t=2:T
        if I(t-1)>=Ib(t)
            tempQ=0;
        else
            if b>10
                tempQ=max(0,min((B(t-1)-a)/v,Sb(t)-I(t-1)));
            else
                tempQ=max(0,Sb(t)-I(t-1));
            end
        end
        if tempQ>1e-1
            delta(t)=1;
        end
        I(t)=I(t-1)+tempQ-realD(t);
        B(t)=B(t-1)+price*min(max(I(t-1),0)+tempQ,realD(t)+max(-I(t-1),0))-a*delta(t)-v*tempQ-h*max(I(t),0)-pai*max(-I(t),0);
%         if b<10&&B(t)<a*delta(1)+v*tempQ
%             B(t)=B(t)-b*(B(t-1)-a*delta(t)-v*tempQ);
%         end
        if b<10&&B(t)<0
            B(t)=B(t)+b*B(t);
        end
    end
    results(i)=B(T);
end
result=-mean(results);

end
