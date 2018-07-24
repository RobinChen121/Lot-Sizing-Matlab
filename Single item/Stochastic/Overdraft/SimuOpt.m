function SimuOpt
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

%% parameters
sampleNum=1;
global a v pai h I0 B0 price b T sceNum;
a=10;v=1;pai=2;h=1;I0=0;B0=5;price=5;b=0.2;sceNum=1000;
T=3;


%% forward computation
result=zeros(sampleNum,11);
for k=1:1
    realD=importdata('benchmark-example.mat');
    %realD=[2,1,2;2,1,1;2,2,2;1,1,2;1,2,1];
    fullSampleNum=length(realD)/1000;
    finalB=zeros(fullSampleNum,1);
    recordQ=zeros(fullSampleNum,3);
    tic
    for n=1:fullSampleNum
        tRealD=realD(n,:);
        for t=1:T  
            if t==1
                iniI=I0;iniB=B0;
            end
            %Q=SampleQ(iniI,iniB,t);
            Q=round(SampleQ(iniI,iniB,t));
            I=Q+iniI-tRealD(t); B=iniB+price*min(Q+max(iniI,0),tRealD(t)+max(-iniI,0))-h*max(0,I)-pai*max(0,-I);
            if Q>=1e-1
                B=B-a-v*Q;
            end
            if B<0
                B=B+b*B;
            end
            iniI=I;iniB=B;recordQ(n,t)=Q;
        end
        finalB(n)=B;
    end
    ttime=toc
    
    fprintf('expected final capital: %.2f\n',mean(finalB));
    
    %% record results
    result(k,1)=a;result(k,2)=v;result(k,3)=price;result(k,4)=h;result(k,5)=pai;
    result(k,6)=b;result(k,7)=dPattern;result(k,8)=I0;result(k,9)=B0;result(k,10)=mean(finalB);
    result(k,11)=ttime;
end

end

%% use sampling to get optimal Q
function qOpt=SampleQ(iniI,iniB,tStart)
global a v pai h price b T sceNum;

%% parameter values    
Qmax=30;

%% import samples
string='example-1000.mat';
fullD=importdata(string);
D=fullD(:,tStart:T);tLength=T-tStart+1;

%% select Q, using Golden search
capitalMax=-1000;qOpt=0;
for t=1:tLength  
    % find optimal q for every t
    eps=0.1;q1=0;q2=0.382*Qmax;q3=0.618*Qmax;q4=Qmax;
    while abs(q2-q3)>eps
        B2=CapitalForEachQ(q2,tStart,tStart+t-1);B3=CapitalForEachQ(q3,tStart,tStart+t-1);
        if B2<B3
            q1=q2;q2=q3;q3=q1+0.618*(q4-q1);
        else
            q4=q3;q3=q2;q2=q1+0.382*(q4-q1);
        end
    end
    q=(q2+q3)/2;endCapital=CapitalForEachQ(q,tStart,tStart+t-1);
    if endCapital<CapitalForEachQ(0,tStart,tStart+t-1)
        endCapital=CapitalForEachQ(0,tStart,tStart+t-1);q=0;
    end
    if endCapital/t>capitalMax
        capitalMax=endCapital;qOpt=q;
    else
        return;
    end
end

    %% compute capital for each q
    function expeFinalCapital=CapitalForEachQ(q,tStart,tEnd)
        I=zeros(sceNum,tEnd-tStart+1);B=zeros(sceNum,tEnd-tStart+1);
        if q<1e-1
            for t1=1:tEnd-tStart+1
                if t1==1
                    I(:,t1)=iniI-D(:,t1);
                    B(:,t1)=iniB+price*min(max(iniI,0),D(:,t1)+max(-iniI,0))-h*max(0,I(:,t1))-pai*max(0,-I(:,t1));
                    B(:,t1)=B(:,t1)-b*max(-iniB,0);
                else
                    I(:,t1)=I(:,t1-1)-D(:,t1);
                    B(:,t1)=B(:,t1-1)+price*min(max(iniI,0),D(:,t1)+max(-iniI,0))-h*max(0,I(:,t1))-pai*max(0,-I(:,t1));
                    B(:,t1)=B(:,t1)-b*max(-B(:,t1),0);
                end
            end
        else
            for t1=1:tEnd-tStart+1
                if t1==1
                    I(:,t1)=iniI+q-D(:,t1);
                    B(:,t1)=iniB+price*min(q+max(iniI,0),D(:,t1)+max(-iniI,0))-h*max(0,I(:,t1))-pai*max(0,-I(:,t1))-a-v*q;
                    B(:,t1)=B(:,t1)-b*max(-B(:,t1),0); % max(a+v*q-B(:,t1-1),0)
                else
                    I(:,t1)=I(:,t1-1)+q-D(:,t1);
                    B(:,t1)=B(:,t1-1)+price*min(q+max(iniI,0),D(:,t1)+max(-iniI,0))-h*max(0,I(:,t1))-pai*max(0,-I(:,t1))-a-v*q;
                    B(:,t1)=B(:,t1)-b*max(-B(:,t1),0); % max(a+v*q-B(:,t1-1),0)
                end
            end
        end
        expeFinalCapital=mean(B(:,tEnd-tStart+1));
    end

end

