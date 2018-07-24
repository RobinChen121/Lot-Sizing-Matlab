function SimulatePiecewise(R)
a=10;v=1;pai=2;h=1;I0=0;B0=14;price=5;
T=8;

meand=[3,4,3,5,4,3,5,4];
sigma=[1,1,1,1,1,1,1,1];

Q=zeros(1,T);B=zeros(1,T);
I=zeros(1,T);delta=zeros(1,T);

if R(1)>0
    Q(1)=min(R(1)-I0,(B0-a)/v);
    delta(1)=0;
end
d=normrnd(meand(1),sigma(1));
I(1)=Q(1)+I0-d;
B(1)=B0+price*min(I0+Q(1),d)-a*delta(1)-v*Q(1)-h*max(I(1),0)-pai*max(-I(1),0);

for t=1:T
    if R(t)>0
        Q(t)=min(R(t)-I(t-1),(B(t-1)-a)/v);
        delta(t)=0;
    end
    d=normrnd(meand(t),sigma(t));
    I(t)=Q(t)+I(t-1)-d;
    B(t)=B(t-1)+price*min(max(I(t-1),0)+Q(t),d+max(-I(t-1),0))-a*delta(t)-v*Q(t)-h*max(I(t),0)-pai*max(-I(t),0);
end

end