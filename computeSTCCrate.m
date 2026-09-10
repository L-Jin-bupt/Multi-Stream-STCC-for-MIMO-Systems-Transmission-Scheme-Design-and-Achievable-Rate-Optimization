function [re] = computeSTCCrate(Q,a,Ns,D,gamma)
re=0;
for d=1:D
    V(d)=0;
    C(d)=0;
    for i=1:Ns
        C(d)=C(d)+log2(1+Q(i,d)*gamma(i));
        V(d)=V(d)+1-1/(1+Q(i,d)*gamma(i))^2;
    end
    re=re+C(d)-a*sqrt(V(d));
end
clear V C

   