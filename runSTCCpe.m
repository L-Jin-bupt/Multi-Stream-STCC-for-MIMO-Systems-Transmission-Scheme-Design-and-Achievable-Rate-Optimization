function re = runSTCCpe(gamma,P,D,a,Ns,Q)
Q=Q(:,1:D);
for i=1:Ns
      G(i,:)=zeros(1,D);
      [~,idx]=max(Q(i,:));
      G(i,idx)=Q(i,idx);
end
G1 = G * (P / sum(G(:)));
S_ini=zeros(Ns,D);
for i=1:Ns
    for d=1:D
        if(G(i,d)>0)
            S_ini(i,d)=1;
        end
    end
end
for d=1:D
    for i=1:Ns
        if(Q(i,d)==0)
            Q(i,d)=10^-4;
        end
    end
end
Q_ini=Q;
ir=1;
re(ir)=computeSTCCrate(G1,a,Ns,D,gamma);
re2(ir)=computeSTCCrate(Q_ini,a,Ns,D,gamma);
fal=0;
beta=100;
while true %outer iteration
     ir=ir+1;
%      cvx_precision best
     cvx_solver sedumi
     cvx_begin 
            variable Q(Ns,D)
            variable S(Ns,D)
            variable p(Ns,1)
            for d=1:D
                gradfac(d)=0;
                temp=0;
                for i=1:Ns
                    temp=temp+1-1/(1+Q_ini(i,d)*gamma(i))^2;
                end
                gradfac(d)=a/sqrt(temp);
            end
            f1=0;f2=0;f3=0;
            for d=1:D
                for i=1:Ns
                    grad(i,d)=abs(gradfac(d)*(gamma(i))/(1+Q_ini(i,d)*gamma(i))^3);
                    f1=f1+a*grad(i,d)*Q(i,d);
                    f2=f2+log(1+Q(i,d)*gamma(i));
                    f3=f3+S(i,d)-2*S_ini(i,d)*S(i,d);
                end
            end
            f3=beta*f3;
            minimize f1-f2+f3
            subject to
            p>=0;
            for d=1:D
                for i=1:Ns
                    Q(i,d)<=P*S(i,d);
                    0<=Q(i,d)<=p(i);
                    Q(i,d)>=p(i)-(1-S(i,d))*P;
                end
            end
            sum(p)<=P;
            for i=1:Ns
                sum(S(i,:))<=1;
            end
            cvx_end
            if strcmp(cvx_status, 'Infeasible')
                fal==1;
                break;
            end
            if isnan(cvx_optval)
                fal==1;
                break;
            end
            S_ini=S;
            Q_ini=Q;
            for i=1:Ns
                [~,idx]=max(S(i,:));
                Sr(i,:)=zeros(1,D);
                Sr(i,idx(1))=S(i,idx(1))/S(i,idx(1));
            end
            re(ir)=abs(computeSTCCrate(p.*Sr,a,Ns,D,gamma));
            re2(ir)=abs(computeSTCCrate(Q,a,Ns,D,gamma));
            if(ir>100||abs(re2(ir)-re2(ir-1))<=0.001)
                break;
            end
end
end