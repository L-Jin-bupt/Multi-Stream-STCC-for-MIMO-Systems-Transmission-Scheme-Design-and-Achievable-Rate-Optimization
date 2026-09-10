function re2 = runSTCCproposed(gamma,P,D,a,Ns,rho,etarho,Q)
ir=1;
Q=Q(:,1:D);
re(ir)=computeSTCCrate(Q,a,Ns,D,gamma);
for i=1:Ns
      G(i,:)=zeros(1,D);
      [~,idx]=max(Q(i,:));
      G(i,idx)=Q(i,idx);
end
G1 = G * (P / sum(G(:)));
re2(ir)=computeSTCCrate(G1,a,Ns,D,gamma);
for d=1:D
    for i=1:Ns
        if(Q(i,d)==0)
            Q(i,d)=10^-8;
        end
    end
end
fal=0;%防止NaN的指示
while true %outer iteration
    ir=ir+1;
    m=1;
    reinner(m)=inf;
    while true %medium iteration 
      m=m+1;
      %% update G
      for i=1:Ns
          G(i,:)=zeros(1,D);
          [~,idx]=max(Q(i,:));
          G(i,idx)=Q(i,idx);
      end
      %% update Q
      k=1;
      reCCCP(m-1,k)=-computeSTCCrate(Q,a,Ns,D,gamma)+rho*norm(Q-G,'fro')^2;
      Q_ini=Q;
      while true %inner
          k=k+1;
%         cvx_precision best
        cvx_solver sedumi
        cvx_begin 
        variable Q(Ns,D)
        for d=1:D
            gradfac(d)=0;
            temp=0;
            for i=1:Ns
                temp=temp+1-1/(1+Q_ini(i,d)*gamma(i))^2;
            end
            gradfac(d)=1/sqrt(temp);
            if (isnan(gradfac(d)))
                fal=1;
                break
            end
        end
        if(fal==1)
            break;
        end
        f1=0;f2=0;
        for d=1:D
            for i=1:Ns
                grad(i,d)=gradfac(d)*(gamma(i))/(1+Q_ini(i,d)*gamma(i))^3;
                f1=f1+a*grad(i,d)*Q(i,d);
                f2=f2+log(1+Q(i,d)*gamma(i));
            end
        end
        f3=rho*square_pos(norm(Q-G,'fro'));
        minimize f1-f2+f3
        subject to
        sum(Q(:))<=P;
        Q>=0;
        cvx_end
        if strcmp(cvx_status, 'Solved')
             fal=0;
        else
            fal=1;
            break;
        end
        if (isnan(cvx_optval)==1)
            fal==1;
            break;
        end
      if(fal==1)
          break;
      end
        reCCCP(m-1,k)=-computeSTCCrate(Q,a,Ns,D,gamma)+rho*norm(Q-G,'fro')^2;
        Q_ini=Q;
        if(reCCCP(m-1,k-1)-reCCCP(m-1,k)<=10^-3)
            break;
        end
        clear grad gradfac
      end
      if(fal==1)
          break;
      end
      reinner(m)=reCCCP(m-1,k);
      if(reinner(m-1)-reinner(m)<=10^-3)
          break;
      end
    end
    if(fal==1)
        break;
    end
    G1 = G * (P / sum(G(:)));
    re(ir)=abs(computeSTCCrate(Q,a,Ns,D,gamma));
    re2(ir)=abs(computeSTCCrate(G1,a,Ns,D,gamma));
    rho=rho*etarho;
    if(fal==1)
        break;
    end
    if(abs(re2(ir)-re2(ir-1))<=10^-3&&abs(re2(ir)-re(ir))<=10^-1)
        break;
    end
    clear reCCCP reinner;
end
end