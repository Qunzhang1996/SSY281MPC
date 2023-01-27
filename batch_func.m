function [k,eigvalue,N] = batch_func(A2,B2,N_input,Q,R,Pf)
for N=1:N_input
    n=size(A2,1);%A is nxn, get n
    p=size(B2,2);%B is nxp, get p
    omega=zeros(N*n,n); %omega is (N*n,n)
    gama=zeros((N+1)*n,N*p);%kaf is (N*n,N*p)
    tempt=eye(n);
    
    for i =1:N%add element by row
        row=i*n+(1:n);
        gama(row,:)=[tempt*B2,gama(row-n,1:end-p)];
        tempt=A2*tempt;
        omega(row,:)=tempt; 
    end
    gama(1:n,:)=[];
    omega(1:n,:)=[];
    S_q=size(Q,1);
    S_r=size(R,1);
    %initial Q_bar
    Q_bar=zeros(N*S_q,N*S_q);
    for i=0:N-1
        Q_bar(i*S_q+1:(i+1)*S_q,i*S_q+1:(i+1)*S_q)=Q;
    end
    Q_bar((N-1)*S_q+1:(N)*S_q,(N-1)*S_q+1:(N)*S_q)=Pf;
    %initial R_bar
    R_bar=zeros(N*S_r,N*S_r);
    for i=0:N-1
        R_bar(i*S_r+1:(i+1)*S_r,i*S_r+1:(i+1)*S_r)=R;
    end
    kb=-inv(gama'*Q_bar*gama+R_bar)*gama'*Q_bar*omega;
    k=kb(1,:);
    eigvalue=eig(A2+B2*k);
    count=0;
    %to get a stable eigvalue
    for i =1:size(A2,1)
        if eigvalue(i)^2<1
            count=count+1;
        end    
    end
    if count==size(A2,1)
        fprintf("This is the result of batch solution\n")
        fprintf("the stable value of batch is %d",N)
        break
    end
end
end