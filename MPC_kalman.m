function [x,x_hat,x_true,u,d_hat] = MPC_kalman(A4,B4,C4,Q,R,Bd_Cd,x0,N,m,pk,x_hat,u)
    Pf=Q;
    Bd=Bd_Cd{m,1};
    Cd=Bd_Cd{m,2};
    n=size(A4,1);
    p=size(Bd,2);
    A_1=[A4 Bd;zeros(p,n),eye(p)];%augmented A matrix
    B_1=[B4;zeros(p,size(B4,2))];%augmented B matrix
    C_1=[C4,Cd];%augmented C matrix
    
    [P1,~]=idare(A_1',C_1',1,1);%IDARE to solve the Raccati
    L=P1*C_1'*inv(C_1*P1*C_1'+eye(2));%kalman gain
    
    if size(A_1,1)<=5
        MSS=[0 0 -1 0 0]';%Mss matrix of case one
    else
        MSS=[0 0 0 0 0;0 0 -1 0 -1]';%Mss matrix of case three
    end

    d_hat=x_hat(5:end,:);
    Zs=MSS*d_hat;%set point of [xs;us]
    delta_x=x_hat(1:4,:)-Zs(1:4,:);

    [Z,~]=URHC(A4,B4,N,Q,R,Pf,delta_x,n);%Z is [delta_x;delta_u]
    
    u=Z((size(A4,1)*N+1),1)+Zs(5,:);% delta_u+us
%     x_hat=x_hat+L(:,2)*(C_1(2,:)*x0-C_1(2,:)*x_hat);%corrrection,x_hat=[x_hat;d_hat]
%     x_hat=A_1*x_hat+B_1*u;%update,x_hat(k+1)=A_1x_hat(k)+B_1u(last time)

    x_hat=x_hat+L*(C_1*x0-C_1*x_hat);%corrrection,x_hat=[x_hat;d_hat]
    x_hat=A_1*x_hat+B_1*u;%update,x_hat(k+1)=A_1x_hat(k)+B_1u(last time)
    
%     x_true=Z(1:4,:);
    x_true=A_1*x0+B_1*u;%Augmented A
    x_true=x_true(1:4,:);%[x1 x2 x3 x4]'
    
    

    x=x_hat(1:4,:);

end
