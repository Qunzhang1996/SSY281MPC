A=[1.0041 0.0100 0 0;0.8281 1.0041 0 -0.0093;
    0.0002 0.0000 1 0.0098;0.0491 0.0002 0 0.9629];
B=[0.0007 0.01;0.1398 1;0.0028 0;0.5605 0];
C=[1 0 0 0;
    0 0 1 0];
ys=[pi/18 -pi]';

%%
%%%%%%%question 1a%%%%%%%%%%%
[ans1_x,ans1_u] = setpoint(A,B,C,ys)
%%
%%%%%%question 1b%%%%%%%%%%%%
B2=[0.01 1 0 0]';C2=[1 0 0 0;0 0 1 0];
[ans2_x,ans2_u] = setpoint(A,B2,C2,ys)
%%
%%%%%question 1c%%%%%%%%%%%%%
B3=[0.0007 0.01;0.1398 1;0.0028 0;0.5605 0];C3=[1 0 0 0];
ys1=pi/18;
[ans3_x,ans3_u] = setpoint(A,B,C3,ys1)
%%
%%%%%%question 2a%%%%%%%%%%%%
A4=[1.0041 0.0100 0 0;0.8281 1.0041 0 -0.0093;
    0.0002 0.0000 1 0.0098;0.0491 0.0002 0 0.9629];
B4=[0.0007 0.1398 0.0028 0.5605]';
C4=[1 0 0 0;0 0 1 0];
Bp=B4;
Cp=[0;1];
Bd_Cd={zeros(4,1),[1;1];zeros(4,2),eye(2);[zeros(4,1) Bp],eye(2)};
A_d=cell(3,1);
B_d=cell(3,1);
C_d=cell(3,1);

for i=1:3
    Bd=Bd_Cd{i,1}
    Cd=Bd_Cd{i,2}
    n=size(A4,1);
    p=size(Bd,2);
    A_d{i}=[A4 Bd;zeros(p,n),eye(p)];
    B_d{i}=[B4;zeros(p,size(B4,2))];
    C_d{i}=[C4,Cd];
    OBSER=obsv(A_d{i},C_d{i});
    unobsv = length(A_d{i}) - rank(OBSER)
    if unobsv ~= 0
        unob_A=A_d{i};
        unob_C=C_d{i};
        eigvalue=eig(A_d{i})
        if all(eigvalue<1) 
            fprintf("detectalbe")
        else
            fprintf("not detectable")
        end
    end
end
%%
%%%%%%question 2b%%%%%%%%%%%%
%model 1
i=1;Q=eye(size(A_d{i}));
R=eye(size(C_d{i},1));
[P1,lamda1]=idare(A_d{i}',C_d{i}',Q,R)
L11=P1*C_d{i}'*inv(C_d{i}*P1*C_d{i}'+R)
%model 3
i=3;Q=eye(size(A_d{i}));R=eye(size(C_d{i},1))
[P3,lamda3]=idare(A_d{i}',C_d{i}',Q,R)
L33=P3*C_d{i}'*inv(C_d{i}*P3*C_d{i}'+R)

%%
%%%%%%question 2c%%%%%%%%%%%%
Cd1=[1;1];Bd1=zeros(4,1);
Cd2=eye(2);Bd2=zeros(4,2);
Cd3=eye(2);Bd3=[zeros(4,1),Bp];
H=[0 1];
Mss1=inv([eye(4)-A4 -B4;H*C4 0])*[Bd1;-H*Cd1]
Mss2=inv([eye(4)-A4 -B4;H*C4 0])*[Bd2;-H*Cd2]
Mss3=inv([eye(4)-A4 -B4;H*C4 0])*[Bd3;-H*Cd3]
%%
%%%%%%question 2d%%%%%%%%%%%%
Q=diag([5 2 0.5 0.1]);
R=0.1;
m=1;N=50;
x=[pi/36 0 0 0 ]';
x_true=[pi/36 0 0 0 ]';
x_hat=[0 0 0 0 0]';
d_hat=0;
u=0;
X_pre=zeros(5,1000);%%observor
X_true=zeros(4,1000);
y=zeros(2,1000);%
d_hat1=zeros(1,1000);

% xe_h(:,1) = [zeros(n,1) ; zeros(nd,1)];
% x(:,1)  = x0;
for tf=1:1000
    if tf<=50
        X_pre(:,tf)=x_hat;
        X_true(:,tf)=x_true;
        d_hat1(:,tf)=d_hat;
        pk=0; 
        y(:,tf)=[1 0 0 0;0 0 1 0]*x_true+[1;1]*pk;
        x0=[x_true;0];
        [x,x_hat,x_true,u,d_hat] = MPC_kalman(A4,B4,C4,Q,R,Bd_Cd,x0,N,m,pk,x_hat,u)
        
    else
        X_pre(:,tf)=x_hat;
        X_true(:,tf)=x_true;
        d_hat1(:,tf)=d_hat;
        pk=0.2;
        y(:,tf)=[1 0 0 0;0 0 1 0]*x_true+[1;1]*pk
        x0=[x_true;pk];
        [x,x_hat,x_true,u,d_hat]  = MPC_kalman(A4,B4,C4,Q,R,Bd_Cd,x0,N,m,pk,x_hat,u);
    end
end
figure(1)
plot(0:999,X_pre(3,:));hold on;plot(0:999,X_true(3,:))
legend("Xpre","Xtrue")
hold off
figure(2)
plot(0:999,y(1,:));hold on;plot(0:999,y(2,:))
legend("y1","y2")
hold off
figure(3)
plot(0:999,d_hat(1,:));
legend("d_hat")
hold off
%%

Q=diag([5 2 0.5 0.1]);
Pf=Q;R=0.1;
m=3;N=50;
x=[pi/36 0 0 0 ]'
x_true=x;
x_hat=[0 0 0 0 0 0]'
u=0;
X_pre1=zeros(6,1000);%%observor
X_true1=zeros(4,1000);%%observor
y11=zeros(2,1000);%
d_hat=[0;0];
d_hat2=zeros(2,1000);
for tf=1:1000
    if tf<=50
        X_pre1(:,tf)=x_hat;
        X_true1(:,tf)=x_true;
        d_hat2(:,tf)=d_hat;
        pk=[0;0];
        y11(:,tf)=[1 0 0 0;0 0 1 0]*x_true+eye(2)*pk;
        x0=[x_true;pk];
        [x,x_hat,x_true,u,d_hat] = MPC_kalman(A4,B4,C4,Q,R,Bd_Cd,x0,N,m,pk,x_hat,u);
       
    else
        X_pre1(:,tf)=x_hat;
        X_true1(:,tf)=x_true;
        d_hat2(:,tf)=d_hat;
        pk=[0.2;0.2];
        y11(:,tf)=[1 0 0 0;0 0 1 0]*x_true+eye(2)*pk;
        x0=[x_true;pk];
        [x,x_hat,x_true,u,d_hat] = MPC_kalman(A4,B4,C4,Q,R,Bd_Cd,x0,N,m,pk,x_hat,u);
    end
end
figure(4)
plot(0:999,X_pre1(3,:));hold on;plot(0:999,X_true1(3,:))
hold on
legend("Xpre","Xtrue")
hold off

figure(5)
plot(0:999,y11(1,:));hold on;plot(0:999,y11(2,:))
legend("y1","y2")
hold off

figure(6)
plot(0:999,d_hat2(1,:));
legend("d_hat")