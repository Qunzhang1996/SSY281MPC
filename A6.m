


%1a
A=[
1.0041 0.0100 0 0;
0.8281 1.0041 0 -0.0093;
0.0002 0.0000 1 0.0098;
0.0491 0.0002 0 0.9629]
B=[0.0007;0.1398;0.0028;0.5605]
%%A'SA-S+Q=0
S=eye(4);
Q=S-A'*S*A
eigvalue=eig(A)
eigvalue_Q=eig(Q)
eig(S-Q)
%%%eig_Q have to be positive,unstable
%%
%b
%x(k+1)=(A-BK)x(k)
K=[114.3879 12.7189 -1.2779 -1.5952]
Q2=S-(A-B*K)'*S*(A-B*K)
eigvalue_Q2=eig(Q2)
eig(A-B*K)
%%%eig_Q have to be positive,unstable
%%
%c
Q3=eye(4)
S3=dlyap((A-B*K)',Q3)
eigvalue_S3=eig(S3)
%%%postive, means S positive definite matrix
%%
%2a
A=[
1.0041 0.0100 0 0;
0.8281 1.0041 0 -0.0093;
0.0002 0.0000 1 0.0098;
0.0491 0.0002 0 0.9629];
B=[0.0007;0.1398;0.0028;0.5605];
Q_2=eye(4);R_2=1;Pf=Q_2
%no affect the result
%%
[K_dp,P0,eigvalue_dp,N] = dp_func(A,B,1000,Q_2,R_2,Pf)
%%
%2C
eig(P0)
Q_c=P0-(A+B*K_dp)'*P0*(A+B*K_dp)
eig(Q_c)
S=dlyap((A+B*K_dp),Q_c)
eigs=eig(S)
eig(A+B*K_dp)
%%
%2d
Pf=idare(A,B,Q_2,R_2)
eig(Pf)
[k_d,~,~,~,] = MPC_batch(A,B,1,Q_2,R_2,Pf)
Q_d=Pf-(A+B*k_d)'*Pf*(A+B*k_d)
eig(Q_d)
%%
%3
%%
%4
clc
clear
A=[2 0.1;0 1.5]
B=[0;1]
Q=eye(2);R=100
pf1=eye(2)
pf2=idare(A,B,Q,R)
x_init=[-0.1 1.3]'
model=LTISystem('A',A,'B',B)
model.u.min = -1;
model.u.max = 1;
model.x.penalty = QuadFunction(Q);
model.u.penalty = QuadFunction(R);
C=eye(2)

Xf_constraint =Polyhedron('A',[C;-C],'b',[0.5;0.5;0.5;0.5]); 
model.x.with('terminalSet');
model.x.terminalSet = Xf_constraint
X = Polyhedron('lb', model.x.min, 'ub', model.x.max);
U = Polyhedron('lb', model.u.min, 'ub', model.u.max);
%%
C_inf=model.invariantSet();

%%
Tset=C_inf.intersect(Xf_constraint)
model.x.terminalSet = Tset;
figure()
plot(Tset)
legend("xf")
XN = Tset;
N=2
for i = 1:N
    pre_XN = model.reachableSet('X', XN, 'U', U, 'direction', 'backward','N',1);
    XN = pre_XN.intersect(X);
end
X0=XN
figure()
plot(X0)
legend("X0")
hold off

%%
%p1
x0 = [-0.1;1.3];
Nsim = 10;


model.x.with('terminalPenalty');
model.x.terminalPenalty = QuadFunction(pf1);
mpc  = MPCController(model, 2);
loop = ClosedLoop(mpc, model);
datap1  = loop.simulate(x0, Nsim);

figure()

plot(0:6,datap1.X(1,:),'LineWidth',2)
hold on
plot(0:6,datap1.X(2,:),'LineWidth',2)
legend('X1','X2')
hold off
%%
figure()
plot(X0,'alpha',0.5,Tset,'alpha',0.3)
hold on
plot(datap1.X(1,:),datap1.X(2,:), 'Linewidth', 2)
scatter(x0(1), x0(2), 'o','LineWidth',4)
scatter(datap1.X(1,end), datap1.X(2,end), '*','LineWidth',4)
title('X_0, X_f and trajectory')
legend('X0','Xf','trajectory','x0','xf')
%%
%p2
x0 = [-0.1;1.3];
Nsim = 20;
model.x.with('terminalPenalty');
model.x.terminalPenalty = QuadFunction(pf2);
mpc  = MPCController(model, N);
loop = ClosedLoop(mpc, model);
data  = loop.simulate(x0, Nsim);

figure()

plot(0:Nsim,data.X(1,:),'LineWidth',2)
hold on
plot(0:Nsim,data.X(2,:),'LineWidth',2)
legend('X1','X2')
hold off
%%
figure()
plot(X0,'alpha',0.5,Tset ,'alpha',0.3)
hold on
plot(data.X(1,:),data.X(2,:), 'Linewidth', 2)
scatter(x0(1), x0(2), 'o','LineWidth',4)
scatter(data.X(1,end), data.X(2,end), '*','LineWidth',4)
title('X_0, X_f and trajectory')
legend('X0','Xf','trajectory','x0','xf')
% plot(data.X(2,end),'LineWidth',2)
% plot(X0,'alpha',0.5)
% hold on
% plot(Tset,'alpha',0.4)
% hold on
%%
function [K_dp,P0,eigvalue_dp,N] = dp_func(A2,B2,N_input,Q,R,Pf)
Pkk = Pf;
for N=1:N_input
    count1=0;
    K_dp= -inv(R+B2'*Pkk*B2)*B2'*Pkk*A2;
    Pk  = Q + A2'*Pkk*A2 - A2'*Pkk*B2*inv(R+B2'*Pkk*B2)*B2'*Pkk*A2;
    Pkk = Pk;
    P0 = Pk;
    eigvalue_dp=eig(A2+B2*K_dp);
    for i =1:size(A2,1)
        if eigvalue_dp(i)^2<1
            count1=count1+1;
        end 
    end
    if count1==size(A2,1)
        fprintf("This is the result of dp solution\n")
        fprintf("the stable value of dp is %d",N);
        break
    end
end
end
function [k,Omega,Gamma,Q_bar,R_bar] = MPC_batch(A2,B2,N_input,Q,R,Pf)
for N=1:N_input
    n=size(A2,1);%A is nxn, get n
    p=size(B2,2);%B is nxp, get p
    Q_bar = blkdiag(kron(eye(N-1),Q),Pf);
    R_bar = kron(eye(N),R);
    Gamma = zeros(n*N,p*N);
    Omega = zeros(n*N,n);    
    for i=0:N-1
        Gamma = Gamma + kron(diag(ones(N-i,1),-i),A2^i*B2);
        Omega = Omega + [zeros(i*n,n); A2^(i+1); zeros((N-1-i)*n,n)];
    end
    kb=-inv(Gamma'*Q_bar*Gamma+R_bar)*Gamma'*Q_bar*Omega;
    k=kb(1,:);
end
end