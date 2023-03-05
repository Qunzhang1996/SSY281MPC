clc
clear all
%Linear MPC design
%(a)
A=[1.2 1;0 1];
B=[0; 1];
model = LTISystem('A', [1.2 1;0 1], 'B', [0; 1]);% stable A
model.x.min = [-15; -15];
model.x.max = [ 15;  15];
model.u.min = -1;
model.u.max =  1;
X = Polyhedron('lb', model.x.min, 'ub', model.x.max);
U = Polyhedron('lb', model.u.min, 'ub', model.u.max);
Q=eye(2);
R=100;
N=4;
model.x.penalty = QuadFunction(Q);
model.u.penalty = QuadFunction(R);
%%%%%%%%%%%%%%%%
%Theorom 11.13%
%%%%%%%%%%%%%%%%%%%%%%%%K from LQR
Tset   = Polyhedron([0 0]);
[K,S,e] = dlqr(A,B,Q,R)
%%%%%%%%%%%%%Pf from dlyap
[K,~,~,~,~] = MPC_batch(A,B,N,Q,R,S)
X111 = dlyap((A-B*K)',(Q+K'*R*K),[],eye(size(A)))
%%%%%%%%%%%%%%%%%%
P = model.LQRPenalty
P.weight

model.x.with('terminalPenalty');
model.x.terminalPenalty = P;

model.x.with('terminalSet');
model.x.terminalSet = Tset;


% Compute the feasible set XN
XN = Tset;
for i = 1:N
    XN = model.reachableSet('X', XN, 'U', U, 'direction', 'backward');
    XN = XN.intersect(X).minHRep();
end
plot(XN)
%%
%(b)
N=[10,15,20];
x0=[7;-4];
for i=1:3
    pred=MPCController(model, N(i));
    ctrl = MPCController(model,N(i));
    [~, ~, data1] = pred.evaluate(x0);
    loop = ClosedLoop(ctrl, model);
    Nsim = 30;
    data2 = loop.simulate(x0, Nsim);
    figure()
    plot(0:N(i),data1.X(1,:))
    hold on
    plot(0:N(i),data1.X(2,:))
    hold on
    plot(0:Nsim,data2.X(1,:))
    hold on
    plot(0:Nsim,data2.X(2,:))
    legend('pred_x1','pred_x2','cl_x1','cl_x2')
    hold off
end
%%
%(c)
N=20
Tset=Polyhedron('lb', [-0.01;-0.01], 'ub', [0.01;0.01]);
model.x.terminalSet = Tset;
model.x.terminalPenalty = model.LQRPenalty;
mpc = MPCController(model, N);
expmpc = mpc.toExplicit();
[~, ~, data3] = expmpc.evaluate(x0);
figure()
expmpc.partition.plot()
%%
%(d)
C_inf=model.invariantSet();
Xf = model.reachableSet('X',C_inf,'U',U,'N',1,'direction','forward')
Xf = Xf.intersect(C_inf);

x01=model.reachableSet('X',Xf,'U',U,'N',1,'direction','backward')
x01= x01.intersect(X).minHRep();

figure()
plot(X,'alpha',0.5,'color','blue',C_inf,'alpha',0.8,'color','yellow',Xf)
legend({'X','C_inf','Xf'})

figure()
plot(X,'alpha',0.5,x01,C_inf,'alpha',0.8,'color','yellow','alpha',0.5)
legend({'X','X0','C_inf'})
%%
%2.Finite time control of a DC motor
%%
%(a)
Ls=1;ds=0.02;jm=0.5;bm=0.1;R=20;kt=10;
p=20;k0=1280.2;jl=50*jm;bl=25;
A=[0 1 0 0;-k0/jl -bl/jl k0/(p*jl) 0;0 0 0 1;k0/(p*jm) 0 -k0/(p^2*jm) -(bm+kt^2/R)/jm];
B=[0 0 0 kt/(R*jm)]';
C=[k0 0 -k0/p 0];
D=[0];
sys=ss(A,B,C,D);
h=0.1;
sys_d=c2d(sys,h);
Ad=sys_d.A
Bd=sys_d.B

%%
%(b)
model = LTISystem('A', sys_d.A, 'B', sys_d.B,'C',C);
Tset = Polyhedron('Ae',[0 1 0 0; 0 0 0 1],'be',[0; 0]);
model.u.min = -200;
model.u.max = 200;
u_c=Polyhedron('lb',model.u.min,'ub',model.u.max);
model.x.with('terminalSet');
model.x.terminalSet = Tset;
xb=[];
ub=[];
yb=[];
x0=[0 2.5 0 75]'
xN=Tset
k(1)=Tset
%%
N=0
% Tset = Polyhedron('Ae',[0 1 0 0; 0 0 0 1],'be',[0; 0]);
while ~Tset.contains(x0)
    N=N+1
    x_pre=model.reachableSet('X',Tset,'U',u_c,'N',1,'direction','backward');
    Tset=x_pre%no constrain here
    k(N+1)=Tset
end
%%
x0=[0 2.5 0 75]'
for i=1:N
    model.x.terminalSet = k(N-i+1);
    [x0,x,u,y] =MPC_Controller(Ad,Bd,model,1,x0);
    disp(i)
    x0
    xb(:,i)=x0
    ub(:,i)=u
end
xb1=[[0 2.5 0 75]',xb]
figure()
plot(0:3,xb1(1,:))
hold on 
plot(0:3,xb1(2,:))
hold on 
plot(0:3,xb1(3,:))
hold on 
plot(0:3,xb1(4,:))
hold on 
legend('x1','x2','x3','x4')
%%
y=C*xb1
figure()
plot(0:3,y)
%%
% %(C)
model = LTISystem('A', sys_d.A, 'B', sys_d.B,'C',C);
Tset = Polyhedron('Ae',[0 1 0 0; 0 0 0 1],'be',[0; 0]);
model.u.min = -200;
model.u.max = 200;
u_c=Polyhedron('lb',model.u.min,'ub',model.u.max);
model.x.with('terminalSet');
model.x.terminalSet = Tset;
X=Polyhedron('A',[C;-C],'b',[150;150]);  
model.x.with('setConstraint');
model.x.setConstraint=X;

% xb=[];
% ub=[];
% yb=[];
x0=[0 2.5 0 75]'
k(1)=Tset
%%
N=0
Tset = Polyhedron('Ae',[0 1 0 0; 0 0 0 1],'be',[0; 0]);
while ~Tset.contains(x0)
    N=N+1;
    x_pre=model.reachableSet('X',Tset,'U',u_c,'N',1,'direction','backward');
    Tset=x_pre.intersect(X);
    k(N+1)=Tset;
end


%%
x0=[0 2.5 0 75]'
for i=1:N
    model.x.terminalSet = k(N-i+1);
    [x0,x,u,y] =MPC_Controller(Ad,Bd,model,1,x0);
%     disp(i)
%     x0
    xc(:,i)=x0
    uc(:,i)=u
end
xc1=[[0 2.5 0 75]',xc]
figure()
plot(0:N,xc1(1,:))
hold on 
plot(0:N,xc1(2,:))
hold on 
plot(0:N,xc1(3,:))
hold on 
plot(0:N,xc1(4,:))
hold on 
legend('x1','x2','x3','x4')
%%
y=C*xc1
figure()
plot(0:N,y)
%%
%d
model = LTISystem('A', sys_d.A, 'B', sys_d.B,'C',C);
model.u.min = -200;
model.u.max = 200;
u_c=Polyhedron('lb',model.u.min,'ub',model.u.max);
model.x.with('terminalSet');
Tset = Polyhedron('lb',[-10;-0.01;-10;-0.01],'ub',[10;0.01;10;0.01]);
model.x.terminalSet = Tset;
X=Polyhedron('A',[C;-C],'b',[150;150]);  
model.x.with('setConstraint');
model.x.setConstraint=X;


xd=[];
ud=[];
yd=[];
x0=[0 2.5 0 75]'
xN=Tset
k(1)=Tset
%%
N=0
while ~Tset.contains(x0)
    N=N+1
    x_pre=model.reachableSet('X',Tset,'U',u_c,'N',1,'direction','backward');
    Tset=x_pre.intersect(X)
    k(N+1)=Tset
end
%%%%%%%%%%%%
%%
N=11
x0=[0 2.5 0 75]'
for i=1:N
    model.x.terminalSet = k(N-i+1);
    [x0,x,u,y] =MPC_Controller(Ad,Bd,model,1,x0);
%     disp(i)
%     x0
    xd(:,i)=x0;
    ud(:,i)=u;
end
xd1=[[0 2.5 0 75]',xd]
xd1(:,12)=[0.3 -0.010 8.33 -0.0100]'
figure()
plot(0:N,xd1(1,:))
hold on 
plot(0:N,xd1(2,:))
hold on 
plot(0:N,xd1(3,:))
hold on 
plot(0:N,xd1(4,:))
hold on 
legend('x1','x2','x3','x4')
%%
y=C*xd1
figure()
plot(0:N,y)
%%
x0=[0 2.5 0 75]'
N=11
i=1
Nsim=20
for m=20:-1:10
    mpc  = MPCController(model, N);
    loop = ClosedLoop(mpc, model);
    N=N-1
    edata = loop.simulate(x0, 1);
    x0=edata.X(:,2);
    xd2(:,i)=x0;
    i=i+1;
end

for i=12:20
    mpc  = MPCController(model, 1);
    loop = ClosedLoop(mpc, model);
    edata = loop.simulate(x0, 1);
    x0=edata.X(:,1);
    xd2(:,i)=x0;
end
xd2=[[0 2.5 0 75]',xd2]
plot(0:Nsim,xd2(1,:), 'Linewidth', 2)
hold on
plot(0:Nsim,xd2(2,:), 'Linewidth', 2)
plot(0:Nsim,xd2(3,:), 'Linewidth', 2)
plot(0:Nsim,xd2(4,:), 'Linewidth', 2)
% stairs(edata.U, 'k--')
legend({'x_1','x_2','x_3','x_4'})
Y=C*xd2;
figure()
plot(0:Nsim,Y, 'Linewidth', 2)
legend({'y'})
%%
function [x0,x,u,y] =MPC_Controller(Ad,Bd,model,N,x0)    
    mpc = MPCController(model,N);
    [~, ~, data4] = mpc.evaluate(x0);
    x=data4.X(:,1);
    u=data4.U(:,1);
    y=data4.Y(:,1);
    x0=Ad*x+Bd*u;
end
function [kb,Omega,Gamma,Q_bar,R_bar] = MPC_batch(A2,B2,N,Q,R,Pf)
for N=1:N
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
    kb=inv(Gamma'*Q_bar*Gamma+R_bar)*Gamma'*Q_bar*Omega;
    kb=kb(1,:);
end
end