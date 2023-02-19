
%%
%%%%%question 3%%%%%%%%%%%%%%%%%%%%
A=[0.4889 0.2939;
   1.0347 -0.7873;
   0.7269 0.8884;
   -0.3034 -1.1471];
B=[-1.0689;
   -0.8095;
   -2.9443;
   1.4384];
Ain_q3=[A,-ones(4,1);-A,-ones(4,1)];%%%F
Bin_q3=[B;-B];%G
f_q3=[zeros(2,1);1];
lb=[-inf(2,1);0];

[x,fval]= linprog(f_q3,Ain_q3,Bin_q3,[],[],lb,[])
%%%%%%dual question%%%%%%%%%%%%%%%%
g=[B;-B];
C=[zeros(2,1);1];
Aeq_q4=[A,-ones(4,1);-A,-ones(4,1)]';
beq_q4=-C;
lb=zeros(8,1);
% ub=zeros(8,1);
options = optimoptions('linprog','Algorithm','dual-simplex');
[x2,fval2,]= linprog(g,[],[],Aeq_q4,beq_q4,lb,[],options)
fval2=-fval2
%%%%%%%%%%%%%%%%%%%%%add salck condition
Z=inv([Ain_q3(3,:);Ain_q3(4,:);Ain_q3(6,:)])*[Bin_q3(3,:);Bin_q3(4,:);Bin_q3(6,:)]

%%
%%%%%question 4%%%%%%%%%%%%%%%%%%%
H=eye(4);
f=[];
F=[-1 0 0 0; 1 0 0 0;0 -1 0 0;0 1 0 0];
G=[0 0 -1 0;0 0 1 0;0 0 0 -1;0 0 0 1];
Ain=[F;G];
Bin=[-2.5 5 0.5 0.5 2 2 2 2]';
Aeq=[1 0 -1 0;0.4 -1 0 1];
x0=1.5;
Beq=[0.4*x0;0]
[x_star1,fval41,~,~,lambda1]=quadprog(H,f,Ain,Bin,Aeq,Beq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
gradf=H*x_star1;
gradg = [Aeq; Ain]';
KKT1 = gradf + gradg*[lambda1.eqlin;lambda1.ineqlin] %first KKT condition
KKT2 = lambda1.ineqlin%the second KKT condition
H=Aeq*x_star1-Beq
KKT3 = Ain*(x_star1)-Bin % the third KKT condition
KKT4 = lambda1.ineqlin.*(Ain*x_star1-Bin) %the fourth KKT condition
%%
H=eye(4);
f=[];
F=[0 0 0 0; 1 0 0 0;0 -1 0 0;0 1 0 0];
G=[0 0 -1 0;0 0 1 0;0 0 0 -1;0 0 0 1];
Ain=[F;G];
Bin=[0 5 0.5 0.5 2 2 2 2]';
Aeq=[1 0 -1 0;0.4 -1 0 1];
x0=1.5;
Beq=[0.4*x0;0];
[x_star2,fval42,~,~,lambda2]=quadprog(H,f,Ain,Bin,Aeq,Beq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
gradf=H*x_star2;
gradg = [Aeq; Ain]';
KKT12 = gradf + gradg*[lambda2.eqlin;lambda2.ineqlin] %first KKT condition
KKT22 = lambda2.ineqlin%the second KKT condition
KKT32 = Ain*(x_star2)-Bin % the third KKT condition
KKT42= lambda2.ineqlin.*(Ain*x_star2-Bin) %the fourth KKT condition
%%
H=eye(4);
f=[];
F=[-1 0 0 0; 0 0 0 0;0 -1 0 0;0 1 0 0];
G=[0 0 -1 0;0 0 1 0;0 0 0 -1;0 0 0 1];
Ain=[F;G];
Bin=[-2.5 0 0.5 0.5 2 2 2 2]';
Aeq=[1 0 -1 0;0.4 -1 0 1];
x0=1.5;
Beq=[0.4*x0;0];
[x_star3,fval43,~,~,lambda3]=quadprog(H,f,Ain,Bin,Aeq,Beq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
gradf=H*x_star3;
gradg = [Aeq; Ain]';
KKT13 = gradf + gradg*[lambda3.eqlin;lambda3.ineqlin] %first KKT condition
KKT23 = lambda3.ineqlin%the second KKT condition
KKT33 = Ain*(x_star3)-Bin % the third KKT condition
KKT43= lambda3.ineqlin.*(Ain*x_star3-Bin) %the fourth KKT condition