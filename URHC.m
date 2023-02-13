function [Z,VN]=URHC(A,B,N,Q,R,Pf,x0,n)
%% This function solves a quadratic program (QP) to compute the optimal control 
% of an unconstrained finite-horizon LQ problem, which is in uncondensed form,
% i.e. both x and u are optimization variables.
%% A and B are the system matrices when x(k+1)=Ax(k)+Bu(k)
%% Q, R, and Pf are the gains in the cost function
%% N is the length of the horizon
%% Z is the vector of optimal variables and VN is the cost function 
%% x0 is the initial condition
% Objective function
n=size(A,1);
p=size(B,2);
Qbar = blkdiag(kron(eye(N-1),Q),Pf);
Rbar = kron(eye(40),R);
H    = blkdiag(Qbar,Rbar);
f    = [];
% Equality constraints
I    = eye(n);
Aeq1 = kron(eye(N),I)+kron(diag(ones(N-1,1),-1),-A);
Aeqb1= kron(eye(40),-B);
Aeqb2=zeros(1*n,40*p);
Aeqb2(:,end)=-B;
Aeq2=[Aeqb1;repmat(Aeqb2,10,1)];
Aeq  = [Aeq1 Aeq2];
beq  = [A*x0;zeros(n*(N-1),1)];
% Aeq2 = kron(eye(N),-B);
% Aeq  = [Aeq1 Aeq2];
% beq  = [A*x0;zeros(n*(N-1),1)];
% (no) Inequality constraints
Ain  = [];
bin  = [];
% Solve QP
[Z,VN,exitflag,output,lambda] = quadprog(2*H,f,Ain,bin,Aeq,beq);
% Check if the problem has been solved to optimality
if exitflag ~= 1
    disp(['The exitflag is ' num2str(exitflag) ', which means that something is wrong.'])
    pause
end
% Note: output returns some more information on how it went for the solver
% Note: lambda returns the Lagrange multipliers. Keep that in mind, it will
% be useful later in the course !