function [ans_x,ans_u] = setpoint(A,B,C,ys)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    n=size(A,1);
    m=size(B,2);
    p=size(C,1);
    if m==p
        fprintf("same output and same input!")
        ans111=inv([eye(4)-A -B;C zeros(p,m)])*[zeros(4,1);ys];
        ans_x=ans111(1:n,:);
        ans_u=ans111(n+1:end,:);
    else if m<p
        fprintf("less inputs than outputs!")
        %%%%objective function%%%%%
        %quad=0.5XHX+fX
        %question=Cz^2-2ys*Cz
        Cz=[C zeros(p,m)];
        H=Cz'*Cz;
        f=-2*ys'*Cz;
        %%%%Equality constraints
        Aeq=[eye(n)-A -B];
        beq=zeros(n,1);
        ans222=quadprog(2*H,f,[],[],Aeq,beq);
        ans_x=ans222(1:n,:);
        ans_u=ans222(n+1:end,:);
    else 
        fprintf("more inputs than outputs!")
        %question=Cz^2-2ys*Cz
        Cz=[C zeros(p,m)];
%         H=Cz'*Cz;
        H=blkdiag((C'*C),eye(2));
        f=-2*ys'*Cz;
        %%%%Equality constraints
        Aeq=[eye(n)-A -B;Cz];
        beq=[zeros(n,1);ys];
        ans222=quadprog(2*H,f,[],[],Aeq,beq);
        ans_x=ans222(1:n,:);
        ans_u=ans222(n+1:end,:);
    end

end