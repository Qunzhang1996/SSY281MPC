function [X_K,uk,X1_k,X2_k,X3_k,UK] = xu_plot(A2,B2,N_input,Q,R,Pf,X0,X1_steps,X2_steps,X3_steps,U_step)
    n=size(A2,1);
    p=size(B2,2);
    X_K=zeros(n,X3_steps);
    X_K(:,1)=X0;
    uk=zeros(p,X3_steps);
    %     u_k=U_kb(1:p,1);
    for k=1:X3_steps
        [kb,~] = MPC_batch(A2,B2,N_input,Q,R,Pf);
%         X_K(:,k)
%         kb
        U_kb=kb*X_K(:,k);
        uk(p,k)=U_kb(1,:);
        X_K(:,k+1)=A2*X_K(:,k)+B2*uk(p,k);
%         ,X1_steps,X2_steps,U_steps
    end
    X1_k=X_K(1,1:X1_steps);
    X2_k=X_K(2,1:X2_steps);
    X3_k=X_K(3,1:X3_steps);
    UK=uk(:,1:U_step);
end
