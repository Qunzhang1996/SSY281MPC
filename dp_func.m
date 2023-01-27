function [K_dp,P0,eigvalue_dp,N] = dp_func(A2,B2,N_input,Q,R,Pf)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for N=1:N_input
    count1=0;
    P=cell(N+1,1);
    P{N+1}=Pf;
    for i=N+1:-1:2
        P{i-1}=Q+A2'*P{i}*A2-A2'*P{i}*B2*inv(R+B2'*P{i}*B2)*B2'*P{i}*A2;
    end
    K_dp=-inv(R+B2'*P{2}*B2)*B2'*P{2}*A2;
    P0=P{1};
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