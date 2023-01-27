function [i,P_next] = dp_iteration(aim_error,init_error,P_last,A2,B2,Q,R)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
i=1;
while (init_error>=aim_error)
    P_next=Q+A2'*P_last*A2-A2'*P_last*B2*inv(R+B2'*P_last*B2)*B2'*P_last*A2;
    init_error=norm(P_next-P_last);
    P_last=P_next;
    i=i+1;
end
fprintf("This is the result of dp_literation solution\n")
fprintf("the  value of literation is %d",i);
end