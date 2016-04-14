% Function_new 1 
% logistic_seq.m
%
% Function to generate a sequence of logistic chaotic coding
%
% Programmed by Yue Zhao
%
function [log]=logistic_seq(para,nd,ml,miu_L,init)
log = zeros(1,para*nd*ml);
log(1) = init;
for i=2:para*nd*ml 
    log(i) = miu_L*log(i-1)*(1-log(i-1));
end
for i=1:para*nd*ml
    if log(i)>=0&&log(i)<0.5
        log(i)=0;
    else log(i)=1;
    end
end

%{
count=0;
for i=1:1536
    if A(i)==1
        count = count+1;
    end
end
%}



