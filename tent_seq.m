% Function_new 1 
% tent_seq.m
%
% Function to generate a sequence of logistic chaotic coding
%
% Programmed by Yue Zhao
%
function [tent]=tent_seq(para,nd,ml,a,init)
tent = zeros(1,para*nd*ml);
tent(1) = init;
for i=2:para*nd*ml 
    if tent(i-1)<0.5&&tent(i-1)>0
        tent(i) = 2*tent(i-1);
    else
        tent(i) = 2-2*tent(i-1);
    end
end
for i=1:para*nd*ml
    if tent(i)>=0&&tent(i)<0.5
        tent(i)=0;
    else tent(i)=1;
    end
end