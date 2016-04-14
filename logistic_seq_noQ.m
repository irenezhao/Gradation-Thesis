function [log]=logistic_seq_noQ(m,n,miu_L,init)
log = zeros(1,m*n);
log(1) = init;
for i=2:m*n 
    log(i) = miu_L*log(i-1)*(1-log(i-1));
end