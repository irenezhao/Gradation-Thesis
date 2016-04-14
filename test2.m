
M = 4; codeRate = 1/2; constlen = 7; k = log2(M); codegen = [171 133];
tblen = 32;     % traceback length
trellis = poly2trellis(constlen, codegen);

%初始化
numSymb = 100;
msg_orig = randsrc(numSymb, 1, 0:1);

%卷积码
[msg_enc_bi] = convenc(msg_orig, trellis);

%维特比译码
msg_dec = vitdec(msg_enc_bi, trellis, tblen, 'cont', 'hard');
%msg_dec = vitdec(msg_demod_bi, trellis, tblen, 'cont', 'hard');

k =0;
for i=1:100;
    if (msg_dec(i)==msg_orig(i))
        k = k+1;
    end
end
        

