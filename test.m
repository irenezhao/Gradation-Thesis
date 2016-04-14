% 2011-3-29 从matlab demos 当中找来的关于卷积编码和维特比译码的实例
figure(1) %
subplot(211) %
EbNo = 4.5:.5:7; linEbNo = 10.^(EbNo(:).*0.1);
M = 4; codeRate = 1/2; constlen = 7; k = log2(M); codegen = [171 133];
tblen = 32;     % traceback length
trellis = poly2trellis(constlen, codegen);
dspec = distspec(trellis, 7);
expVitBER = bercoding(EbNo, 'conv', 'hard', codeRate, dspec);
semilogy(EbNo, expVitBER, 'g');  xlabel('Eb/No (dB)'); ylabel('BER');
title('Performance for R=1/2, K=7 Conv. Code and QPSK with Hard Decision');
grid on; axis([4 8 10e-7 10e-3]); legend('Union Bound', 0);
figure(2) %
subplot(311) %
numSymb = 100; numPlot = 20;
Nsamp = 4;      % oversampling rate
EbNoDemo = 3; EsN0 = EbNoDemo + 10*log10(k);
seed = [654321 123456];
rand('state', seed(1));  randn('state', seed(2));
msg_orig = randsrc(numSymb, 1, 0:1);
stem(0:numPlot-1, msg_orig(1:numPlot),'bx');
axis([ 0 numPlot -0.2 1.2]);  xlabel('时间'); ylabel('幅度');
title('卷积编码前的二进制符号' );
legend off
subplot(312) %
[msg_enc_bi] = convenc(msg_orig, trellis);
numEncPlot = numPlot / codeRate; tEnc = (0:numEncPlot-1) * codeRate;
stem(tEnc, msg_enc_bi(1:length(tEnc)),'rx');
axis([min(tEnc) max(tEnc) -0.2 1.2]);  xlabel('时间'); ylabel('幅度');
title('卷积编码后的二进制符号' );
figure(3) %
subplot(211) %
randn('state', seed(2));
msg_enc = bi2de(reshape(msg_enc_bi, ...
   size(msg_enc_bi,2)*k,size(msg_enc_bi,1) / k)');
grayencod = bitxor(0:M-1, floor((0:M-1)/2));
msg_gr_enc = grayencod(msg_enc+1);
msg_tx = modulate(modem.pskmod(M, pi/4), msg_gr_enc);
msg_tx = rectpulse(msg_tx, Nsamp);
msg_rx = awgn(msg_tx, EsN0-10*log10(1/codeRate)-10*log10(Nsamp));
numModPlot = numEncPlot * Nsamp ./ k;
tMod = (0:numModPlot-1) ./ Nsamp .* k;
plot(tMod, real(msg_tx(1:length(tMod))),'c-', ...
   tMod, imag(msg_tx(1:length(tMod))),'m-');
axis([ min(tMod) max(tMod) -1.5 1.5]);  xlabel('时间'); ylabel('幅度');
title('QPSK基带调制后的编码符号');
legend('同相分量', '正交分量', 0);
subplot(212) %
msg_rx_int = intdump(msg_rx, Nsamp);
msg_gr_demod = demodulate(modem.pskdemod(M, pi/4), msg_rx_int);
[dummy graydecod] = sort(grayencod); graydecod = graydecod - 1;
msg_demod = graydecod(msg_gr_demod+1)';
msg_demod_bi = de2bi(msg_demod,k)'; msg_demod_bi = msg_demod_bi(:);
stem(tEnc, msg_enc_bi(1:numEncPlot),'rx'); hold on;
stem(tEnc, msg_demod_bi(1:numEncPlot),'bo'); hold off;
axis([0 numPlot -0.2 1.2]);
xlabel('时间'); ylabel('幅度'); title('解调符号' );
figure(2) %
subplot(313) %
msg_dec = vitdec(msg_demod_bi, trellis, tblen, 'cont', 'hard');
stem(0:numPlot-1, msg_orig(1:numPlot), 'rx'); hold on;
stem(0:numPlot-1, msg_dec(1+tblen:numPlot+tblen), 'bo'); hold off;
axis([0 numPlot -0.2 1.2]);  xlabel('时间'); ylabel('幅度');
title('译码符号' );
% cla; %
figure(1) %
subplot(212) %
load('vitsimresults.mat');
semilogy(EbNo, expVitBER, 'g', EbNo, ratio, 'b*-');
xlabel('Eb/No (dB)'); ylabel('BER');
title('Performance for R=1/2, K=7 Conv. Code and QPSK with Hard Decision');
axis([4 8 10e-7 10e-3]); legend('Union Bound', 'Simulation Results', 0); grid on;