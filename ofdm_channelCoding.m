% Program 4-1
% ofdm.m
%
% Simulation program to realize OFDM transmission system
%
% programmed by T.Yamamura and H.Harada
%

%********************** preparation part ***************************

para=128;   % Number of parallel channel to transmit (points)
para_channelCoding = 128*2;
fftlen=128; % FFT length
fftlen_channelCoding = 128*2;
noc=128;    % Number of carrier
nd=6;       % Number of information OFDM symbol for one loop
ml=2;       % Modulation level : QPSK
sr=250000;  % Symbol rate
br=sr.*ml;  % Bit rate per carrier
gilen=32;   % Length of guard interval (points)
ebn0=3;     % Eb/N0 3 pre

miu_L=3.56995; % logistic映射 分型参数，3.56994...<=miu<4
init = 0.4; %混沌初值

%************************** main loop part **************************

nloop=100;  % Number of simulation loops

noe = 0;    % Number of error data
nod = 0;    % Number of transmitted data
eop=0;      % Number of error packet
nop=0;      % Number of transmitted packet

for iii=1:nloop

%************************** transmitter *********************************

%************************** Data generation **************************** 

seldata=rand(1,para*nd*ml)>0.5;  %  rand : built in function


%************************** 混沌加密 ************************************
log = logistic_seq(para,nd,ml,miu_L,init);
seldata_L = mod(seldata+log,2);%模2相加

%************************** 卷积编码器 *************************************

M = 4; codeRate = 1/2; constlen = 7; k = log2(M); codegen = [171 133];
tblen = 32;     % traceback length
trellis = poly2trellis(constlen, codegen);
[msg_enc_bi] = convenc(seldata_L, trellis);



%****************** Serial to parallel conversion ***********************

%paradata=reshape(seldata,para,nd*ml); %  reshape : built in function
%paradata=reshape(seldata_L,para,nd*ml); %  reshape : built in function

paradata=reshape(msg_enc_bi,para_channelCoding,nd*ml); %  reshape : built in function


%************************** QPSK modulation ***************************** 

[ich,qch]=qpskmod(paradata,para_channelCoding,nd,ml);
kmod=1/sqrt(2); %  sqrt : built in function
ich1=ich.*kmod;
qch1=qch.*kmod;

%************************* 相位旋转 **************************************

grid on;
hold on;

%******************* IFFT ************************

x=ich1+qch1.*1i;
y=ifft(x);      %  ifft : built in function
ich2=real(y);   %  real : built in function
qch2=imag(y);   %  imag : built in function

%********* Guard interval insertion **********

[ich3,qch3]= giins(ich2,qch2,fftlen_channelCoding,gilen,nd);
fftlen2=fftlen_channelCoding+gilen;

%********* Attenuation Calculation *********

spow=sum(ich3.^2+qch3.^2)/nd./para;  %  sum : built in function
attn=0.5*spow*sr/br*10.^(-ebn0/10);
attn=sqrt(attn);

%***************************  Receiver  *****************************
%***************** AWGN addition ********* 

[ich4,qch4]=comb(ich3,qch3,attn);

%****************** Guard interval removal *********

[ich5,qch5]= girem(ich4,qch4,fftlen2,gilen,nd);

%******************  FFT  ******************

rx=ich5+qch5.*1i;
ry=fft(rx);   	% fft : built in function
ich6=real(ry);	% real : built in function
qch6=imag(ry);	% imag : built in function


%***************** demoduration *******************

ich7=ich6./kmod;
qch7=qch6./kmod;
[demodata]=qpskdemod(ich7,qch7,para_channelCoding,nd,ml);   


plot(ich7,qch7,'g.');
grid on;
hold on;



%**************  Parallel to serial conversion  *****************

demodata1=reshape(demodata,1,para_channelCoding*nd*ml);


%************************ 维特比译码器 ********************************

demodata1 = vitdec(demodata1, trellis, tblen, 'cont', 'hard');

%************************** 混沌解密 ************************************
demodata1 = mod(demodata1+log,2);%模2相加




%************************** Bit Error Rate (BER) ****************************

% instantaneous number of error and data

noe2=sum(abs(demodata1-seldata));  %  sum : built in function
nod2=length(seldata);  %  length : built in function

% cumulative the number of error and data in noe and nod

noe=noe+noe2;
nod=nod+nod2;

% calculating PER

if noe2~=0  
   eop=eop+1;
else
   eop=eop;
end   
   eop;
   nop=nop+1;
   

%fprintf('%d\t%e\t%d\n',iii,noe2/nod2,eop);  %  fprintf : built in function
   
end

%********************** Output result ***************************

per=eop/nop;
ber=noe/nod;

fprintf('%f\t%e\t%e\t%d\t\n',ebn0,ber,per,nloop);
fid = fopen('BERofdm.dat','a');
fprintf(fid,'%f\t%e\t%e\t%d\t\n',ebn0,ber,per,nloop);
fclose(fid);

%******************** end of file ***************************
