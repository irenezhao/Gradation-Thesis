
% programmed by Yue Zhao
% reference Dashuang Chen, T.Yamamura and H.Harada
%
function [ber] = ofdm_viterbi(ebn0)


ConvEncoder=comm.ConvolutionalEncoder(...
        'TerminationMethod','Terminated');
    Viterbi=comm.ViterbiDecoder('InputFormat','Hard',...
        'TerminationMethod','Terminated');

%********************** preparation part ***************************

para=128;   % Number of parallel channel to transmit (points)
fftlen=512; % FFT length
noc=128;    % Number of carrier
nd=6;       % Number of information OFDM symbol for one loop
ml=2;       % Modulation level : QPSK
sr=250000;  % Symbol rate
br=sr.*ml;  % Bit rate per carrier
gilen=32;   % Length of guard interval (points)
%ebn0=3;     % Eb/N0 3 pre

miu_L=3.56995; % logistic映射 分型参数，3.56994...<=miu<4
init = 0.4; %losgistic混沌初值
a = 0.5;     %三角帐篷映射

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
%log = tent_seq(para,nd,ml,a,init);
seldata_L = mod(seldata+log,2);%模2相加


%************************* 卷积码

seldata_L = seldata_L';

seldata_encoded = ConvEncoder.step(seldata_L);
seldata_encoded = seldata_encoded';



%****************** Serial to parallel conversion ***********************

para2 = length(seldata_encoded)/(nd*ml);
%paradata=reshape(seldata,para,nd*ml); %  reshape : built in function
paradata=reshape(seldata_encoded,para2,nd*ml); %  reshape : built in function

%************************** QPSK modulation ***************************** 

[ich,qch]=qpskmod(paradata,para2,nd,ml);
kmod=1/sqrt(2); %  sqrt : built in function
ich1=ich.*kmod;
qch1=qch.*kmod;


%------------- data mapping (DC=0) -----------

[ich1,qch1]=crmapping_512(ich1,qch1,fftlen,nd);


%******************* IFFT ************************

x=ich1+qch1.*1i;
y=ifft(x);      %  ifft : built in function
ich2=real(y);   %  real : built in function
qch2=imag(y);   %  imag : built in function

%********* Guard interval insertion **********

[ich3,qch3]= giins(ich2,qch2,fftlen,gilen,nd);
fftlen2=fftlen+gilen;





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

%*************** demapping
[ich6,qch6]=crdemapping_512(ich6,qch6,fftlen,nd);

%***************** demoduration *******************

ich7=ich6./kmod;
qch7=qch6./kmod;
[demodata]=qpskdemod(ich7,qch7,para2,nd,ml);   

%绘制星座图
if 0
    plot(ich7,qch7,'g.');
    grid on;
    hold on;
end


%**************  Parallel to serial conversion  *****************

demodata1=reshape(demodata,1,para2*nd*ml);

%************** Viterbi decoding**********************************
demodata1 = demodata1';
demodata1_decoded = Viterbi.step(demodata1);
demodata1_decoded = demodata1_decoded(1:length(seldata));

demodata1_decoded = demodata1_decoded';

%************************** 混沌解密 ************************************
demodata1 = mod(demodata1_decoded+log,2);%模2相加




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

%fprintf('%f\t%e\t%e\t%d\t\n',ebn0,ber,per,nloop);
fid = fopen('BERofdm.dat','a');
fprintf(fid,'%f\t%e\t%e\t%d\t\n',ebn0,ber,per,nloop);
fclose(fid);

%******************** end of file ***************************
