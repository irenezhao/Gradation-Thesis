% Program 4-2
% ofdm_fading.m
%
% Simulation program to realize OFDM transmission system
% (under one path fading)
% 
% One-path Rayleigh fading leigh(非标准，标准->书p184)
%
% programmed by T.Yamamura and H.Harada
%
function [ber] = ofdm_fading_chaos(ebn0)
%********************** preparation part ***************************

para=128;   % Number of parallel channel to transmit (points)
fftlen=128; % FFT length
noc=128;    % Number of carrier
nd=6;       % Number of information OFDM symbol for one loop
ml=2;       % Modulation level : QPSK
sr=250000;  % Symbol rate
br=sr.*ml;  % Bit rate per carrier
gilen=32;   % Length of guard interval (points)
%ebn0=10;    % Eb/N0

miu_L=3.56995; % logistic映射 分型参数，3.56994...<=miu<4
init = 0.4; %混沌初值

%******************* Fading initialization ********************
% If you use fading function "sefade", you can initialize all of parameters.
% Otherwise you can comment out the following initialization.
% The detailed explanation of all of valiables are mentioned in Program 2-8.

% Time resolution

tstp=1/sr/(fftlen+gilen); 

% Arrival time for each multipath normalized by tstp
% If you would like to simulate under one path fading model, you have only to set 
% direct wave.

itau = [0];

% Mean power for each multipath normalized by direct wave.
% If you would like to simulate under one path fading model, you have only to set 
% direct wave.
dlvl = [0];

% Number of waves to generate fading for each multipath.
% In normal case, more than six waves are needed to generate Rayleigh fading
n0=[6];

% Initial Phase of delayed wave
% In this simulation four-path Rayleigh fading are considered.
th1=[0.0];

% Number of fading counter to skip 
itnd0=nd*(fftlen+gilen)*10;

% Initial value of fading counter
% In this simulation one-path Rayleigh fading are considered.
% Therefore one fading counter are needed.
  
itnd1=[1000];

% Number of directwave + Number of delayed wave
% In this simulation one-path Rayleigh fading are considered
now1=1;        

% Maximum Doppler frequency [Hz]
% You can insert your favorite value
fd=320;       

% You can decide two mode to simulate fading by changing the variable flat
% flat     : flat fading or not 
% (1->flat (only amplitude is fluctuated),0->nomal(phase and amplitude are fluctutated)
flat =1;

%************************** main loop part **************************

nloop=500;  % Number of simulation loops

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

%****************** Serial to parallel conversion ***********************

paradata=reshape(seldata_L,para,nd*ml); %  reshape : built in function

%************************** QPSK modulation ***************************** 

[ich,qch]=qpskmod(paradata,para,nd,ml);
kmod=1/sqrt(2); %  sqrt : built in function
ich1=ich.*kmod;
qch1=qch.*kmod;

%******************* IFFT ************************

x=ich1+qch1.*i;
y=ifft(x);      %  ifft : built in function
ich2=real(y);   %  real : built in function
qch2=imag(y);   %  imag : built in function

%********* Gurad interval insertion **********

[ich3,qch3]= giins(ich2,qch2,fftlen,gilen,nd);
fftlen2=fftlen+gilen;

%********* Attenuation Calculation *********

spow=sum(ich3.^2+qch3.^2)/nd./para;  %  sum : built in function
attn=0.5*spow*sr/br*10.^(-ebn0/10);
attn=sqrt(attn);

%********************** Fading channel **********************

% Generated data are fed into a fading simulator
[ifade,qfade]=sefade(ich3,qch3,itau,dlvl,th1,n0,itnd1,now1,length(ich3),tstp,fd,flat);
  
% Updata fading counter
itnd1 = itnd1+ itnd0;

%***************************  Receiver  *****************************
%***************** AWGN addition ********* 

[ich4,qch4]=comb(ifade,qfade,attn);

%****************** Guard interval removal *********

[ich5,qch5]= girem(ich4,qch4,fftlen2,gilen,nd);

%******************  FFT  ******************

rx=ich5+qch5.*i;
ry=fft(rx);   	% fft : built in function
ich6=real(ry);	% real : built in function
qch6=imag(ry);	% imag : built in function

%***************** demoduration *******************

ich7=ich6./kmod;
qch7=qch6./kmod;
[demodata]=qpskdemod(ich7,qch7,para,nd,ml);   

%**************  Parallel to serial conversion  *****************

demodata1=reshape(demodata,1,para*nd*ml);

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

%fprintf('%f\t%e\t%e\t%d\t\n',ebn0,ber,per,nloop);
fid = fopen('BERofdmfad.dat','a');
fprintf(fid,'%f\t%e\t%e\t%d\t\n',ebn0,ber,per,nloop);
fclose(fid);

%******************** end of file ***************************
