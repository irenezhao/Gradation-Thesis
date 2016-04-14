% Program 4-9
% ofdmci.m
%
% Simulation program to realize OFDM transmission system
% Simulate effect of interfarence noise
%
% Programmed by T.Yamamura and H.Harada
%
% GI CE GI data GI data...(data 6symbols)
% 

%********************** preparation part ***************************

para=52;      % Number of parallel channel to transmit (points)
fftlen=64;    % FFT length
noc=53;       % Number of carriers
nd=6;         % Number of information OFDM symbol for one loop
knd=1;        % Number of known channel estimation (CE) OFDM symbol
ml=2;         % Modulation level : QPSK
sr=250000;    % OFDM symbol rate (250 ksyombol/s)
br=sr.*ml;    % Bit rate per carrier
gilen=16;     % Length of guard interval (points)
ebn0=1000;    % Eb/N0

%---------------------- fading initialization ----------------------

tstp=1/sr/(fftlen+gilen); % Time resolution
itau=[0];       % Arrival time for each multipath normalized by tstp 
dlvl1=[0];      % Mean power for each multipath normalized by direct wave.
n0=[6];	        % Number of waves to generate fading n0(1),n0(2)
th1=[0.0];      % Initial Phase of delayed wave
itnd1=[1000];   % set fading counter        	
now1=1;         % Number of directwave + Number of delayed wave
fd=150;         % Maximum Doppler frequency
flat=0;         % Flat or not (see ofdm_fading.m)
itnd0=nd*(fftlen+gilen)*10; % Number of fading counter to skip 

%----------------- interference wave initialization --------------------

ci=10;           % C/I ratio 
ml2=2;           % modulation level

itau2=[0];
dlvl2=[0];
n02=[6];
th2=[0.0];
itnd2=[10000+floor(rand(1)*10)*1000];       	
now2=1;
fd2=fd;
flat2=0;
itnd02=nd*(fftlen+gilen)*300; % Number of fading counter to skip 

%% store all parameters in one matrix "fadingpara"

fadingpara=zeros(8,length(itau2));
fadingpara(1,:)=itau2;
fadingpara(2,:)=dlvl2;
fadingpara(3,:)=n02;
fadingpara(4,:)=th2;
fadingpara(5,:)=itnd2;
fadingpara(6,:)=now2;
fadingpara(7,:)=fd2;
fadingpara(8,:)=flat2;

%************************** main loop part **************************

nloop=1000;  % Number of simulation loops

noe = 0;    % Number of error data
nod = 0;    % Number of transmitted data
eop=0;      % Number of error packet
nop=0;      % Number of transmitted packet

%************************** transmitter *****************************
for iii=1:nloop
   
seldata=rand(1,para*nd*ml)>0.5;  %  DC=0

paradata=reshape(seldata,para,nd*ml); %size(51  *  nd*ml)

%-------------- ml modulation ---------------- 

[ich,qch]=qpskmod(paradata,para,nd,ml);
kmod=1/sqrt(2);
ich=ich.*kmod;
qch=qch.*kmod;

% CE data generation
kndata=zeros(1,fftlen);
kndata0=2.*(rand(1,52)>0.5)-1;
kndata(2:27)=kndata0(1:26);
kndata(39:64)=kndata0(27:52);
ceich=kndata; % CE:BPSK
ceqch=zeros(1,64);

%------------- data mapping (DC=0) -----------

[ich1,qch1]=crmapping(ich,qch,fftlen,nd);

ich2=[ceich.' ich1]; % I-channel transmission data
qch2=[ceqch.' qch1]; % Q-channel transmission data

%------------------- IFFT  -------------------

x=ich2+qch2.*i;
y=ifft(x);
ich3=real(y);
qch3=imag(y);

%---------- Gurad interval insertion ---------

fftlen2=fftlen+gilen;
[ich4,qch4]= giins(ich3,qch3,fftlen,gilen,nd+1);

%---------- Attenuation Calculation ----------
spow=sum(ich4.^2+qch4.^2)/nd./52;
attn=0.5*spow*sr/br*10.^(-ebn0/10);
attn=sqrt(attn);


%********************** fading channel ******************************   
%If you would like to simulate performance under fading, please remove "*"
%from the following four sentenses
[ifade,qfade,ramp,rcos,rsin]=sefade(ich4,qch4,itau,dlvl1,th1,n0,itnd1,now1,length(ich4),tstp,fd,flat);
itnd1 = itnd1+itnd0;  % Updata fading counter
ich4=ifade;
qch4=qfade;

%%% interference wave addition
% interference
[iintw,qintw]=interwave(ci,spow,ml2,length(ich4),tstp,fadingpara);
itnd2 = itnd2+itnd02;
fadingpara(5,:)=itnd2;
ich4=ich4+iintw;
qch4=qch4+qintw;

%***************************  Receiver  *****************************
%--------------- AWGN addition --------------- 
[ich5,qch5]=comb(ich4,qch4,attn);

%----Perfect fading compensation for one path fading ----
%If you would like to simulate performance under perfect compensation, please remove "*"
%from the following four sentenses
%ifade2=1./ramp.*(rcos(1,:).*ich5+rsin(1,:).*qch5);
%qfade2=1./ramp.*(-rsin(1,:).*ich5+rcos(1,:).*qch5);
%ich5=ifade2;
%qch5=qfade2;

%----------- Guard interval removal ----------

[ich6,qch6]= girem(ich5,qch5,fftlen2,gilen,nd+1);

%------------------  FFT  --------------------

rx=ich6+qch6.*i;
ry=fft(rx);
ich7=real(ry);
qch7=imag(ry);

%-------------- Fading compensation by CE symbol --------------
%
%If you would like to simulate performance under CE-based compensation, please remove "*"
%in this area
%

% preparation known CE data
ce=1;
ice0=ich2(:,ce);
qce0=qch2(:,ce);

% taking CE data out of received data
ice1=ich7(:,ce);
qce1=qch7(:,ce);

% calculating reverse rotation 
iv=real((1./(ice1.^2+qce1.^2)).*(ice0+i.*qce0).*(ice1-i.*qce1));
qv=imag((1./(ice1.^2+qce1.^2)).*(ice0+i.*qce0).*(ice1-i.*qce1));

% matrix for reverse rotation
ieqv1=[iv iv iv iv iv iv iv];
qeqv1=[qv qv qv qv qv qv qv];

% reverse rotation
icompen=real((ich7+i.*qch7).*(ieqv1+i.*qeqv1));
qcompen=imag((ich7+i.*qch7).*(ieqv1+i.*qeqv1));
ich7=icompen;
qch7=qcompen;

%---------- CE symbol removal ----------------

ich8=ich7(:,knd+1:nd+1);
qch8=qch7(:,knd+1:nd+1);

%---------- DC and pilot data removal --------

[ich9,qch9]=crdemapping(ich8,qch8,fftlen,nd);

%----------------- demoduration --------------

ich10=ich9./kmod;
qch10=qch9./kmod;
[demodata]=qpskdemod(ich10,qch10,para,nd,ml);   

%--------------  error calculation  ----------

demodata1=reshape(demodata,1,para*nd*ml);
noe2=sum(abs(demodata1-seldata));
nod2=length(seldata);

% calculating PER
if noe2~=0
   eop=eop+1;
else
   eop=eop;
end   
   eop;
   nop=nop+1;
   
% calculating BER
noe=noe+noe2;
nod=nod+nod2;

fprintf('%d\t%e\t%d\n',iii,noe2/nod2,eop);

end
per=eop/nop;
ber=noe/nod;

%********************** Output result ***************************

fprintf('%f\t%e\t%e\t%d\t%d\n',ci,ber,per,nloop,fd);
  
fid = fopen('BERofdmci.dat','a');
fprintf(fid,'%f\t%e\t%e\t%d\t%d\n',ci,ber,per,nloop,fd);
fclose(fid);

%******************** end of file ***************************



  
 
