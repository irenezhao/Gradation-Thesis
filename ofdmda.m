% Program 4-5
% ofdmda.m
%
% Simulation program to realize OFDM transmission system
%
% Programmed by T.Yamamura and H.Harada
% Prease refer Program 4-2
%
% GI data GI data ....(data 6symbols)
%


%%%%%%%%%%%%% preparation part %%%%%%%%%%%%%%%%%

para=52;     % Number of parallel channel to transmit (points)
fftlen=64;   % FFT length
noc=53;      % Number of carriers
nd=6;        % Number of information OFDM symbol for one loop
ml=2;        % Modulation level : QPSK
sr=250000;   % OFDM symbol rate (250 ksyombol/s)
br=sr.*ml;   % Bit rate per carrier
gilen=16;    % Length of guard interval (points)
ebn0=3;      % Eb/N0

%%%%%%%%%%%%% fading initialization %%%%%%%%%%%

tstp=1/sr/(fftlen+gilen); % Time resolution
itau=[0];       % Arrival time for each multipath normalized by tstp 
dlvl1=[0];      % Mean power for each multipath normalized by direct wave.
n0=[6];	        % Number of waves to generate fading n0(1),n0(2)
th1=[0.0];      % Initial Phase of delayed wave
itnd1=[1000];   % set fading counter        	
now1=1;         % Number of directwave + Number of delayed wave
fd=150;         % Maximum Doppler frequency
flat=0;         % Flat or not (see ofdm_fading.m)
itnd0=nd*(fftlen+gilen)*20; % Number of fading counter to skip 

%************************** main loop part **************************

nloop=1000;  % Number of simulation loops

noe = 0;    % Number of error data
nod = 0;    % Number of transmitted data
eop=0;      % Number of error packet
nop=0;      % Number of transmitted packet

%************************** transmitter *****************************
for iii=1:nloop
   
seridata=rand(1,para*nd*ml)>0.5;  %  rand : built in function
paradata=reshape(seridata,para,nd*ml); %  reshape : built in function

%-------------- ml modulation ---------------- 

[ich,qch]=qpskmod(paradata,para,nd,ml);
kmod=1/sqrt(2); %  sqrt : built in function
ich=ich.*kmod;
qch=qch.*kmod;

%------------- data mapping (DC=0) -----------

[ich1,qch1]=crmapping(ich,qch,fftlen,nd);

%------------------- IFFT  -------------------

x=ich1+qch1.*i;
y=ifft(x);      %  ifft : built in function
ich2=real(y);   %  real : built in function
qch2=imag(y);   %  imag : built in function

%---------- Gurad interval insertion ---------
fftlen2=fftlen+gilen;
[ich4,qch4]= giins(ich2,qch2,fftlen,gilen,nd);

%---------- Attenuation Calculation ----------
spow=sum(ich4.^2+qch4.^2)/nd./para;  %  sum : built in function
attn=0.5*spow*sr/br*10.^(-ebn0/10);
attn=sqrt(attn);


%********************** fading channel ******************************   
%If you would like to simulate performance under fading, please remove "*"
%from the following four sentenses
%[ifade,qfade,ramp,rcos,rsin]=sefade(ich4,qch4,itau,dlvl1,th1,n0,itnd1,now1,length(ich4),tstp,fd,flat);
%itnd1 = itnd1+itnd0;  % Updata fading counter
%ich4=ifade;
%qch4=qfade;


%***************************  Receiver  *****************************
%--------------- AWGN addition --------------- 
[ich5,qch5]=comb(ich4,qch4,attn);

%----1 Path Rayleigh Fading perfect compensation ----
%If you would like to simulate performance under perfect compensation, please remove "*"
%from the following four sentenses
%ifade2=1./ramp.*(rcos(1,:).*ich5+rsin(1,:).*qch5);
%qfade2=1./ramp.*(-rsin(1,:).*ich5+rcos(1,:).*qch5);
%ich5=ifade2;
%qch5=qfade2;

%----------- Guard interval removal ----------
[ich6,qch6]= girem(ich5,qch5,fftlen2,gilen,nd);

%------------------  FFT  --------------------
rx=ich6+qch6.*i;
ry=fft(rx);   %  fft : built in function
ich7=real(ry);
qch7=imag(ry);

%----------------- demoduration --------------
[ich8,qch8]=crdemapping(ich7,qch7,fftlen,nd);
ich9=ich8./kmod;
qch9=qch8./kmod;
[demodata]=qpskdemod(ich9,qch9,para,nd,ml);   


%--------------  error calculation  ----------
demodata1=reshape(demodata,1,para*nd*ml);
noe2=sum(abs(demodata1-seridata));  %  sum : built in function
nod2=length(seridata);  %  length : built in function

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

fprintf('%d\t%e\t%d\n',iii,noe2/nod2,eop);  %  fprintf : built in function
   
end
per=eop/nop;
ber=noe/nod;

%********************** Output result ***************************

fprintf('%f\t%e\t%e\t%d\t%d\n',ebn0,ber,per,nloop,fd);
fid = fopen('BERofdmda.dat','a');
fprintf(fid,'%f\t%e\t%e\t%d\t%d\t\n',ebn0,ber,per,nloop,fd);
fclose(fid);

%******************** end of file ***************************