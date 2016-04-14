%
% README Chapter 4
%
% by Hiroshi Harada
%
% If you have any bugs and questions in our simulation programs, please e-mail
% to harada@ieee.org. We try to do our best to answer your questions.
%

In this directory, we can find the sixteen files. The relationship between file name and the number of program written in the book is shown in as follows.

Program4-1   ofdm.m
Program4-2   ofdm_fading.m
Program4-3   giins.m
Program4-4   girem.m
Program4-5   ofdmda.m
Program4-6   crmapping.m
Program4-7   crdemapping.m
Program4-8   ofdmce.m
Program4-9   ofdmci.m
Program4-10  interwave.m
Program2-4   comb.m
Program2-5   fade.m
Program2-6   sefade.m
Program2-7   delay.m
Program3-9   qpskmod.m
Program3-10  qpskdemod.m

If you would like to try to use the above programs by using MATLAB. First of all, please copy all of files to your created adequate directory. Then, you start to run MATLAB and you can see the following command prompt in the command window.

>>

Next, you can go to the directory that have all of programs in this section by using change directory (cd) command. If you copy all of files to /matlabR12/work/chapter4, you only type the following command.

>> cd /matlabR12/work/chapter4

In this directory, we can find five main functions, ofdm.m, ofdm_fading.m, ofdmda.m, ofdmce.m, and ofdmci.m.

##################################################
(I) Simulation of "ofdm.m"

This program simulates bit error rate (BER) and packet error rate (PER) performances under AWGN environment.

(1) Set paremeters

First of all, we set simulation parameters in "ofdm.m".

para=128;   % Number of parallel channel to transmit (points)
fftlen=128; % FFT length
noc=128;    % Number of carrier
nd=6;       % Number of information OFDM symbol for one loop
ml=2;       % Modulation level : QPSK
sr=250000;  % OFDM symbol rate (256 ksyombol/s)
br=sr.*ml;  % Bit rate per carrier
gilen=32;   % Length of guard interval (points)
ebn0=3;     % Eb/N0
nloop=100;  % Number of simulation loops

(2) Type just the following command

>> clear
>> ofdm

(3) Then, you can see the following simulation result on your command window.
(example)
3.000000	3.694661e-002	1.000000e+000	100	

where first number 3 is Eb/No, second number 3.694661e-002 is the BER performance, third number 1.000000e+000 is PER performance, and fourth number 100 is Number of simulation loops. And, simulation result is stored in the file (BERofdm.dat) .

##################################################
(II) Simulation of "ofdm_fading.m"

This program simulates bit error rate (BER) and packet error rate (PER) performances under fading environment.

(1) Set paremeters

First of all, we set simulation parameters in "ofdm_fading.m".

para=128;   % Number of parallel channel to transmit (points)
fftlen=128; % FFT length
noc=128;    % Number of carrier
nd=6;       % Number of information OFDM symbol for one loop
ml=2;       % Modulation level : QPSK
sr=250000;  % OFDM symbol rate
br=sr.*ml;  % Bit rate per carrier
gilen=32;   % Length of guard interval (points)
ebn0=10;    % Eb/N0

% Time resolution
tstp=1/sr/(fftlen+gilen); 

% Arrival time for each multipath normalized by tstp
itau = [0];

% Mean power for each multipath normalized by direct wave.
dlvl = [0];

% Number of waves to generate fading for each multipath.
n0=[6];

% Initial Phase of delayed wave
th1=[0.0];

% Number of fading counter to skip 
itnd0=nd*(fftlen+gilen)*10;

% Initial value of fading counter
itnd1=[1000];

% Number of directwave + Number of delayed wave
now1=1;        

% Maximum Doppler frequency [Hz]
fd=320;       

% You can decide two mode to simulate fading by changing the variable flat
% flat     : flat fading or not 
% (1->flat (only amplitude is fluctuated),0->nomal(phase and amplitude are fluctutated)
flat =1;

nloop=500;  % Number of simulation loops


(2) Type just the following command

>> clear
>> ofdm_fading

(3) Then, you can see the following simulation result on your command window.
(example)
10.000000	2.153255e-002	4.660000e-001	500		

where first number 10.000000 is Eb/No, second number 2.153255e-002 is the BER performance, third number 4.660000e-001 is PER performance, and fourth number 500 is Number of simulation loops. And, simulation result is stored in the file (BERofdmfad.dat) .

##################################################
(III) Simulation of "ofdmda.m"

This program simulates bit error rate (BER) and packet error rate (PER) performances under AWGN and fading environments. In this simulation, there are six OFDM
symbols in one packet, and number of FFT length is 64 and number of parallel channel is 52.

(1) Set paremeters

First of all, we set simulation parameters in "ofdmda.m".

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

nloop=1000;   % Number of simulation loops

(2) Type just the following command

>> clear
>> ofdmda

(3) Then, you can see the following simulation result on your command window.
(example)
3.000000	3.704968e-002	1.000000e+000	1000	150

where first number 3 is Eb/No, second number 3.704968e-002 is the BER performance, third number 1.000000e+000 is PER performance, fourth number 1000 is number of simulation loops, and fifth number 150 is Doppler frequency. And, simulation result is stored in the file (BERofdmda.dat) . 

(4) The default mode of "ofdmda.m" is BER and PER under AWGN environment. Therefore, the fifth number in the simulation result is nothing to mean. If you would like to simulate BER and PER performances under fading environment. You must remove "*" from the top of the following sentences in OFDMda.m.

%[ifade,qfade,ramp,rcos,rsin]=sefade(ich4,qch4,itau,dlvl1,th1,n0,itnd1,now1,length(ich4),tstp,fd,flat);
%itnd1 = itnd1+itnd0;  % Updata fading counter
%ich4=ifade;
%qch4=qfade;

Then, type the following command again.

>> clear
>> ofdmda

Then, you can see the following simulation result on your command window.
(example when ebn0=10, flat=0, fd=150, and the others are same of the first simulation)

10.000000	4.986442e-001	9.320000e-001	1000	150

(5) The above simulation results show that BER=0.5. This means we need compensation of fading. If you would like to simulate the effect of perfect compensation under fading environment. You must remove "*" from the top of the following sentences in OFDMda.m.

%[ifade,qfade,ramp,rcos,rsin]=sefade(ich4,qch4,itau,dlvl1,th1,n0,itnd1,now1,length(ich4),tstp,fd,flat);
%itnd1 = itnd1+itnd0;  % Updata fading counter
%ich4=ifade;
%qch4=qfade;

%ifade2=1./ramp.*(rcos(1,:).*ich5+rsin(1,:).*qch5);
%qfade2=1./ramp.*(-rsin(1,:).*ich5+rcos(1,:).*qch5);
%ich5=ifade2;
%qch5=qfade2;

This means that we obtain the information of fading from the fading simulator,
and multiply the complex reciprocal by the received signal.

Then, type the following command again.

>> clear
>> ofdmda

Then, you can see the following simulation result on your command window.
(example when ebn0=10, flat=0, fd=150, and the others are same of the first simulation)

10.000000	2.121795e-002	4.020000e-001	1000	150

From the above result, you can see the effect of the fading compesation.

##################################################
(IV) Simulation of "ofdmce.m"

This program simulates bit error rate (BER) and packet error rate (PER) performances under AWGN and fading environments. In this simulation, there are one OFDM channel estimation (CE) symbols and six OFDM symbols in one packet, and number of FFT length is 64, number of parallel channel is 52 and number of carriers is 53.

(1) Set paremeters

First of all, we set simulation parameters in "ofdmce.m".

%********************** preparation part ***************************

para=52;    % Number of parallel channel to transmit (points)
fftlen=64;  % FFT length
noc=53;     % Number of carriers
nd=6;       % Number of information OFDM symbol for one loop
knd=1;      % Number of known channel estimation (CE) OFDM symbol
ml=2;       % Modulation level : QPSK
sr=250000;  % OFDM symbol rate (250 ksyombol/s)
br=sr.*ml;  % Bit rate per carrier
gilen=16;   % Length of guard interval (points)
ebn0=3;     % Eb/N0

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

(2) Type just the following command

>> clear
>> ofdmce

(3) Then, you can see the following simulation result on your command window.
(example)
3.000000	4.970833e-002	1.000000e+000	1000	150

where first number 3 is Eb/No, second number 4.970833e-002 is the BER performance, third number 1.000000e+000 is PER performance, fourth number 1000 is number of simulation loops, and fifth number 150 is Doppler frequency. And, simulation result is stored in the file (BERofdmce.dat) . 

(4) The default mode of "ofdmda.m" is BER and PER under AWGN environment. Therefore, the fifth number in the simulation result is nothing to mean. If you would like to simulate BER and PER performances under fading environment. You must remove "*" from the top of the following sentences in OFDMda.m.

%[ifade,qfade,ramp,rcos,rsin]=sefade(ich4,qch4,itau,dlvl1,th1,n0,itnd1,now1,length(ich4),tstp,fd,flat);
%itnd1 = itnd1+itnd0;  % Updata fading counter
%ich4=ifade;
%qch4=qfade;

Then, type the following command again.

>> clear
>> ofdmce

Then, you can see the following simulation result on your command window.
(example when ebn0=10, flat=0, fd=150, and the others are same of the first simulation)

10.000000	4.981266e-001	9.360000e-001	1000	150

(5) The above simulation results show that BER=0.5. This means we need compensation of fading. If you would like to simulate the effect of perfect compensation under fading environment. You must remove "*" from the top of the following sentences in OFDMda.m.

%[ifade,qfade,ramp,rcos,rsin]=sefade(ich4,qch4,itau,dlvl1,th1,n0,itnd1,now1,length(ich4),tstp,fd,flat);
%itnd1 = itnd1+itnd0;  % Updata fading counter
%ich4=ifade;
%qch4=qfade;

%ifade2=1./ramp.*(rcos(1,:).*ich5+rsin(1,:).*qch5);
%qfade2=1./ramp.*(-rsin(1,:).*ich5+rcos(1,:).*qch5);
%ich5=ifade2;
%qch5=qfade2;

This means that we obtain the information of fading from the fading simulator,
and multiply the complex reciprocal by the received signal.

Then, type the following command again.

>> clear
>> ofdmce

Then, you can see the following simulation result on your command window.
(example when ebn0=10, flat=0, fd=150, and the others are same of the first simulation)

10.000000	2.455449e-002	4.660000e-001	1000	150

From the above result, you can see the effect of the fading compesation.

(6) Instead of perfect fading compensation, we can use CE OFDM symbol to 
estimate fading and the estimation result is used to compensate the signal fluctuation of fading. If you would like to simulate the effect of CE-symbol under fading environment. You must remove "*" from the top of the following sentences in OFDMda.m.

%[ifade,qfade,ramp,rcos,rsin]=sefade(ich4,qch4,itau,dlvl1,th1,n0,itnd1,now1,length(ich4),tstp,fd,flat);
%itnd1 = itnd1+itnd0;  % Updata fading counter
%ich4=ifade;
%qch4=qfade;

%ce=1;
%ice0=ich2(:,ce);
%qce0=qch2(:,ce);

%ice1=ich7(:,ce);
%qce1=qch7(:,ce);

%iv=real((1./(ice1.^2+qce1.^2)).*(ice0+i.*qce0).*(ice1-i.*qce1));
%qv=imag((1./(ice1.^2+qce1.^2)).*(ice0+i.*qce0).*(ice1-i.*qce1));

%ieqv1=[iv iv iv iv iv iv iv];
%qeqv1=[qv qv qv qv qv qv qv];

%icompen=real((ich7+i.*qch7).*(ieqv1+i.*qeqv1));
%qcompen=imag((ich7+i.*qch7).*(ieqv1+i.*qeqv1));
%ich7=icompen;
%qch7=qcompen;

Then, you can see the following simulation result on your command window.
(example when ebn0=10, flat=0, fd=150, and the others are same of the first simulation)

10.000000	4.992468e-002	6.590000e-001	1000	150

From the above result, you also see the effect of the fading compesation based on CE OFDM symbol.

##################################################
(V) Simulation of "ofdmci.m"

(I) Simulation of "ofdmci.m"

This program simulates bit error rate (BER) and packet error rate (PER) performances when interference wave exists under fading environment.In this simulation, there are one OFDM channel estimation (CE) symbols and six OFDM symbols in one packet, and number of FFT length is 64, number of parallel channel is 52 and number of carriers is 53.

(1) Set paremeters

First of all, we set simulation parameters in "ofdmci.m".

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
nloop=1000;  % Number of simulation loops

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

Moreover, we must set simulation parameters for the interference wave.

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


(2) Type just the following command

>> clear
>> ofdmci

(3) Then, you can see the following simulation result on your command window.

10.000000	6.388462e-002	3.150000e-001	1000	150

where first number 10 is Eb/N0, second number 6.388462e-002 is the BER performance, third number 3.150000e-001 is PER performance, and fourth number 1000 is Number of simulation loops. And, simulation result is stored in the file (BERofdmci.dat) .

Default simulation used CE-based fading compensation. But you can do the same whe used perfect fading compensation by using the same procedure in the simulation of "ofdmce.m".

By changing the value of Eb/N0 (variable ebn0), you can obtain the graph that shows the relationship between Eb/N0 and BER or Eb/N0 and PER and that can been seen in the figures of the book.

********** end of file ********** 

