% Function 4-10
% interwave.m
%
% Function to add interference wave
%
% Programmed by T.Yamamura and H.Harada
%

function [iout,qout]=interwave(ci,spow,ml,nsamp,tstp,fadingpara);

%****************** variables *************************
% ci         : Carrier to interference ratio
% spow       : Power of desired signal
% ml         : Modulation level
% nsamp      : Number of samples
% tstp       : Time resolution
% fadingpara : Fading parameter
% iout       : Output Ich signal
% qout       : Output Qch signal
% *****************************************************

itau=fadingpara(1,:);
dlvl1=fadingpara(2,:);
n0=fadingpara(3,:);
th1=fadingpara(4,:);
itnd1=fadingpara(5,:);
now1=fadingpara(6,:);
fd=fadingpara(7,:);
flat=fadingpara(8,:);

if ci < 40;

%%%%%%%%%%%%% preparation part %%%%%%%%%%%%%%%%
%%% frame format
para=52;
fftlen=64;
noc=53;  % the number of carrier
nd=6;    % the number of information sysmbol
knd=1;   % the number of known data symbol

sr=250000;  % symbol rate
br=sr.*ml;  % bit rate per carrier
gilen=16;   % the length of guard interval

%%% Set CE data load
kndata=zeros(1,fftlen);
kndata0=2.*(rand(1,52)>0.5)-1;
kndata(2:27)=kndata0(1:26);
kndata(39:64)=kndata0(27:52);

%%% Simulation start
%%% fading initialization


%%%%%%%%%%%%%%%%%%%%%%%%%%%% transmitter %%%%%%%%%%%%%%%%%%%%%%%%%%%
   
seridata=rand(1,para*nd*ml)>0.5;  %  DC=0

paradata=reshape(seridata,para,nd*ml); %size(51  *  nd*ml)

%%% ml modulation 

[ich,qch]=qpskmod(paradata,para,nd,ml);
kmod=1/sqrt(2);
ich=ich.*kmod;
qch=qch.*kmod;

% CE modulation

ceich=kndata; % CE:BPSK
ceqch=zeros(1,64);

%%% data mapping (DC=0)

[ich2,qch2]=crmapping(ich,qch,fftlen,nd);

% addition of pilot carrier and CE symbol

ich22=[ceich.' ich2]; % I-channel transmission data
qch22=[ceqch.' qch2]; % Q-channel transmission data

%%% IFFT

x=ich22+qch22.*i;
y=ifft(x);
ich3=real(y);
qch3=imag(y);

%%% Gurad interval insertion

% guard interval insertion
[ich5,qch5]= giins(ich3,qch3,fftlen,gilen,nd+1);

%%% fading Calculation   
[ifade2,qfade2,ramp,rcos,rsin]=sefade(ich5,qch5,itau,dlvl1,th1,n0,itnd1,now1,length(ich5),tstp,fd,flat);

%%% C/I reduction
spowintw=sum(ich5.^2+qch5.^2)/(nd)/52;
rint=spow/spowintw*10^(-ci/10);
iout=ifade2.*sqrt(rint);
qout=qfade2.*sqrt(rint);
   
else

iout=zeros(1,nsamp);
qout=zeros(1,nsamp);

end

%******************** end of file ***************************

