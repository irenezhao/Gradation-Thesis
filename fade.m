% Program 2-5
% fade.m
%
% Generate Rayleigh fading
% 
% Programmed by H.Harada	 	

function [iout,qout,ramp,rcos,rsin]=fade(idata,qdata,nsamp,tstp,fd,no,counter,flat)

%****************** variables *************************
% idata  : input Ich data     
% qdata  : input Qch data     
% iout   : output Ich data
% qout   : output Qch data
% ramp   : Amplitude contaminated by fading
% rcos   : Cosine value contaminated by fading
% rsin   : Cosine value contaminated by fading
% nsamp  : Number of samples to be simulated       
% tstp   : Minimum time resolution                    
% fd     : maximum doppler frequency               
% no     : number of waves in order to generate fading   
% counter  : fading counter                          
% flat     : flat fading or not 
% (1->flat (only amplitude is fluctuated),0->nomal(phase and amplitude are fluctutated)    
%******************************************************

if fd ~= 0.0  
    ac0 = sqrt(1.0 ./ (2.0.*(no + 1)));   % power normalized constant(ich)
    as0 = sqrt(1.0 ./ (2.0.*no));         % power normalized constant(qch)
    ic0 = counter;                        % fading counter
 
    pai = 3.14159265;   
    wm = 2.0.*pai.*fd;
    n = 4.*no + 2;
    ts = tstp;
    wmts = wm.*ts;
    paino = pai./no;                        

    xc=zeros(1,nsamp);
    xs=zeros(1,nsamp);
    ic=[1:nsamp]+ic0;

  for nn = 1: no
	  cwn = cos( cos(2.0.*pai.*nn./n).*ic.*wmts );
	  xc = xc + cos(paino.*nn).*cwn;
	  xs = xs + sin(paino.*nn).*cwn;
  end

  cwmt = sqrt(2.0).*cos(ic.*wmts);
  xc = (2.0.*xc + cwmt).*ac0;
  xs = 2.0.*xs.*as0;

  ramp=sqrt(xc.^2+xs.^2);   
  rcos=xc./ramp;
  rsin=xs./ramp;

  if flat ==1
    iout = sqrt(xc.^2+xs.^2).*idata(1:nsamp);    % output signal(ich)
    qout = sqrt(xc.^2+xs.^2).*qdata(1:nsamp);    % output signal(qch)
  else
    iout = xc.*idata(1:nsamp) - xs.*qdata(1:nsamp);   % output signal(ich)
    qout = xs.*idata(1:nsamp) + xc.*qdata(1:nsamp);   % output signal(qch)
  end

else  
  iout=idata;
  qout=qdata;
end

% ************************end of files***********************************
