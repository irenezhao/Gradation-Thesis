% Program 2-6
% sefade.m
%
% This function generates frequency selecting fading
% 
% Programmed by H.Harada	 	
% 

function[iout,qout,ramp,rcos,rsin]=sefade(idata,qdata,itau,dlvl,th,n0,itn,n1,nsamp,tstp,fd,flat)

%****************** variables *************************
% idata  input Ich data     
% qdata  input Qch data     
% iout   output Ich data
% qout   output Qch data
% ramp   : Amplitude contaminated by fading
% rcos   : Cosine value contaminated by fading
% rsin   : Cosine value contaminated by fading
% itau   : Delay time for each multipath fading
% dlvl   : Attenuation level for each multipath fading
% th     : Initialized phase for each multipath fading
% n0     : Number of waves in order to generate each multipath fading
% itn    : Fading counter for each multipath fading
% n1     : Number of summation for direct and delayed waves 
% nsamp   : Total number od symbols
% tstp   : Mininum time resolution
% fd   : Maxmum doppler frequency
% flat     flat fading or not 
% (1->flat (only amplitude is fluctuated),0->nomal(phase and amplitude are fluctutated)   
%******************************************************

iout = zeros(1,nsamp);
qout = zeros(1,nsamp);

total_attn = sum(10 .^( -1.0 .* dlvl ./ 10.0));

for k = 1 : n1 

	atts = 10.^( -0.05 .* dlvl(k));

	if dlvl(k) >= 40.0 
	       atts = 0.0;
	end

	theta = th(k) .* pi ./ 180.0;	

	[itmp,qtmp] = delay ( idata , qdata , nsamp , itau(k));
	[itmp3,qtmp3,ramp,rcos,rsin] = fade (itmp,qtmp,nsamp,tstp,fd,n0(k),itn(k),flat);
	
  iout = iout + atts .* itmp3 ./ sqrt(total_attn);
  qout = qout + atts .* qtmp3 ./ sqrt(total_attn);

end
% ************************end of file***********************************
