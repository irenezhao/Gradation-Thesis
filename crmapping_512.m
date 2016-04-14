% Function 4-6
% crmapping.m
%
% Function to set data on subcarrier
% (DC=0)
% 
% Programmed by T.Yamamura and H.Harada
%

function [iout,qout]=crmapping_512(idata,qdata,fftlen,nd);

%****************** variables *************************
% idata     : Input Ich data
% qdata     : Input Qch data
% iout      : Output Ich data
% qout      : Output Qch data
% fftlen    : Length of FFT (points)
% nd        : Number of OFDM symbols
% *****************************************************

iout=zeros(fftlen,nd);
qout=zeros(fftlen,nd);

iout(100:356,:)=idata(1:257,:);
qout(100:356,:)=qdata(1:257,:);


%******************** end of file ***************************
