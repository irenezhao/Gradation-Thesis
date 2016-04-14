% Function 4-6
% crmapping.m
%
% Function to set data on subcarrier
% (DC=0)
% 
% Programmed by T.Yamamura and H.Harada
%

function [iout,qout]=crmapping(idata,qdata,fftlen,nd);

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

iout(2:27,:)=idata(1:26,:);
qout(2:27,:)=qdata(1:26,:);
iout(39:64,:)=idata(27:52,:);
qout(39:64,:)=qdata(27:52,:);

%******************** end of file ***************************
