% Function 4-7
% crdemapping.m
%
% Function to separate data from the subcarrier
% (DC=0)
%
% MATLAB version
% programmed by T.Yamamura
%

function [iout,qout]=crdemapping_512(idata,qdata,fftlen,nd);

%****************** variables *************************
% idata     : Input Ich data
% qdata     : Input Qch data
% iout      : Output Ich data
% qout      : Output Qch data
% fftlen    : Length of FFT (points)
% nd        : Number of OFDM symbols
% *****************************************************

iout(1:257,:)=idata(100:356,:);
qout(1:257,:)=qdata(100:356,:);


%******************** end of file ***************************