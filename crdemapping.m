% Function 4-7
% crdemapping.m
%
% Function to separate data from the subcarrier
% (DC=0)
%
% MATLAB version
% programmed by T.Yamamura
%

function [iout,qout]=crdemapping(idata,qdata,fftlen,nd);

%****************** variables *************************
% idata     : Input Ich data
% qdata     : Input Qch data
% iout      : Output Ich data
% qout      : Output Qch data
% fftlen    : Length of FFT (points)
% nd        : Number of OFDM symbols
% *****************************************************

iout(1:26,:)=idata(2:27,:);
qout(1:26,:)=qdata(2:27,:);
iout(27:52,:)=idata(39:64,:);
qout(27:52,:)=qdata(39:64,:);

%******************** end of file ***************************