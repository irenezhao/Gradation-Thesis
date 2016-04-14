% Program 2-4
% comb.m
%
% Generate additive white gausian noise
% 
% Programmed by H.Harada	 	
% 

function [iout,qout] = comb (idata,qdata,attn)

%****************** variables *************************
% idata : input Ich data
% qdata : input Qch data
% iout   output Ich data
% qout   output Qch data
% attn : attenuation level caused by Eb/No or C/N
%******************************************************

iout = randn(1,length(idata)).*attn;
qout = randn(1,length(qdata)).*attn;
 
iout = iout+idata(1:length(idata));
qout = qout+qdata(1:length(qdata));

% ************************end of file***********************************
