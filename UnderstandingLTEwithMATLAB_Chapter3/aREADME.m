% zREADME.m 
% Instructions regarding how to run MATLAB experiments in this directory
% (UnderstandingLTEwithMATLAB_Chapter3)
%
% This folder contains a series of MATLAB functions that showcase 
% modulation, scrambling and coding functionalities used in the LTE standard 
% for DLSCH and PDSCH as presented in chapter 3 of the "Understanding LTE with
% MATLAB"

%
% How to run the demos:
% 1. To perform experiment 1, type chap3_ex01 at the MATLAB command prompt. 
% Notice how simulated and theoretical BER performaqnce measures are
% computed and compared for a simple transceiver composed of a modulator
% and demodulator with an AWGN channel. 
% The experiment shows how to use System objects of
% Communications System Toolbox to set up a basic experiment and verify the
% results.
%
% To perform experiments 2 to 5, use the BERTOOL of the Communications
% System Toolbox as follows:
% type bertool at the MATLAB command prompt. 
% Go to Monte Carlo tab, 
% In the "Simulation MATLAB file ..." edit-box choose either of the
% following main functions
% chap3_ex02_qpsk.m
% chap3_ex03_qpsk_viterbi.m
% chap3_ex04_qpsk_viterbi_soft.m
% chap5_ex05_qpsk_turbo.m
% In each case, set a range of EbNo values, set the BER variable name to ber, and
% set typical values for parameters Number of errors & Number of bits. 
% For better results, the experiments have to be long enough, 
% which usually means parametrers Number of errors or Number of bits
% need to be larger than 1e4 and 1e7 respectively. 
