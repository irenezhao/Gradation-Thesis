function [ber, bits]=chap3_ex04_qpsk_viterbi_soft(EbNo, maxNumErrs, maxNumBits)
%% Initializations
persistent Modulator AWGN DeModulator BitError ConvEncoder Viterbi Quantizer
if isempty(Modulator)
    Modulator      = comm.QPSKModulator('BitInput',true);
    AWGN             = comm.AWGNChannel;
    DeModulator =  comm.QPSKDemodulator('BitOutput',true,...
        'DecisionMethod','Log-likelihood ratio',...
        'VarianceSource', 'Input port');
    BitError           = comm.ErrorRate;
    ConvEncoder=comm.ConvolutionalEncoder(...
        'TerminationMethod','Terminated');
    Viterbi=comm.ViterbiDecoder(...
        'InputFormat','Soft',...
        'SoftInputWordLength', 4,...
        'OutputDataType', 'double',...
        'TerminationMethod','Terminated');
    Quantizer=dsp.ScalarQuantizerEncoder(...
        'Partitioning',  'Unbounded',...
        'BoundaryPoints', -7:7,...
        'OutputIndexDataType','uint8');
end
%% Constants
FRM=2048;
M=4; k=log2(M); codeRate=1/2;
snr = EbNo + 10*log10(k) + 10*log10(codeRate);
noise_var = 10.^(-snr/10);
AWGN.EbNo=snr;
%% Processsing loop modeling transmitter, channel model and receiver
numErrs = 0; numBits = 0;results=zeros(3,1);
while ((numErrs < maxNumErrs) && (numBits < maxNumBits))
    % Transmitter
    u                   = randi([0 1], FRM,1);                                      % Random bits generator
    encoded     = ConvEncoder.step(u);                                     % Convolutional encoder
    mod_sig     = Modulator.step(encoded);                              % QPSK Modulator
    % Channel
    rx_sig          =  AWGN.step(mod_sig);                                  % AWGN channel
    % Receiver
    demod        =  DeModulator.step(rx_sig, noise_var);          % Soft-decision QPSK Demodulator
    llr                 = Quantizer.step(-demod);                              % Quantize Log-Likelihood Ratios    
    decoded     = Viterbi.step(llr);                                              % Viterbi decoder with LLRs
    y                   = decoded(1:FRM);                                          % Compute output bits
    results        = BitError.step(u, y);                                          % Update BER
    numErrs     = results(2);
    numBits      = results(3);
end
%% Clean up & collect results
ber = results(1); bits= results(3);
reset(BitError);