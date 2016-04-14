function [ber,err,bits]=chap3_ex03_qpsk_viterbi(EbNo)
%% Initializations
persistent Modulator AWGN DeModulator BitError ConvEncoder Viterbi
if isempty(Modulator)
    Modulator      = comm.QPSKModulator('BitInput',true);
    AWGN             = comm.AWGNChannel;
    DeModulator =  comm.QPSKDemodulator('BitOutput',true);
    BitError           = comm.ErrorRate;
    ConvEncoder=comm.ConvolutionalEncoder(...
        'TerminationMethod','Terminated');
    Viterbi=comm.ViterbiDecoder('InputFormat','Hard',...
        'TerminationMethod','Terminated');
end
%% Constants
FRM=2048;
M=4; k=log2(M); codeRate=1/2;
snr = EbNo + 10*log10(k) + 10*log10(codeRate);
AWGN.EbNo=snr;
%% Processsing loop modeling transmitter, channel model and receiver
% numErrs = 0; numBits = 0;
results=zeros(3,1);
% while ((numErrs < maxNumErrs) && (numBits < maxNumBits))
    % Transmitter
    for i=1:2
    u                   = randi([0 1], FRM,1);                % Random bits generator
    encoded     = ConvEncoder.step(u);               % Convolutional encoder
    mod_sig     = Modulator.step(encoded);       % QPSK Modulator
    % Channel
    rx_sig          =  AWGN.step(mod_sig);           % AWGN channel
    % Receiver
    demod        =  DeModulator.step(rx_sig);      % QPSK Demodulator
    decoded     = Viterbi.step(demod);                % Viterbi decoder
    y                   = decoded(1:FRM);                    % Compute output bits
    results        = BitError.step(u, y);       % Update BER

    end
%    numErrs     = results(2);
%    numBits      = results(3);
% end
%% Clean up & collect results
ber = results(1); bits= results(3);
err=results(2);
reset(BitError);