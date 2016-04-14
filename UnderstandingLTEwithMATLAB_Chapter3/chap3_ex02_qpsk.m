function [ber, bits]=chap3_ex02_qpsk(EbNo, maxNumErrs, maxNumBits)
%% Initializations
persistent Modulator AWGN DeModulator BitError
if isempty(Modulator)
    Modulator      = comm.QPSKModulator('BitInput',true);
    AWGN             = comm.AWGNChannel;
    DeModulator =  comm.QPSKDemodulator('BitOutput',true);
    BitError           = comm.ErrorRate;
end
%% Constants
FRM=2048;
M=4; k=log2(M); 
snr = EbNo + 10*log10(k);
AWGN.EbNo=snr;
%% Processsing loop modeling transmitter, channel model and receiver
numErrs = 0; numBits = 0;results=zeros(3,1);
while ((numErrs < maxNumErrs) && (numBits < maxNumBits))
    % Transmitter
    u                   = randi([0 1], FRM,1);                   % Random bits generator
    mod_sig     = Modulator.step(u);                    % QPSK Modulator
    % Channel
    rx_sig          =  AWGN.step(mod_sig);             % AWGN channel
    % Receiver
    demod        =  DeModulator.step(rx_sig);      % QPSK Demodulator
    y                   = demod(1:FRM);                         % Compute output bits
    results        = BitError.step(u, y);                       % Update BER
    numErrs     = results(2);
    numBits      = results(3);
end
%% Clean up & collect results
ber = results(1); bits= results(3);
reset(BitError);