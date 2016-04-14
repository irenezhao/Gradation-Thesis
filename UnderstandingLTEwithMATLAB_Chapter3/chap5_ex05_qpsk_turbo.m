function [ber, bits]=chap5_ex05_qpsk_turbo(EbNo, maxNumErrs, maxNumBits)
%% Constants
FRM=2048;
Trellis=poly2trellis(4, [13 15], 13);
Indices=randperm(FRM);
M=4;k=log2(M);
R= FRM/(3* FRM + 4*3);
snr = EbNo + 10*log10(k) + 10*log10(R);
noise_var = 10.^(-snr/10);
%% Initializations
persistent Modulator AWGN DeModulator BitError TurboEncoder TurboDecoder 
if isempty(Modulator)
    Modulator      = comm.QPSKModulator('BitInput',true);
    AWGN             = comm.AWGNChannel;
    DeModulator =  comm.QPSKDemodulator('BitOutput',true,...
        'DecisionMethod','Log-likelihood ratio',...
        'VarianceSource', 'Input port');
    BitError           = comm.ErrorRate;
    TurboEncoder=comm.TurboEncoder(...
        'TrellisStructure',Trellis,...
        'InterleaverIndices',Indices);
    TurboDecoder=comm.TurboDecoder(...
        'TrellisStructure',Trellis,...
        'InterleaverIndices',Indices,...
        'NumIterations',6);
end
%% Processsing loop modeling transmitter, channel model and receiver
AWGN.EbNo=snr;
numErrs = 0; numBits = 0;results=zeros(3,1);
while ((numErrs < maxNumErrs) && (numBits < maxNumBits))
    % Transmitter
    u                  = randi([0 1], FRM,1);                                       % Random bits generator
    encoded      = TurboEncoder.step(u);                                   % Turbo Encoder
    mod_sig      = Modulator.step(encoded);                             % QPSK Modulator
    % Channel
    rx_sig          =  AWGN.step(mod_sig);                                  % AWGN channel
    % Receiver
    demod         =  DeModulator.step(rx_sig, noise_var);          % Soft-decision QPSK Demodulator
    decoded      = TurboDecoder.step(-demod);                         % Turbo Decoder
    y                  = decoded(1:FRM);                                            % Compute output bits
    results         = BitError.step(u, y);                                          % Update BER
    numErrs       = results(2);
    numBits       = results(3);
end
%% Clean up & collect results
ber = results(1); bits= results(3);
reset(BitError);