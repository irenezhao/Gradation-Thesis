%% Constants
FRM=2048;
MaxNumErrs=200;MaxNumBits=1e7;
EbNo_vector=0:10;BER_vector=zeros(size(EbNo_vector));
%% Initializations
Modulator      = comm.QPSKModulator('BitInput',true);
AWGN             = comm.AWGNChannel;
DeModulator =  comm.QPSKDemodulator('BitOutput',true);
BitError           = comm.ErrorRate;
%% Outer Loop computing Bit-error rate as a function of EbNo
for EbNo = EbNo_vector
    snr = EbNo + 10*log10(2);
    AWGN.EbNo=snr;
    numErrs = 0; numBits = 0;results=zeros(3,1);
    %% Inner loop modeling transmitter, channel model and receiver for each EbNo
    while ((numErrs < MaxNumErrs) && (numBits < MaxNumBits))
        % Transmitter
        u             = randi([0 1], FRM,1);                  % Generate random bits
        mod_sig = step(Modulator,   u);                % QPSK Modulator
        % Channel
        rx_sig  = step(AWGN,        mod_sig);        % AWGN channel
        % Receiver
        y =       step(DeModulator, rx_sig);           % QPSK Demodulator
        results = step(BitError,    u, y);                  % Update BER
        numErrs = results(2);
        numBits = results(3);
    end
    % Compute BER
    ber = results(1); bits= results(3);
    %% Clean up & collect results
    reset(BitError);
    BER_vector(EbNo+1)=ber;
end
%% Visualize results
EbNoLin = 10.^(EbNo_vector/10);  
theoretical_results = 0.5*erfc(sqrt(EbNoLin));
semilogy(EbNo_vector, BER_vector)
grid;title('BER vs. EbNo - QPSK modulation');
xlabel('Eb/No (dB)');ylabel('BER');hold;
semilogy(EbNo_vector,theoretical_results,'dr');hold;
legend('Simulation','Theoretical');