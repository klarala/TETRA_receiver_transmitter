
clear all; close all;

do_disturb = 0;

text_in = 'This is the tetra signal';

% Parameters
fs = 102400; %2560000;     % sampling frequency 102.4 kHz or 2.56 MHz - written from the file
fsymb = 18000;      % symbol frequency 18 kHz
modtype = 'DQPSK';  % 2PAM, 4PAM, 8PAM, BPSK, QPSK, DQPSK, 8PSK, 4QAM, 16QAM
Nhead = 19;         % STS header length (number of carrier states)
Ndata = 129;        % data length - searched (number of carrier states)
Nframe = 255;       % frame length, including Nhead and Ndata
rctype = 'sqrt';    % PSF: 'sqrt': raised cosine filter type for TX and RX 
r = 0.35;           % PSF: filter roll-off factor
K = 5;              % number of samples per symbol, at present fs/fsymb
                    % but we can change it
Ns = 10;            % PSF: symbols per PS filter 
fcar = 434000000;   % carrier frequency in Hz

% STS 19 header - known
numHead = [ 3 0 0 1 2 1 3 0 3 2 2 1 3 0 0 1 2 1 3 ];

%chan = [ 0.5, -0.25, 0.1, -0.1 ]; % channel impulse response in baseband symbol-spaced
chan = [ 1 ];                      % perfect channel
%SNR=160;  chanG=1; chanPh=0;  carDF=0.0000; carDPh=0;  ADCdt=0;           % No disturb
SNR=40;  chanG=0.25; chanPh=pi/7; carDF=0.0002; carDPh=pi/11; ADCdt=0.5; % Disturb
Mdecim = 1;        % decimation order: 1, 2, 3, 4, 6, 8, 12, 24
Mdelay = 0;        % decimation delay: in samples before decimation  

Npsf = K*Ns+1; Mpsf = (Npsf-1)/2;
[IQcodes, Nstates, Nbits, R ] = IQdef( modtype );       % take carrier IQ codes
% IQk of Header
IQkHead= numbers2IQ( numHead, modtype, IQcodes );       % calculate IQ states
% IQk of Data
numbers = text2numbers( text_in, Nbits );
numbers = [numbers', zeros(129 - size(numbers,1),1)'];
IQkData = numbers2IQ( numbers', modtype, IQcodes );     % IQ state numbers to IQ values
figure;
subplot(211); stem(real(IQkData),'b'); grid; xlabel('k'); title('I(k)');
subplot(212); stem(imag(IQkData),'r'); grid; xlabel('k'); title('Q(k)'); pause
% Numbers ALL, IQk ALL
%num = [ numData numHead numData numHead numData ]; % ALL transmitted IQ numbers
IQk = [ IQkHead  IQkData  IQkHead  IQkData ];  % ALL transmitted IQ states
% IQn of Header only (pulse shaping)
%IQnHead = IQ2psf( IQkHead, K, Ns, r, rctype );         % IQn of Header   
% IQn of everything (pulse shaping)
[IQn, hpsf ] = IQ2psf( IQk, K, Ns, r, rctype );         % IQn of ALL

    N = length( IQn ); n = Npsf : N-Npsf+1; ns = Npsf : K : N-Npsf+1;
    figure; plot( real(IQn(n)), imag(IQn(n)), real(IQn(ns)), imag(IQn(ns)),'ro','MarkerFaceColor','red'); grid; title('TX: Q(n) = f( I(n) )'); pause
    n = n + floor(K/2)+1;
    figure
    subplot(121); plot( reshape(real(IQn(n)),K,length(n)/K),'b'); xlabel('n'); title('TX: Eye diagram for I(n)'); grid;
    subplot(122); plot( reshape(imag(IQn(n)),K,length(n)/K),'r'); xlabel('n'); title('TX: Eye diagram for Q(n)'); grid; pause

data = IQn;
    
% Optional frequency UP conversion in TX plus channel simulation  
if( 0 )
   if( length( chan ) > 1 )         % when more than one channel/filter weight
      chan = resample(chan,K,1);    % upsampling channel impulse response
      figure; plot(chan); title('h(n)'); pause
   end 
   df = 0; dphi = 0;  % CFO added in the base-band OR df = carDF*fs; dphi = carDPh;
   IQn = IQupchandown( IQn, fcar, fs, chan, df, dphi );
end

% Addition of disturbances in the base-band
if( do_disturb )
    IQn = IQdisturb( IQn, SNR, chanG, chanPh, carDF, carDPh, ADCdt, Npsf );
    
        N = length( IQn ); n = Npsf : N-Npsf+1; ns = Npsf : K : N-Npsf+1;
        figure; plot( real(IQn(n)), imag(IQn(n)), real(IQn(ns)), imag(IQn(ns)),'ro','MarkerFaceColor','red'); grid; title('RX: Q(n) = f( I(n) )'); pause
        n = n + floor(K/2)+1;
        figure
        subplot(121); plot( reshape(real(IQn(n)),K,length(n)/K),'b'); xlabel('n'); title('RX: Eye diagram for I(n)'); grid;
        subplot(122); plot( reshape(imag(IQn(n)),K,length(n)/K),'r'); xlabel('n'); title('RX: Eye diagram for Q(n)'); grid; pause
end

IQn2 = repmat(IQn, 1, 100);
IQn_real = real(IQn2) /(max(abs(real(IQn))));
IQn_imag = imag(IQn2) / (max(abs(imag(IQn))));
IQn_new = [IQn_real, IQn_imag];

audiowrite('signal.wav',IQn_new,fs);


