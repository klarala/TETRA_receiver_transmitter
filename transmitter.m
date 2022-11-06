
clear all; close all;

do_disturb = 0;

%load('TETRA_423.4125MHz_noise_-22.mat');  K=10;  % files: -22,-10,0,10,22.mat % AWGN (noisy) channel
%load('TETRA_423.4125MHz_flat.mat');  K=5;        % flat channel 
%load('TETRA_423.4125MHz_tu50.mat');  K=8;        % typical urban channel at 50 kph
%load('TETRA_423.4125MHz_ht200.mat'); K=?;        % hilly terrain channel at 200 kph

%IQn = y.'; clear y; fs   % sampling frequency is written from the MAT file 

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
fcar = 145000000;   % carrier frequency in Hz

% Signal resampling
%fsnew = ceil(K)*fsymb;                 % new sampling frequency
%[UP,DOWN] = rat( fsnew/fs );           % UP/DOWN interpolation/decimation orders
%I = resample( real(IQn), UP, DOWN );   % resampling I(n)
%Q = resample( imag(IQn), UP, DOWN );   % resampling Q(n)
%IQn = I + j*Q; clear I Q;              % combining I(n)+j*Q(n)
%fs = fsnew; K = ceil(K);               % setting new values

% STS 19 header - known
numHead = [ 3 0 0 1 2 1 3 0 3 2 2 1 3 0 0 1 2 1 3 ];
% Data to be detected after the 1st header in files -22, -10, 0, 10, 22.mat:
% Length = 129 = 15+108+1+5
numData = [ ...
  2,3,3,2,2,3,1,3,2,1,1,3,0,2,2,1,0,2,3,3,3,3,0,1,0,0,0,1,0,3,0,2,3,...
  3,1,0,2,3,2,1,2,1,2,3,3,3,2,2,3,2,1,1,2,0,0,3,1,2,1,3,2,0,0,2,2,3,...
  1,3,2,3,0,2,1,0,3,1,3,0,1,1,3,1,2,2,2,1,1,2,0,2,1,3,2,0,1,2,1,3,1,...
  3,0,2,1,2,1,0,1,1,3,3,1,3,2,2,2,1,1,3,1,0,2,2,1,3,2,3,1,3,0 ];

%chan = [ 0.5, -0.25, 0.1, -0.1 ]; % channel impulse response in baseband symbol-spaced
chan = [ 1 ];                      % perfect channel
%SNR=160;  chanG=1; chanPh=0;  carDF=0.0000; carDPh=0;  ADCdt=0;           % No disturb
SNR=40;  chanG=0.25; chanPh=pi/7; carDF=0.0002; carDPh=pi/11; ADCdt=0.5; % Disturb
Mdecim = 1;        % decimation order: 1, 2, 3, 4, 6, 8, 12, 24
Mdelay = 0;        % decimation delay: in samples before decimation  

Npsf = K*Ns+1; Mpsf = (Npsf-1)/2;
[IQcodes, Nstates, Nbits, R ] = IQdef( modtype );       % take carrier IQ codes
% IQk of Header
%[numHead, Nhead ] = modtype2header( modtype );          % take header IQ numbers
IQkHead= numbers2IQ( numHead, modtype, IQcodes );       % calculate IQ states
% IQk of Data
%numData = floor( Nstates*(rand(Ndata,1)-10*eps) );      % generate random IQ numbers
IQkData = numbers2IQ( numData, modtype, IQcodes );      % calculate IQ states
% Numbers ALL, IQk ALL
%num = [ numData numHead numData numHead numData ]; % ALL transmitted IQ numbers
IQk = [ IQkData  IQkHead  IQkData  IQkHead  IQkData ];  % ALL transmitted IQ states
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

% Signal and its spectrum
%figure; n = 1:2500;
%subplot(211); plot( real(IQn(n)) ); grid; title('I(n)');
%subplot(212); plot( real(IQn(n)) ); grid; title('Q(n)'); pause
%figure;
%pwelch(IQn,2048,2048-1024,2048,fs,'centered'); pause

