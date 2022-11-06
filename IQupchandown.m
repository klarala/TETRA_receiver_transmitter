function IQn = IQupchandown( IQn, fcar, fs, chan, df, dphi )
% IQ signal frequency-UP and DOWN conversion with low-pass filtering in the RX

N = length(IQn); n = 0:N-1;                 % signal length, sample indexes

% UP - quadrature modulator
y = real(IQn).*cos(2*pi*fcar/fs*n) - imag(IQn).*sin(2*pi*fcar/fs*n);

% CHANNEL
if( length(chan) > 1 ) y = conv( y, chan, 'same' ); 
else                   y = chan(1) * y;
end    

% DOWN  
%IQn = 2* y .* exp( -j*( 2*pi*(fcar+df)/fs*n + dphi ) );

% Filtering of the (2*fcar) component
%if(1) % 0/1
 % L = 250;                                             % choice of filter order
  %hLP = fir1(L,fcar/(fs/2),kaiser(L+1,10));            % design of filter weigths
  %IQn = conv( IQn, hLP ); IQn = IQn( L/2+1: end-L/2);  % low-pass filtering
%end

% figure; spectrogram(IQn,2048,2048-256,2048,fs); pause

