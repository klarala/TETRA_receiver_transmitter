function IQn = IQdisturb( IQn, SNR, chanG, chanPh, carDF, carDPh, ADCdt, Nskip )
% Disturbances - noise, channel gain and phase shift, RX carrier/ADC offsets 
% ADCdt - AD converter time offset in fraction of sampling period, e.g. 0.5

N = length( IQn );

% Noise - noise addition with SNR 
  s = IQn(Nskip : end-Nskip+1);
  scale = sqrt( sum( s.*conj(s) ) / (2*length(s)*10^(SNR/10)) ); clear s
  IQn = IQn + scale * (randn(1,N)+j*randn(1,N));
% IQn = awgn( IQn, SNR, 'measured' ); % alternative noise addition

% Channel - equivalent channel in the BB
  IQn = IQn .* ( chanG * exp(j*chanPh ) );       
      
% Carrier and ADC errors 
  IQn = IQn .* exp(j*(2*pi*carDF*(0:length(IQn)-1)+carDPh));
  
% ADC sampling moments error
if( ADCdt ~= 0 )
    IQn(1,1:N)=interp1( [0:N-1], IQn(1,1:N), [ADCdt : 1 : ADCdt+N-1], 'spline' );
end


