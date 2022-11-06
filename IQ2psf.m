function [ sig, hpsf ] = IQ2psf( IQk, Ksample, Nsymb, r, rctype )
% Pulse shaping:
% from vector of many [I(k),Q(k)] carrier states (k=1,2,3,...)
% to transmitted signal IQ(n) coding sequence of many modulation symbols

Npsf = Nsymb*Ksample+1; Mpsf=(Npsf-1)/2;             % PSF filter length and its half
IQ0 = zeros( 1, length(IQk)*Ksample ); IQ0(1:Ksample:end) = IQk;  % zero insertion
hpsf = firrcos( Npsf-1, 1/(2*Ksample), r, 1,'rolloff',rctype);    % 'normal' or 'sqrt'
sig = conv( IQ0, Ksample*hpsf );                                  % pulse shaping
