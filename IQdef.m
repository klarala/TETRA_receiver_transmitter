function [ IQcoder, Nstates, Nbits, R ] = IQdef( modtype )
% Costellation tables for different digital modulations
switch modtype
case {'2PAM'},  Nbits=1; IQcoder = [ -1 1 ];
case {'4PAM'},  Nbits=2; IQcoder = [ -3 -1 3 1 ]; 
case {'8PAM'},  Nbits=3; IQcoder = [ -7 -5 -1 -3 7 5 1 3 ];
case {'BPSK'},  Nbits=1; IQcoder = exp(j*pi*[ 1 0 ]);
case {'QPSK'},  Nbits=2; IQcoder = exp(j*pi/4*[ 5 3 7 1 ]);
case {'8PSK'},  Nbits=3; IQcoder = exp(1i*pi/8*[ 11 9 5 7 13 15 3 1 ]);
case {'DQPSK'}, Nbits=2; IQcoder = exp(j*pi/4*[ 1 3 7 5 ]); 
case {'4QAM'},  Nbits=2; IQcoder = [-1-j,-1+j,+1-j,+1+j ];
case {'16QAM'}, Nbits=4; IQcoder = [ -3-3*j, -3-1*j, -3+3*j, -3+1*j, ...
                                     -1-3*j, -1-1*1i, -1+3*j, -1+1*j, ...
                                     +3-3*j, +3-1*j, +3+3*j, +3+1*j, ...
                                     +1-3*j, +1-1*j, +1+3*j, +1+1*j];
otherwise disp('Unknown modulation type'); return; 
end
Nstates = 2^Nbits;
if( strcmp(modtype,'4QAM') || strcmp(modtype,'16QAM') ) R=sqrt(2); else R=1; end
return
