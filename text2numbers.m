function [numbers, bitsnum, bitschar] = text2numbers( text, Nbits )
% text to IQ state numbers conversion

text

bitschar = dec2bin( double(text), 8 );       % text array, letters in rows,'0'/'1'
[rows,cols] = size( bitschar );              % matrix size
N = rows*cols;                               % number of all bits
work = reshape( bitschar',[ 1, N] )';        % bits in one column
Nadd = Nbits-rem(N,Nbits);                   % lacking bits for the last state
for k=1:Nadd, work=[work;'0']; end           % appending '0' bits at the end
bitsnum = reshape( work',[ Nbits, (N+Nadd)/Nbits] )'; % bits of all states
numbers = bin2dec( bitsnum );                % state numbers: from 0 to 2^Nbits-1
return