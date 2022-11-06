function text = numbers2text( numbers, Nbits )
% IQ state numbers to text conversion

text = dec2bin( numbers, Nbits );                % state numbers as strings of bits
[rows,cols] = size( text );                      % size of matrix of characters? 
text = reshape( text',[ rows*cols, 1] )';        % one big stream of chars '0'/'1'
N=length(text); N=N-rem(N,8); text = text(1:N);  % remove appended bits
text = reshape( text',[ 8, N/8 ] )';             % strings of bytes
text = strcat( char( bin2dec( text) )' );        % conversion to text
return