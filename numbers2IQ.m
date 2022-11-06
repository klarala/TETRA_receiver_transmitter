function IQk = numbers2IQ( numbers, modtype, IQstates )
% State numbers to IQ values

if( isequal( modtype, 'DQPSK' ) )  % differential coding only for DQPSK
     IQk(1) = exp(j*0);            % initial IQ state
     for k = 1:length(numbers)-1   % loop
         IQk(k+1) = IQk(k) * IQstates( numbers(k)+1 );   % next IQ state
     end
else IQk = IQstates( numbers+1 );
end
end