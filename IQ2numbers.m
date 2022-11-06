function numbers = IQ2numbers( IQ, modtype) 
% from [I,Q] values to carrier state numbers for many input IQ pairs

N = length(IQ);

if( isequal(modtype,'2PAM')  || isequal(modtype,'4PAM')  || isequal(modtype,'8PAM') || ...
    isequal(modtype,'QPSK')  || isequal(modtype,'DQPSK') || isequal(modtype,'8PSK') || ...
    isequal(modtype,'4QAM')  || isequal(modtype,'16QAM') )   
    
    if( isequal(modtype,'2PAM') )   % checking only sign of I
        k = 1;                      % 
        for ns = 1 : N
            if(real(IQ(ns)) < 0)  numbers(k)=0; else numbers(k)=1; end
            k = k + 1;
        end   
    end
    if( isequal(modtype,'4PAM') )   % checking only I
        k = 1;                      % in predefined intervals: -2, 0, 2
        for ns = 1 : N
            I = real( IQ(ns) );
            if(      I < -2          ) numbers(k) = 0;   % lower the -2
            elseif( -2 <= I && I < 0 ) numbers(k) = 1;   % from -2 to 0
            elseif(  0 <= I && I < 2 ) numbers(k) = 3;   % from  0 to 2
            else                       numbers(k) = 2;   % greater then 2
            end
            k = k + 1;
        end   
    end
    if( isequal(modtype,'8PAM') )   % checking only I
        k = 1;                      % in predefined intervals: -4, -2, 0, 2, 4
        for ns = 1 : N
            I = real( IQ(ns) );
            if    (            I < -6 )  numbers(k) = 0;   % lower than -6
            elseif( -6 <= I && I < -4 )  numbers(k) = 1;   % from -6 to -4
            elseif( -4 <= I && I < -2 )  numbers(k) = 3;   % from -4 to -2
            elseif( -2 <= I && I <  0 )  numbers(k) = 2;   % from -2 to  0
            elseif(  0 <= I && I <  2 )  numbers(k) = 6;   % from  0 to  2
            elseif(  2 <= I && I <  4 )  numbers(k) = 7;   % from  2 to  4
            elseif(  4 <= I && I <  6 )  numbers(k) = 5;   % from  4 to  6
            elseif(  6 <= I           )  numbers(k) = 4;   % greater than 6
            end
            k = k + 1;
        end   
    end
    if( isequal(modtype,'4QAM') | isequal(modtype,'QPSK') ) % checking I and Q
        k = 1;                                              % positions in +/- quadrants
        for ns = 1 : N
            I = real( IQ(ns) );
            Q = imag( IQ(ns) );
            if( I>0 && Q>0 )  numbers(k) = 3; end        % +/+ right-up
            if( I>0 && Q<0 )  numbers(k) = 2; end        % +/- right-down
            if( I<0 && Q>0 )  numbers(k) = 1; end        % -/+ left-up
            if( I<0 && Q<0 )  numbers(k) = 0; end        % -/- left-down
            if( I==0 || Q==0 ) disp('ZERO'); pause; end
            k = k + 1;
        end
    end
    if( isequal(modtype,'DQPSK') ) % checking I and Q phase shifts in +/- quadrants
        IQdiff = IQ(2:end) .* conj( IQ(1:end-1));      % finding phase shift
        k = 1;                                 
        for ns = 1 : length(IQdiff)-1
            I = real( IQdiff(ns) );
            Q = imag( IQdiff(ns) );
            if( I>0 && Q>0 )  numbers(k) = 0; end        % +/+ right-up
            if( I<0 && Q>0 )  numbers(k) = 1; end        % -/+ left-up      
            if( I>0 && Q<0 )  numbers(k) = 2; end        % +/- right-down
            if( I<0 && Q<0 )  numbers(k) = 3; end        % -/- left-down  
            if( I==0 || Q==0 ) disp('ZERO'); pause; end
            k = k + 1;
        end
    end
    if( isequal(modtype,'8PSK') )     % checking I and Q
        k = 1;                        % positions in circle octants
        for ns = 1 : N
            I = real( IQ(ns) );
            Q = imag( IQ(ns) );
            if( I>0 && Q>0 )
                if( I > Q)    numbers(k)=7; else numbers(k)=6; end % +/+ right-up      
            end    
            if( I>0 && Q<0 )
                if( I > -Q )  numbers(k)=5; else numbers(k)=4; end % +/- right-down
            end
            if( I<0 && Q>0 )
                if( -I > Q )  numbers(k)=3; else numbers(k)=2; end % -/+ left-up
            end
            if( I<0 && Q<0 )
                if( -I > -Q ) numbers(k)=1; else numbers(k)=0; end % -/- left-dow
            end  
            if( I==0 || Q==0 ) disp('ZERO'); pause; end
            k = k + 1;
        end
    end   
    if( isequal(modtype,'16QAM') )     % checking I and Q
        k = 1;                         % first two higher (MSB) bits of I
        for ns = 1 : N                 % then  two lower  (LSB) bits of Q
            I = real( IQ(ns) );                                    % II 
            if(      I < -2          ) numbers(k) = 0;             % 00
            elseif( -2 <= I && I < 0 ) numbers(k) = 4;             % 01
            elseif(  0 <= I && I < 2 ) numbers(k) = 12;            % 11
            else                       numbers(k) = 8;             % 10
            end
            k = k + 1;
        end
        k = 1;
        for ns = 1 : N
            Q = imag( IQ(ns) );                                    % II + QQ
            if(           Q < -2     )  numbers(k) = numbers(k)+0; % ?? + 00
            elseif( -2 <= Q && Q < 0 )  numbers(k) = numbers(k)+1; % ?? + 01
            elseif(  0 <= Q && Q < 2 )  numbers(k) = numbers(k)+3; % ?? + 11
            else                        numbers(k) = numbers(k)+2; % ?? + 10
            end  
            k = k + 1;
        end
    end
  
else
    disp('Demodulation is not supported!');
end
end

