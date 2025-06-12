
function deciding_result = deciding(bits_in, modulation)

full_len = length(bits_in);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BPSK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(findstr(modulation, 'BPSK'))
   sig_temp=[-1 1];
	for i=1:full_len
        for ii=1:2
            if ii==1
	           dot_result=dot(bits_in(1,i),sig_temp(1,ii));
            else
               dot_result=[dot_result dot(bits_in(1,i),sig_temp(1,ii))];
            end
        end

        dot_result_max=max(real(dot_result));
        
        if real(dot_result(1,1))==dot_result_max
           deciding_result_temp=sig_temp(1,1); 
           
        elseif real(dot_result(1,2))==dot_result_max
           deciding_result_temp=sig_temp(1,2); 
                     
        end

        if i==1
           deciding_result=deciding_result_temp;   
        else
           deciding_result=[deciding_result deciding_result_temp];
        end
	end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% QPSK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ~isempty(findstr(modulation, 'QPSK'))
    sig_temp=exp(j*[-3/4*pi 3/4*pi 1/4*pi -1/4*pi]);
	for i=1:full_len
        for ii=1:4
            if ii==1
	           dot_result=dot(bits_in(1,i),sig_temp(1,ii));
            else
               dot_result=[dot_result dot(bits_in(1,i),sig_temp(1,ii))];
            end
        end

        dot_result_max=max(real(dot_result));
        
        if real(dot_result(1,1))==dot_result_max
           deciding_result_temp=sig_temp(1,1); 
           
        elseif real(dot_result(1,2))==dot_result_max
           deciding_result_temp=sig_temp(1,2); 
           
        elseif real(dot_result(1,3))==dot_result_max
           deciding_result_temp=sig_temp(1,3); 
           
        elseif real(dot_result(1,4))==dot_result_max
           deciding_result_temp=sig_temp(1,4); 
           
        end

        if i==1
           deciding_result=deciding_result_temp;   
        else
           deciding_result=[deciding_result deciding_result_temp];
        end
	end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 8PSK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ~isempty(findstr(modulation, '8PSK'))
    sig_temp=exp(j*[1/8*pi 3/8*pi 7/8*pi 5/8*pi -1/8*pi -3/8*pi -7/8*pi -5/8*pi]);
	for i=1:full_len
        for ii=1:8
            if ii==1
	           dot_result=dot(bits_in(1,i),sig_temp(1,ii));
            else
               dot_result=[dot_result dot(bits_in(1,i),sig_temp(1,ii))];
            end
        end

        dot_result_max=max(real(dot_result));
        
        if real(dot_result(1,1))==dot_result_max
           deciding_result_temp=sig_temp(1,1); 
           
        elseif real(dot_result(1,2))==dot_result_max
           deciding_result_temp=sig_temp(1,2); 
           
        elseif real(dot_result(1,3))==dot_result_max
           deciding_result_temp=sig_temp(1,3); 
           
        elseif real(dot_result(1,4))==dot_result_max
           deciding_result_temp=sig_temp(1,4); 
           
        elseif real(dot_result(1,5))==dot_result_max
           deciding_result_temp=sig_temp(1,5); 
           
        elseif real(dot_result(1,6))==dot_result_max
           deciding_result_temp=sig_temp(1,6); 
           
        elseif real(dot_result(1,7))==dot_result_max
           deciding_result_temp=sig_temp(1,7); 
           
        elseif real(dot_result(1,8))==dot_result_max
           deciding_result_temp=sig_temp(1,8); 
           
        end

        if i==1
           deciding_result=deciding_result_temp;   
        else
           deciding_result=[deciding_result deciding_result_temp];
        end
	end

% elseif ~isempty(findstr(modulation, '16QAM'))
%    % generates 16QAM symbols
%    m=1;
%    for k=-3:2:3
%       for l=-3:2:3
%          table(m) = (k+j*l)/sqrt(10); % power normalization
%          m=m+1;
%       end;
%    end;
%    table=table([0 1 3 2 4 5 7 6 12 13 15 14 8 9 11 10]+1); % Gray code mapping pattern for 8-PSK symbols
%    inp=reshape(bits_in,4,full_len/4);
%    mod_symbols=table([8 4 2 1]*inp+1);  % maps transmitted bits into 16QAM symbols
%    
%    % 64-QAM modulation
% elseif ~isempty(findstr(modulation, '64QAM'))
%    % generates 64QAM symbols
%    m=1;
%    for k=-7:2:7
%       for l=-7:2:7
%          table(m) = (k+j*l)/sqrt(42); % power normalization
%          m=m+1;
%       end;
%    end;
%    table=table([[ 0  1  3  2  7  6  4  5]...
%          8+[ 0  1  3  2  7  6  4  5]... 
%          24+[ 0  1  3  2  7  6  4  5]...
%          16+[ 0  1  3  2  7  6  4  5]...
%          56+[ 0  1  3  2  7  6  4  5]...
%          48+[ 0  1  3  2  7  6  4  5]...
%          32+[ 0  1  3  2  7  6  4  5]...
%          40+[ 0  1  3  2  7  6  4  5]]+1);
%    
%    inp=reshape(bits_in,6,full_len/6);
%    
%    mod_symbols=table([32 16 8 4 2 1]*inp+1);  % maps transmitted bits into 64QAM symbol
else
   error('Unimplemented deciding');
end


