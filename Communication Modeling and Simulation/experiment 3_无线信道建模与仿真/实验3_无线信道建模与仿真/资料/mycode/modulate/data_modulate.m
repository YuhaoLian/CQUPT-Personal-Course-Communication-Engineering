
function mod_symbols = data_modulate(bits_in, modulation)

full_len = length(bits_in);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BPSK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(findstr(modulation, 'BPSK'))
   table=[-1 1];
   inp=bits_in;
   mod_symbols=table(inp+1);
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% QPSK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ~isempty(findstr(modulation, 'QPSK'))
   table=exp(j*[-3/4*pi 3/4*pi 1/4*pi -1/4*pi]);  % generates QPSK symbols
   table=table([0 1 3 2]+1); % Gray code mapping pattern for QPSK symbols
   inp=reshape(bits_in,2,full_len/2);
   mod_symbols=table([2 1]*inp+1);  % maps transmitted bits into QPSK symbols
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 8PSK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ~isempty(findstr(modulation, '8PSK'))
   table=exp(j*[1/8*pi 3/8*pi 7/8*pi 5/8*pi -1/8*pi -3/8*pi -7/8*pi -5/8*pi]);  % generates 8PSK symbols
   inp=reshape(bits_in,3,full_len/3);
   mod_symbols=table([1 2 4]*inp+1);  % maps transmitted bits into 8PSK symbols
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 16QAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ~isempty(findstr(modulation, '16QAM'))
   % generates 16QAM symbols
   m=1;
   for k=-3:2:3
      for l=-3:2:3
         table(m) = (k+j*l)/sqrt(10); % power normalization
         m=m+1;
      end;
   end;
   table=table([0 1 3 2 4 5 7 6 12 13 15 14 8 9 11 10]+1); % Gray code mapping pattern for 8-PSK symbols
   inp=reshape(bits_in,4,full_len/4);
   mod_symbols=table([8 4 2 1]*inp+1);  % maps transmitted bits into 16QAM symbols
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 64QAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ~isempty(findstr(modulation, '64QAM'))
   % generates 64QAM symbols
   m=1;
   for k=-7:2:7
      for l=-7:2:7
         table(m) = (k+j*l)/sqrt(42); % power normalization
         m=m+1;
      end;
   end;
   table=table([[ 0  1  3  2  7  6  4  5]...
         8+[ 0  1  3  2  7  6  4  5]... 
         24+[ 0  1  3  2  7  6  4  5]...
         16+[ 0  1  3  2  7  6  4  5]...
         56+[ 0  1  3  2  7  6  4  5]...
         48+[ 0  1  3  2  7  6  4  5]...
         32+[ 0  1  3  2  7  6  4  5]...
         40+[ 0  1  3  2  7  6  4  5]]+1);
   
   inp=reshape(bits_in,6,full_len/6);
   
   mod_symbols=table([32 16 8 4 2 1]*inp+1);  % maps transmitted bits into 64QAM symbol
else
   error('Unimplemented modulation');
end


