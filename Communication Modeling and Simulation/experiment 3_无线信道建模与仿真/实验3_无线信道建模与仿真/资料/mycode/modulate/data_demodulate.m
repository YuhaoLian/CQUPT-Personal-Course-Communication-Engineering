

function [soft_bits_out]  = data_demodulate(rx_symbols, sim_options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BPSK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(findstr(sim_options, 'BPSK'))
    soft_bits = (sign(rx_symbols)+1)'/2;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% QPSK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif ~isempty(findstr(sim_options, 'QPSK'))
     ModType=2;
     DemodData = zeros(1,length(rx_symbols)*ModType);
	 DemodData(1:2:end) = real(rx_symbols);
	 DemodData(2:2:end) = imag(rx_symbols);
     soft_bits = (sign(DemodData)+1)'/2;
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 8PSK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif ~isempty(findstr(sim_options, '8PSK'))
     table=[0 0 0 exp(j*pi*1/8) 0 0 1 exp(j*pi*3/8) 0 1 1 exp(j*pi*5/8) 0 1 0 exp(j*pi*7/8) 1 1 0 exp(-j*pi*7/8) 1 1 1 exp(-j*pi*5/8) 1 0 1 exp(-j*pi*3/8) 1 0 0 exp(-j*pi*1/8)];
     table2=reshape(table,4,8);
     for ii=1:1:length(rx_symbols)
         for jj=1:8
             if rx_symbols(1,ii)==table2(4,jj)
                soft_bits(1,(ii-1)*3+1)=table2(3,jj);    
                soft_bits(1,(ii-1)*3+2)=table2(2,jj);  
                soft_bits(1,(ii-1)*3+3)=table2(1,jj);  
             end
         end
     end
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 16QAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ~isempty(findstr(sim_options, '16QAM'))	
   %soft_bits = rx_qam16_demod(rx_symbols);
   soft_bits = zeros(4*size(rx_symbols,1), size(rx_symbols,2));  % Each symbol consists of 4 bits

   bit0 = real(rx_symbols);
   bit2 = imag(rx_symbols);

   bit1 = 2/sqrt(10)-(abs(real(rx_symbols)));
   bit3 = 2/sqrt(10)-(abs(imag(rx_symbols)));

   soft_bits(1:4:size(soft_bits,1),:) = bit0;
   soft_bits(2:4:size(soft_bits,1),:) = bit1;
   soft_bits(3:4:size(soft_bits,1),:) = bit2;
   soft_bits(4:4:size(soft_bits,1),:) = bit3;
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 64QAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ~isempty(findstr(sim_options, '64QAM'))
   %soft_bits = rx_qam64_demod(rx_symbols);
   soft_bits = zeros(6*size(rx_symbols,1), size(rx_symbols,2));  % Each symbol consists of 6 bits
   bit0 = real(rx_symbols);
   bit3 = imag(rx_symbols);

   bit1 = 4/sqrt(42)-abs(real(rx_symbols));
   bit4 = 4/sqrt(42)-abs(imag(rx_symbols));


   for m=1:size(rx_symbols,2)
      for k=1:size(rx_symbols,1)
         if abs(4/sqrt(42)-abs(real(rx_symbols(k,m)))) <= 2/sqrt(42)  % bit is one
            bit2(k,m) = 2/sqrt(42) - abs(4/sqrt(42)-abs(real(rx_symbols(k,m))));
         elseif abs(real(rx_symbols(k,m))) <= 2/sqrt(42) % bit is zero, close to real axis
            bit2(k,m) = -2/sqrt(42) + abs(real(rx_symbols(k,m)));
         else
            bit2(k,m) = 6/sqrt(42)-abs(real(rx_symbols(k,m))); % bit is zero 
         end;
      
         if abs(4/sqrt(42)-abs(imag(rx_symbols(k,m)))) <= 2/sqrt(42)  % bit is one
            bit5(k,m) = 2/sqrt(42) - abs(4/sqrt(42)-abs(imag(rx_symbols(k,m))));
         elseif abs(imag(rx_symbols(k,m))) <= 2/sqrt(42) % bit is zero, close to real axis
            bit5(k,m) = -2/sqrt(42) + abs(imag(rx_symbols(k,m)));
         else
            bit5(k,m) = 6/sqrt(42)-abs(imag(rx_symbols(k,m)));
         end;
      end;
   end;

   soft_bits(1:6:size(soft_bits,1),:) = bit0;
   soft_bits(2:6:size(soft_bits,1),:) = bit1;
   soft_bits(3:6:size(soft_bits,1),:) = bit2;
   soft_bits(4:6:size(soft_bits,1),:) = bit3;
   soft_bits(5:6:size(soft_bits,1),:) = bit4;
   soft_bits(6:6:size(soft_bits,1),:) = bit5;
   
else
   error('Undefined modulation');
end

soft_bits_out = soft_bits(:)';

