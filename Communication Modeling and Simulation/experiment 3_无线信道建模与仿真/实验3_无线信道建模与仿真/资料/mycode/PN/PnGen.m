function [PnSeqOut, ShiftRegState] = PnGen(CharPoly,ShiftRegState, NumChips)
% Usage: [PnSeqOut, ShiftRegState] = PnGen(CharPoly,ShiftRegState, NumChips)
%
% This function acts as a shift register sequence generator
% Inputs:
%   CharPoly      = Characteristic Polynomial of PN code
%   ShiftRegState = State of the shift register before clocking. 
%                   This Stats must be kept every time
%   NumChips      = Number of chips to compute of PN sequence
%
% Outputs:
%   PnSeqOut      = Vector holding NumChips of PN sequence
%   ShiftRegState = New State of the shift Register
% 
% Example: Create a m sequence that repeats every 2^3-1 bits
%   ShiftRegState = [ 1 0 0 ];        
%   CharPoly      = [ 1 0 1 ]; %[C1 C2 C3]  
%                 = 1 + C1*X + C2*X^2 +C3*X^3 = 1+X+X^3
%
% Original author: N.S. Correal
% For academic use only
% Application Example:
%   
%    PN_RegStateI = [zeros(1,14) 1];  
%    PN_RegStateQ = [zeros(1,14) 1];  
%    PN_PolynI = [0 1 0 0 0 1 1 1 0 1 0 0 0 0 1];  
%    PN_PolynQ = [0 0 1 1 1 0 0 0 1 1 1 1 0 0 1];
%    PN_SequenceI = zeros(1,length(TxSig_Repeat));
%    PN_SequenceQ = zeros(1,length(TxSig_Repeat));
%    [PN_SequenceI,PN_RegStateI] = PnGen(PN_PolynI,PN_RegStateI,length(You_need));
%    [PN_SequenceQ,PN_RegStateQ] = PnGen(PN_PolynQ,PN_RegStateQ,length(You_need)); 
%    PN_Complex = PN_SequenceI + j*PN_SequenceQ;    % get a complex PN sequence.




 if length(ShiftRegState) ~= length(CharPoly)
   error('Characteristic Polynomial size mismatch');
 end
 
 if sum(ShiftRegState)==0
   error('All zeros state is Not allowed')
 end

PnSeqOut = zeros(1,NumChips);
OutputBit = length(ShiftRegState);  % index to ShR output

for chip = 1:NumChips
  Feedback = rem(sum(CharPoly.*ShiftRegState),2); % Feedback Bit
  PnSeqOut(chip) = ShiftRegState(OutputBit); 
  ShiftRegState = [Feedback, ShiftRegState(1:(OutputBit-1))];
end

% Turn Unipolar code to Polar format
PnSeqOut =  2*PnSeqOut-1;
  
