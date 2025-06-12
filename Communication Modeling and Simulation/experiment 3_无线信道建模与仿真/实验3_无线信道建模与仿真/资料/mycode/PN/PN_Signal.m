function [PN_Complex]=PNSignal(Input_signal_number)

PN_RegStateI = [zeros(1,14) 1];  
PN_RegStateQ = [zeros(1,14) 1];  
PN_PolynI = [0 1 0 0 0 1 1 1 0 1 0 0 0 0 1];
PN_PolynQ = [0 0 1 1 1 0 0 0 1 1 1 1 0 0 1];

PN_SequenceI = zeros(1,Input_signal_number);
PN_SequenceQ = zeros(1,Input_signal_number);

[PN_SequenceI,PN_RegStateI] = PN_Generator(PN_PolynI,PN_RegStateI,Input_signal_number);
[PN_SequenceQ,PN_RegStateQ] = PN_Generator(PN_PolynQ,PN_RegStateQ,Input_signal_number);

% PN_Sequence_I=PN_SequenceI;
% PN_Sequence_Q=PN_SequenceQ;

PN_Complex = (PN_SequenceI + j*PN_SequenceQ)/sqrt(2);
% for i=1:sub_numbe
%     PN_Complex_1(i+timing_offset)=PN_Complex1(i);
% end
% 
% PN_Complex_1=PN_Complex1;
%PN_Product=0;
%for i=1:sub_number
%    PN_Product=PN_Product+PN_Sequence(i)*conj(PN_Sequence1(i));
%end
%PN_Product  
