function [PN_Complex_1]=pn_user1(simulate_number,timing_offset)
PN_RegStateI = [zeros(1,14) 1];  
PN_RegStateQ = [zeros(1,14) 1];  
PN_PolynI = [0 1 0 0 0 1 1 1 0 1 0 0 0 0 1];
%PN_PolynI=[0 0 1 0 1];
PN_PolynQ = [0 0 1 1 1 0 0 0 1 1 1 1 0 0 1];
TxSig_Repeat(1:simulate_number)=1;
PN_SequenceI = zeros(1,length(TxSig_Repeat));
PN_SequenceQ = zeros(1,length(TxSig_Repeat));
[PN_SequenceI,PN_RegStateI] = PnGen(PN_PolynI,PN_RegStateI,length(TxSig_Repeat));
[PN_SequenceQ,PN_RegStateQ] = PnGen(PN_PolynQ,PN_RegStateQ,length(TxSig_Repeat));
%PN_Complex = PN_SequenceI + j*PN_SequenceQ;
PN_Sequence_I=PN_SequenceI;
PN_Sequence_Q=PN_SequenceQ;

PN_Complex1 = PN_Sequence_I + j*PN_Sequence_Q;
for i=1:sub_number-timing_offset
    PN_Complex_1(i+timing_offset)=PN_Complex1(i);
end
if timing_offset==0
   PN_Complex_1=PN_Complex1;
else
   PN_Complex_1(1:timing_offset)=1;
end
%PN_Product=0;
%for i=1:sub_number
%    PN_Product=PN_Product+PN_Sequence(i)*conj(PN_Sequence1(i));
%end
%PN_Product  
