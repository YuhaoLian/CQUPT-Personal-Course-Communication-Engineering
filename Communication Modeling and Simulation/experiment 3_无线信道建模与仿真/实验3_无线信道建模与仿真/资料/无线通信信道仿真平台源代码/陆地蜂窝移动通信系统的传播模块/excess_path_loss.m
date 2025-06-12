function  [LHATA, uHATA]=excess_path_loss(area_id, hB, hT)
 
if area_id==1 % small/medium cities
	
    if hB == 30	
        if hT ==1
            LHATA = -9.22;
            uHATA = 15.22;
        elseif hT ==4
            LHATA = -17.02;
            uHATA = 15.22;	
        elseif hT ==7
            LHATA = -24.82;
            uHATA = 15.22;	
        elseif hT ==10
            LHATA = -32.62;
            uHATA = 15.22;
        else
            LHATA = 0;
            uHATA = 0;
        end
	elseif hB == 50 
        if hT ==1
            LHATA = -7.93;
            uHATA = 13.77;	
        elseif hT ==4
            LHATA = -15.73;
            uHATA = 13.77;	
        elseif hT ==7
            LHATA = -23.53;
            uHATA = 13.77;	
        elseif hT ==10
            LHATA = -31.33;
            uHATA = 13.77;
        else
            LHATA = 0;
            uHATA = 0;
        end
	elseif hB == 100
        if hT ==1
            LHATA = -6.17;
            uHATA = 11.80;
        elseif hT ==4
            LHATA = -13.97;
            uHATA = 11.80;
        elseif hT ==7
            LHATA = -21.77;
            uHATA = 11.80;
        elseif hT ==10
            LHATA = -29.57;
            uHATA = 11.80;
        else
            LHATA = 0;
            uHATA = 0;
        end
	elseif hB == 200
        if hT ==1
            LHATA = -4.42;
            uHATA = 9.83;	
        elseif hT ==4
            LHATA = -12.22;
            uHATA = 9.83;	
        elseif hT ==7
            LHATA = -20.02;
            uHATA = 9.83;	
        elseif hT ==10
            LHATA = -27.82;
            uHATA = 9.83;
        else
            LHATA = 0;
            uHATA = 0;
        end
	else
        LHATA = 0;
        uHATA = 0;
	end
    
elseif area_id==2 %  large cities   
	
    if hB == 30	
        if hT ==1
			LHATA = -9.19;
			uHATA = 15.22;
        elseif hT ==4
			LHATA = -14.48;
			uHATA = 15.22;
        elseif hT ==7
			LHATA = -17.27;
			uHATA = 15.22;
        elseif hT ==10
			LHATA = -19.24;
			uHATA = 15.22;
        else
            LHATA = 0;
            uHATA = 0;
        end
	elseif hB == 50 
        if hT ==1
			LHATA = -7.90;
			uHATA = 13.77;
        elseif hT ==4
			LHATA = -13.18;
			uHATA = 13.77;
        elseif hT ==7
			LHATA = -15.97;
			uHATA = 13.77;
        elseif hT ==10
			LHATA = -17.95;
			uHATA = 13.77;
        else
            LHATA = 0;
            uHATA = 0;
        end
	elseif hB == 100
        if hT ==1
			LHATA = -6.15;
			uHATA = 11.80;
        elseif hT ==4
			LHATA = -11.43;
			uHATA = 11.80;
        elseif hT ==7
			LHATA = -14.22;
			uHATA = 11.80;
        elseif hT ==10
			LHATA = -16.19;
			uHATA = 11.80;
        else
            LHATA = 0;
            uHATA = 0;
        end
	elseif hB == 200
        if hT ==1
			LHATA = -4.39;
			uHATA = 9.83;	
        elseif hT ==4
			LHATA = -9.67;
			uHATA = 9.83;		
        elseif hT ==7
			LHATA = -12.46;
			uHATA = 9.83;	
        elseif hT ==10
			LHATA = -14.44;
			uHATA = 9.83;
        else
            LHATA = 0;
            uHATA = 0;
        end
	else
        LHATA = 0;
        uHATA = 0;
	end

elseif area_id==3 %  suburban
    
    if hB == 30	
        LHATA = -20.72;
		uHATA = 15.22;
	elseif hB == 50
		LHATA = -19.43;
		uHATA = 13.77;
	elseif hB == 100
		LHATA = -17.67;
		uHATA = 11.80;
	elseif hB == 200
		LHATA = -15.92;
		uHATA = 9.83;
    else
        LHATA = 0;
        uHATA = 0;
	end
    
elseif area_id==4 %   open areas
    
    if hB == 30	
        LHATA = -39.47;
		uHATA = 15.22;
	elseif hB == 50
		LHATA = -38.18;
		uHATA = 13.77;
	elseif hB == 100
		LHATA = -36.42;
		uHATA = 11.80;
	elseif hB == 200
		LHATA = -34.67;
		uHATA = 9.83;
    else
        LHATA = 0;
        uHATA = 0;
	end
    
else
    LHATA = 0;
    uHATA = 0;
end    
    