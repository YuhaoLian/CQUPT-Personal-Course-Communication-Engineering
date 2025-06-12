% File: channel_maskcallback3.m
% Matlab Simulink IEEE 802.15.3a channel model by Tim Becker
% San Jose State University, Winter 2004
% Advisor: Robert Morelos-Zaragoza

chantype = get_param(gcb,'chantype');
T = get_param(gcb,'T');
nT = get_param(gcb,'nT');
SPF = get_param(gcb,'SPF');

switch chantype
    case 'CM1'
        LOSflag = 'on';
        LAMBDA =   0.0233; lambda = 2.5;
        GAMMA = 7.1; gamma = 4.3;
        sigma1 = 3.3941; sigma2 = 3.3941; std_shadow = 3;
        set_param(gcb,'MaskValues',{T,nT,SPF,chantype,'on','0.0233','2.5','7.1','4.3','3.3941','3.3941','3'});
        set_param(gcb,'MaskEnables',{'on','on','on','on','off','off','off','off','off','off','off','off'});
    case 'CM2'
        LOSflag = 'off';
        LAMBDA = 0.4; lambda = 0.5;
        GAMMA = 5.5; gamma = 6.7;
        sigma1 = 3.3941; sigma2 = 3.3941; std_shadow = 3;
        set_param(gcb,'MaskValues',{T,nT,SPF,chantype,'off','0.4','0.5','5.5','6.7','3.3941','3.3941','3'});
        set_param(gcb,'MaskEnables',{'on','on','on','on','off','off','off','off','off','off','off','off'});
    case 'CM3'
        LOSflag = 'off';
        LAMBDA = 0.0667; lambda = 2.1;
        GAMMA = 14; gamma = 7.9;
        sigma1 = 3.3941; sigma2 = 3.3941; std_shadow = 3;
        set_param(gcb,'MaskValues',{T,nT,SPF,chantype,'off','0.0667','2.1','14','7.9','3.3941','3.3941','3'});
        set_param(gcb,'MaskEnables',{'on','on','on','on','off','off','off','off','off','off','off','off'});
    case 'CM4'
        LOSflag = 'off';
        LAMBDA = 0.0667; lambda = 2.1;
        GAMMA = 24; gamma = 12;
        sigma1 = 3.3941; sigma2 = 3.3941; std_shadow = 3;
        set_param(gcb,'MaskValues',{T,nT,SPF,chantype,'off','0.0667','2.1','24','12','3.3941','3.3941','3'});
        set_param(gcb,'MaskEnables',{'on','on','on','on','off','off','off','off','off','off','off','off'});
    otherwise      % enter custom parameters
        set_param(gcb,'MaskEnables',{'on','on','on','on','on','on','on','on','on','on','on','on'});
end