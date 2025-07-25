function Efilt = rx_filt(Elink,Axis,tx)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : RX_FILT - Receiver side filtering
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-09-10
%   Version         : 1.0
%
% ----- Main idea -----
%   Filter the input field to select the channel using the OPT_FILT
%   Then isolate the wanted channel from others using the ELEC_FILT.
%
% ----- INPUTS -----
%   ELINK   (structure) containing the Fields to be filtered
%               - LAMBDA    [nm]        wavelength
%               - FIELD     [sqrt(mW)]  Normalised electric fields
%   AXIS    (structure)     See SET_AXIS function for details
%   TX      (structure)     See SET_TX function for details
%
% ----- BIBLIOGRAPHY -----
%   Functions           : Nyquist_WDM_DP_QPSKand16QAM_Vfinal
%   Author              : Yves JAOUEN
%   Author contact      : yves.jaouen@telecom-paris.fr
%   Date                : Unknown
%   Title of program    : Code DP_QPSKand16QAM
%   Code version        : 2022
%   Type                : Optical simulator toolbox - source code
%   Web Address         : NA
% ---------------------------------------------

    Eoptfilt    = opt_filt(Elink,Axis);         % Optical filtering
    Efilt       = elec_filt(Eoptfilt,Axis,tx);  % Electrical filtering

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%    NESTED FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----------------------------------------------------
function Efilt = opt_filt(Ein,Axis)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : OPT_FILT
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-09-09
%   Version         : 1.0
%
% ----- Main idea -----
%   Filter the input field to select the channel.
%   Note that this version only selects the center channel.
%   So if an even number of channels are multipled, you will
%   get only noise.
%
% ----- INPUTS -----
%   EIN:    (structure) containing the Fields to be filtered
%               - LAMBDA    [nm]        Wavelength
%               - FIELD     [sqrt(mW)]  Normalised electric fields
% ----- BIBLIOGRAPHY -----
%   Functions           : Nyquist_WDM_DP_QPSKand16QAM_Vfinal
%   Author              : Yves JAOUEN
%   Author contact      : yves.jaouen@telecom-paris.fr
%   Date                : Unknown
%   Title of program    : Code DP_QPSKand16QAM
%   Code version        : 2022
%   Type                : Optical simulator toolbox - source code
%   Web Address         : NA
% ---------------------------------------------

    Elink_XY        = transpose(Ein.field);
    
    m_elec          = 20;  
    Bopt            = 75;
    Hopt            = exp(-(fftshift(Axis.freq)/(Bopt/2)).^(2*m_elec));

    n_polar = size(Ein.field,2);
    if n_polar == 1
        Efilt_XY        = ifft(Hopt .* fft(Elink_XY));
    else
        Efilt_XY(1,:)   = ifft(Hopt .* fft(Elink_XY(1,:)));
        Efilt_XY(2,:)   = ifft(Hopt .* fft(Elink_XY(2,:)));
    end

    Efilt.lambda    = Ein.lambda;
    Efilt.field     = transpose(Efilt_XY);

%-----------------------------------------------------
function Eelecfilt = elec_filt(Eoptfilt,Axis,tx)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : ELEC_FILT
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-09-09
%   Version         : 1.0
%
% ----- Main idea -----
%   Filter the input optical filtered field to ensure getting only
%   the channel.
%
% ----- INPUTS -----
%   EOPTFILT:    (structure) containing the Fields to be filtered
%               - LAMBDA    [nm]        Wavelength
%               - FIELD     [sqrt(mW)]  Normalised electric fields
% ----- BIBLIOGRAPHY -----
%   Functions           : Nyquist_WDM_DP_QPSKand16QAM_Vfinal -
%                         RRCfilter
%   Author              : Yves JAOUEN
%   Author contact      : yves.jaouen@telecom-paris.fr
%   Date                : Unknown
%   Title of program    : Code DP_QPSKand16QAM
%   Code version        : 2022
%   Type                : Optical simulator toolbox - source code
%   Web Address         : NA
% ---------------------------------------------

    n_polar = size(Eoptfilt.field,2);
    if n_polar == 1

        Xin                 = transpose(Eoptfilt.field);
        Eelecfilt.lambda    = tx.lambda;
        Eelecfilt.field     = transpose(RRC_filter(Xin, Axis.symbrate, tx.rolloff, Axis.fs));
    
    else
        [Ex,Ey]     = sep_XYfields(Eoptfilt);
        Xin         = transpose(Ex.field);
        Yin         = transpose(Ey.field);
    
        Xout.lambda = tx.lambda;
        Yout.lambda = tx.lambda;
    
        Xout.field  = transpose(RRC_filter(Xin, Axis.symbrate, tx.rolloff, Axis.fs));
        Yout.field  = transpose(RRC_filter(Yin, Axis.symbrate, tx.rolloff, Axis.fs));
    
        Eelecfilt   = merge_XYfields(Xout,Yout);
    end


%-----------------------------------------------------
function f = RRC_filter(Ax, symbrate, rolloff, Bo)
    
% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : RRC_FILT
%   Author          : Yves JAOUEN
%   Institution     : Telecom Paris
%   Email           : yves.jaouen@telecom-paris.fr
%   Date            : 2021-02-01
%   Version         : 1.0
%
% ----- MAIN IDEA -----
% ----- INPUTS -----
% ----- OUTPUTS -----
% ----- BIBLIOGRAPHY -----
% -----------------------
%   Articles
%   Author              :   Junyi WANG, Chongjin XIE, Zhongqi PAN
%   Title               :   Generation of Spectrally Efficient Nyquist-WDM
%                           QPSK Signals Using Digital FIR or FDE Filters
%                           at Transmitters
%   Jounal              :   Journal of Lightwave Technology
%   Volume - N°         :   30 - 23 
%   Date                :   2012-01-23
%   DOI                 :   10.1109/JLT.2012.2226207
% ---------------------------------------------
    
    Nsamp   = length(Ax);                  
    f       = (-Nsamp/2+1:Nsamp/2)/Nsamp*Bo;
    
    [~,b] = find(f >= - (1-rolloff)*symbrate/2);   f_low_negative  = b(1,1);
    [~,b] = find(f >= - (1+rolloff)*symbrate/2);   f_high_negative = b(1,1);
    [~,b] = find(f >= (1-rolloff)*symbrate/2);     f_low_positive  = b(1,1);
    [~,b] = find(f >= (1+rolloff)*symbrate/2);     f_high_positive = b(1,1);
    
    
    H_RC = 1+cos(pi/(rolloff*symbrate).*(abs(f)-(1-rolloff)*symbrate/2));
    H_RC = H_RC/2;

    H_RC(1,1:f_high_negative)               = 0.0;
    H_RC(1,f_high_positive:end)             = 0.0;
    H_RC(1,f_low_negative:f_low_positive)   = 1.0;
    
    H_RRC   = sqrt(H_RC);
    f       = ifft(fft(Ax) .* fftshift(H_RRC));
