function Rx = set_rx(varargin)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : SET_DSP - Digital Signal Processing
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-10-03
%   Version         : 1.5
%
% ----- Main idea -----
%   Set the DSP structure from the laser transmitter parameters and
%   optional dsp parameters
% 
% ----- INPUTS -----
%   DSP:        (strcture) containing DSP parameters - OPTIONAL
%               - WR        (bool)[]    = 1 if decision step
%               - CPC_AVG   (scalar)[]  window averaging length
%
% ----- OUTPUTS -----
%   DSP         (structure) containing the fields:
%               - CDC       (bool)[]    Chromatic Dispersion Compensation
%               - PMDC      (bool)[]    Polarisation Mode DC
%               - CPC       (bool)[]    Carrier Phase C
%               - CPC_AVG   (scalar)[]  Window averaging length
%               - WR        (bool)[]    Waveform Reconstruction
%
% ----- BIBLIOGRAPHY -----
%   Category            : Book
%   Author              : GANG-DING Peng
%   Title               : Handbook of Optical Fibers
%   Author contact      : NA
%   Date                : 2019
%   Editor              : Springer
%   DOI                 : 10.1007/978-981-10-7087-7
%   ISBN                : 978-981-10-7085-3
%   Pages/Equations     : 176 / 50-51
%   Title of program    : NA
%   Code version        : NA
%   Type                : NA
%   Web Address         : NA
% ---------------------------------------------

    %% MAINTENANCE
    argnames = string(nargin);
    for k=1:nargin
        argnames(k)  = inputname(k);
    end

    is_arg_missing('Axis',argnames);
    is_arg_missing('dsp',argnames);
    is_arg_missing('tx',argnames);
    is_arg_missing('las',argnames);
    is_arg_missing('ft',argnames);
    is_arg_missing('amp',argnames);

    Axis= varargin{argnames == 'Axis'};
    dsp = varargin{argnames == 'dsp'};
    tx  = varargin{argnames == 'tx'};
    las = varargin{argnames == 'las'};
    ft  = varargin{argnames == 'ft'};
    amp = varargin{argnames == 'amp'};

%     %% ORIGINAL OPTILUX PARAMETERS
%     Rx.ebw          = 0.5;            % electrical filter bandwidth normalized to symbrate
%     Rx.eftype       = 'rootrc';       % electrical filter type
%     Rx.epar         = tx.rolloff;     % electrical filter extra parameters
%     Rx.modformat    = tx.modfor;      % modulation format
%     Rx.obw          = Inf;            % optical filter bandwidth normalized to symbrate
%     Rx.oftype       = 'gauss';        % optical filter type
%     Rx.sync.type    = 'da';           % time-recovery method
%     Rx.type         = 'bin';          % binary pattern

    %% ANALOG TO DIGITAL CONVERTER
    Rx.ADC.Nsps     = 2;


    %% CHROMATIC DISPERSION COMPENSATION
    length_pd       = tx.pd/ft.disp*1e3;
    length_link     = amp.Nspan*ft.length;
    Rx.CDC.len_tot  = length_pd+length_link;
    
%%% uncomment if want to downsample in CDC step
%     Rx.CDC.Nsps_b4  = Axis.Nsps;
%     Rx.CDC.Nsps_fnl = 2;

    if sum(strcmp(argnames,'rx')) == 1
        rx = varargin{argnames == 'rx'};
        if isfield(Rx.CDC,'method') == 0
            rx.CDC.method = "fibre";
        else
            rx.CDC.method = "FIR";
        end
    else
        Rx.CDC.method = "fibre";
    end


    %% CARRIER PHASE COMPENSATION
    if sum(strcmp(argnames,'rx')) == 1
        rx = varargin{argnames == 'rx'};
        if isfield(rx,'CPC') == 0
            Rx.CPC.navg = 51;
        else
            Rx.CPC.navg = rx.CPC.navg;
        end
    else
        Rx.CPC.navg = 51;
    end

    %% POLARISATION MODE DISPERSION COMPENSATION
    if isfield(dsp,'pmdc') == 1
        CMAparams.R         = [1,1];
        CMAparams.mu        = 1e-3;

        lll     = ft.length*amp.Nspan;
        Ncd     = 2*pi*abs(ft.beta2)*lll*(tx.Nbps*Axis.symbrate*1e9)^2;
        Npmd    = 2*pi*abs(ft.beta2)*sqrt(lll)*(tx.Nbps*Axis.symbrate*1e9)^2;
        Ntaps   = ceil((Ncd+Npmd)/10);
        % round up to next dozen
        Ntaps   = ceil(Ntaps/10)*10;
        if Ntaps < 20
            Ntaps = 21;
        end
        % ensure odd number or taps
        if mod(Ntaps,2) ==0
            Ntaps = Ntaps+1;
        end

        CMAparams.taps      = Ntaps;
        CMAparams.txpolars  = las.n_polar;
        CMAparams.phizero   = 0;        % rotation angle between tx and rx field
        CMAparams.eps       = 1e-3;

        CMAparams.plot      = 0;
        Rx.CMA              = CMAparams;
    end

    %% NON LINEAR NOISE ESTIMATION
    if isfield(dsp,"nlne") == 1
        Rx.NLN.plot     = 1;
        Rx.NLN.stats    = "var";
    end

end