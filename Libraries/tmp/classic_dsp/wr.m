function [patbinhat,Ewr] = wr(E,pat,patbin,Axis,las,tx,rx)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : WR - WAVE RECONSTRUCTION
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-06-03
%   Version         : 1.0
%
% ----- Main idea -----
%   Reconstruction of the signal from the received signal
%
% ----- INPUTS -----
%   E:      [structure] containing the time/frequency axis informations
%               - LAMBDA    [scalar] [nm]       wavelength
%               - FIELD     [array] [sqrt(mW)]  normalised electric fields
%   PAT:    [array] symbols sent
%   PATBIN: [array] encoded sent symbols
%   AXIS:   [structure] containing the time/frequency axis informations
%               - NSYMB         []      number of symbols transmitted
%               - SYMBRATE      [Gbaud] symbol rate
%   LAS:    [structure] containing the laser parameters
%               - N_POLAR       []          number of polarisations
%               - PLIN          [mW]        power of the signal
%               - LAM           [nm]        wavelength
%               - LINEWIDTH     [GHz]       
%               - n0            [dB/GHz]    noise power spectral density
%   TX:     [structure] containing the transmitter parameters
%               - MODFOR        []          modulation format
%               - PULSE_SHAPE   []
%   RX:     [structure] containing the receiver parameters
%               - EFTYPE        [string]    electrical filter (LPF) type (see MYFILTER).
%               - EBW           []          LPF bandwidth normalized to SYMBRATE.
%               - EPAR          []          electrical filter extra parameters (optional)
%               - SYNC.TYPE     [string]    time recovery method. 'da': data-aided
%               - SYNC.INTERP   [string]    interpolation method to recovery time delay
%               - MODFORMAT     [string]    modulation format, e.g., '16qam','qpsk',etc.
%
% ----- OUTPUTS -----
%  PATBINHAT:   [array] estimated sent symbols
%  EWR:         [structure] containing the laser field(s)
%               - LAMBDA    [scalar] [nm]       wavelength
%               - FIELD     [array] [sqrt(mW)]  normalised electric fields
%
% ----- BIBLIOGRAPHY -----
%   Functions : OPTILUX v2021
%   Author              : Paolo SERENA
%   Author contact      : serena@tlc.unipr.it
%   Date                : 2021
%   Title of program    : Optilux
%   Code version        : 2021
%   Type                : Optical simulator toolbox - source code
%   Web Address         : https://optilux.sourceforge.io/
% ---------------------------------------------
    % n_polar = 2;
    if size(E.field,2) == 2*size(E.lambda,1)
 
        % symbol estimation
        akhat      = rxdsp(E.field,Axis.symbrate,patbin,rx);
        patbinhat  = samp2pat(akhat,tx.modfor,rx);
        
        clear('E')

        % Ideal Electric field
        E       = lasersource(las.Plin,las.lam)

        % split in two orthogonal polarizations
        [Ex,Ey] = pbs(E);

        % construction of electric modulating signal for modulator with the
        % estimated pattern
        [elecx, normx]  = digitalmod(patbinhat(:,1:2),tx.modfor,Axis.symbrate,...
                                    tx.pulse_shape,tx);
        [elecy, normy]  = digitalmod(patbinhat(:,3:4),tx.modfor,Axis.symbrate,...
                                    tx.pulse_shape,tx);
        
        % modulation of the optical fields
        Ewrx = iqmodulator(Ex, elecx,struct('norm',normx));
        Ewry = iqmodulator(Ey, elecy,struct('norm',normy));

        % combine the two polarizations creating a PDM signal
        Ewr  = pbc(Ewrx,Ewry);

    else
        akhat       = rxdsp(E.field,Axis.symbrate,patbin,rx);
        patbinhat   = samp2pat(akhat,tx.modfor,rx);
        E           = lasersource(las.Plin,las.lam,struct(...
                        'pol','single',...
                        'linewidth',10^(-Inf),...
                        'n0',-Inf));

        [elec,norm] = digitalmod(patbinhat,tx.modfor,Axis.symbrate,...
                                tx.pulse_shape,tx);
        Ewr         = iqmodulator(E,elec,struct('norm',norm));
    end

%     disp("  WR finished")
end