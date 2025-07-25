function [Ein,tx] = TX(Axis,las,tx)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : TX - TRANSMITER
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2023-01-12
%   Version         : 1.7.2
%
% ----- Main idea -----
%   Generate the modulated laser field before sending into the channel
%
% ----- INPUTS -----
%   AXIS:   (structure) containing the time/frequency axis informations
%               - NSYMB       (scalar)[]        number of symbols
%               - SYMBRATE    (scalar)[Gbaud]   symbol rate
%
%   LAS:    (structure) containing the laser parameters
%               - N_POLAR     (scalar)[]        number of polarisations
%               - PLIN        (array)[mW]       power of the signal
%               - LAM         (scalar)[nm]      wavelength
%               - LINEWIDTH   (scalar)[GHz]       
%               - n0          (scalar)[dB/GHz]  noise pow spectral density
%
%   TX:     (structure) containing the transmitter parameters
%               - MODFOR     (string)[]         modulation format
%               - PULSE_SHAPE(string)[]
%
% ----- OUTPUTS -----
%  PAT:     (array) symbols sent
%  PATBIN:  (array) encoded sent symbols
%  EIN:     (structure) containing the laser field(s)
%               - LAMBDA    (scalar)[nm]        wavelength
%               - FIELD     (array)[sqrt(mW)]   normalised electric fields
%
% ----- BIBLIOGRAPHY -----
%   Functions           : LASERSOURCE-DATAPATTERN-DIGILTALMOD-
%                         IQMODULATOR-PBC/S
%   Author              : Paolo SERENA
%   Author contact      : serena@tlc.unipr.it
%   Date                : 2021
%   Title of program    : Optilux
%   Code version        : 2021
%   Type                : Optical simulator toolbox - source code
%   Web Address         : https://optilux.sourceforge.io/
% ---------------------------------------------
    
        if tx.seed ~= 0
            rng(1);
        end
    
        switch tx.modfor
            case 'ook'
                ncol = 1;
            case 'bpsk'
                ncol = 1;
            case 'qpsk'
                ncol = 2;
            otherwise
                nn   = str2double(tx.modfor(1:3));
                if isempty(nn) == 1|| isnan(nn)
                    nn = str2double(tx.modfor(1:2));
                    if isempty(nn) == 1|| isnan(nn)
                        nn = str2double(tx.modfor(1));
                    end
                end
                ncol = log2(nn);
        end
    
        if las.n_polar == 1
            tx.patbin.x     = zeros(Axis.Nsymb,ncol);
            tx.pat.x        = zeros(Axis.Nsymb,tx.Nch);
        else
            tx.patbin.x.x1  = zeros(Axis.Nsymb,1);
            tx.patbin.x.x2  = zeros(Axis.Nsymb,1);
            tx.patbin.y.y1  = zeros(Axis.Nsymb,1);
            tx.patbin.y.y2  = zeros(Axis.Nsymb,1);
            tx.pat.x        = zeros(Axis.Nsymb,tx.Nch);
            tx.pat.y        = zeros(Axis.Nsymb,tx.Nch);
        end
    
        if las.n_polar == 1
            % 1 polar
                E   = lasersource(las.Plin,las.lam,tx.spac,tx.Nch,struct(...
                            'pol','single',...
                            'linewidth',las.linewidth,...
                            'n0',las.n0));
    
                for k=1:tx.Nch
                    [tx.pat,tx.patbin]  = datapattern(Axis.Nsymb,'rand',struct('format',tx.modfor));
                    [elec, norm]        = digitalmod(tx.patbin,tx.modfor,Axis.symbrate,tx.pulse_shape,tx);
                    E                   = iqmodulator(E,elec,struct('nch',k,'norm',norm));
                end
                
                Ein = multiplexer(E);
    %             plot(Axis.freq,20*log10(abs(fftshift(fft(Ein.field)))))
        else
            % 2 polar
                E   = lasersource(las.Plin,las.lam,tx.spac,tx.Nch, ...
                            struct('linewidth',las.linewidth,...
                            'n0',las.n0));
    
                % split in two orthogonal polarizations
                [Ex,Ey] = pbs(E);
    
                % now separately modulate the two polarizations
                for k = 1:tx.Nch
                    [patx,patbinx]  = datapattern(Axis.Nsymb,'rand',struct('format',tx.modfor));
                    [paty,patbiny]  = datapattern(Axis.Nsymb,'rand',struct('format',tx.modfor));
                    
                    [elecx, normx]  = digitalmod(patbinx,tx.modfor,Axis.symbrate,tx.pulse_shape,tx);
                    [elecy, normy]  = digitalmod(patbiny,tx.modfor,Axis.symbrate,tx.pulse_shape,tx);
        
                    Ex  = iqmodulator(Ex,elecx,struct('nch',k,'norm',normx));
                    Ey  = iqmodulator(Ey,elecy,struct('nch',k,'norm',normy));
                end
                
    %             plotfield(Ex,'--p-')
    
                % combine the two polarizations creating a PDM signal
                E               = pbc(Ex,Ey);
                Ein             = multiplexer(E);
    
                tx.patbin.x.x1  = patbinx(:,1);
                tx.patbin.x.x2  = patbinx(:,2);
                tx.patbin.y.y1  = patbiny(:,1);
                tx.patbin.y.y2  = patbiny(:,2);
                tx.pat.x        = patx;
                tx.pat.y        = paty;
        end
%     
%         ft.alphadB  = 0;        % [dB/km]       attenuation
%         ft.lambda   = las.lam;  % [nm]          wavelength
%         ft.length   = 1e3;      % [m]           length of the fibre
%         ft.slope    = 0;        % [ps/nm²/km]   slope
%         ft.pmdpar   = 0;        % [ps/sqrt(km)] 
%         ft.aeff     = 80;       % [µm²]         effective area
%         ft.gf       = 0;        % [1/W/km]      non linear param
%         ft.n2       = 0;        % [m²/W]        nonlinear index
%         ft.coupling  = "none";
%     
%         ft.disp      = tx.pd;                % [ps/nm/km]
% 
%     
%     
%         if tx.pd ~= 0
% %             sprintf('pd = %i',tx.pd)
%             assert(nargin == 4,'Missing the FT structure to apply PreDistorsion')
%             ft  = varargin{1};
%             Epd = pd(Ein,ft,tx.pd);
%         else
%             Epd = Ein;
%         end