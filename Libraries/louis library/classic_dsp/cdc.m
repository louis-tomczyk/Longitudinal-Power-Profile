function Ecdc = cdc(Ein,ft,RX)

    % ---------------------------------------------
    % ----- INFORMATIONS -----
    %   Function name   : CDC - CHROMATIC DISPERSION COMPENSATION
    %   Author          : louis tomczyk
    %   Institution     : Telecom Paris
    %   Email           : louis.tomczyk.work@gmail.com
    %   Date            : 2023-01-06
    %   Version         : 1.7
    %
    % ----- Main idea -----
    %   Chromatic compensation using ideal Dispersion Compensating Fiber (DCF).
    %   The chromatically dispersed field is propagated into a DCF having the
    %   exact opposite dispersion characteristics (D [ps/nm/km], S [ps/nm²/km])
    %   but without any (non) linear effects affecting the propagation.
    %
    % ----- INPUTS -----
    %   EIN:    [structure] containing the Fields to be compensated
    %           structure should be organised as:
    %               - LAMBDA [nm]: wavelength
    %               - FIELD [sqrt(mW)]: normalised electric fields
    %   FT:     [structure] fiber parameters
    %           structure should at least contain:
    %               - ALPHADB   [dB/km]         Power attenuation
    %               - DISP      [ps/nm/km]      Dispersion
    %               - SLOPE     [ps/nm²/km]     Slope of the dispersion
    %               - LENGTH    [m]             Length
    %               - n2        [W²/m]          Non linear index
    %               - pmdpar    [ps/sqrt(km)]   Polarization Mode Dispersion
    %   DL:     [scalar] [m] length of the fiber you want CDC
    %
    % ----- OUTPUTS -----
    %  ECDC:    [structure] Chromatic Dispersion Compensated field
    %           structure containing
    %               - LAMBDA [scalar] [nm]: wavelength
    %               - FIELD [array] [sqrt(mW)]: normalised electric fields
    %
    % ----- BIBLIOGRAPHY -----
    %   Functions   : FIBER - PBC - PBS
    %   Author              : Paolo SERENA
    %   Author contact      : serena@tlc.unipr.it
    %   Date                : 2021
    %   Title of program    : Optilux
    %   Code version        : 2021
    %   Type                : Optical simulator toolbox - source code
    %   Web Address         : https://optilux.sourceforge.io/
    % ---------------------------------------------
    
    if isstruct(RX) == 0
        rx.CDC.len_tot = RX;
    else
        rx = RX;
    end

    % if no dispersion or no length, then return the input
    if ft.disp == ft.slope && ft.disp == 0 || rx.CDC.len_tot == 0
        
        Ecdc = Ein;
        
    else

        % DCF properties
        fc.type     = "ideal DCF";
        fc.length   = rx.CDC.len_tot;
        fc.lambda   = ft.lambda;

        fc.alphadB  = 0;
        fc.alphaLin = 0;

        fc.disp     = -ft.disp;
        fc.slope    = -ft.slope;
        fc.beta2    = -ft.beta2;
        fc.beta3    = -ft.beta3;
        
        fc.coupling = 'none';
        fc.pmdpar   = 0;
        fc.nplates  = 1;
        fc.ismanakov= false;
        
        fc.n2       = 0;
        fc.gf       = 0;
        fc.aeff     = ft.aeff;

        % n_polar = 2;
        if size(Ein.field,2) == 2*size(Ein.lambda,1)

            % split the polarisations
            [Einx,Einy]= sep_XYfields(Ein);

            % compensate each polarisation
            Ecdcx  = fiber(Einx,fc);
            Ecdcy  = fiber(Einy,fc);
    
            % recombine both polarisations
            Ecdc    = merge_XYfields(Ecdcx,Ecdcy);
            
        % n_polar = 1;
        else
            Ecdc   = fiber(Ein,fc);
        end

    end