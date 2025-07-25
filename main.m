% ---------------------------------------------
% ----- INFORMATIONS -----
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Date            : 2025-07-22
%
% ----- Main idea -----
%   Simulating TX, PROPAGATION and RX side of a (almost) complete optical
%   telecommunication system
%
%   - transmitter side
%   - propagation
%   - coherent detection - partial
%   - power profile estimation
%   - anomaly detection
%
% ----- BIBLIOGRAPHY -----
%   Functions           : Optilux 2009 and 2022 versions
%   Author              : Paolo SERENA
%   Author contact      : serena@tlc.unipr.it
%   Date                : 2009 and 2022
%   Title of program    : Optilux
%   Code version        : 2009 and 2022
%   Type                : Optical simulator toolbox - source code
%   Web Address         : v2009 --- Partage Zimbra "Phd Louis Tomczyk"/
%                               Optilux/Optilux v2009
%                         v2022 --- https://optilux.sourceforge.io/
% -----------------------
%   Articles
%   Authors             : Takahito TANIMURA - Setsuo YOSHIDA -
%                         Kazuyuki TAJIMA - Shoichiro ODA -
%                         Takeshi HOSHIDA
%   Title               : Fiber-Longitudinal Anomaly Position Identifica
%                         tion Over Multi-Span Transmission Link Out of
%                         Receiver-end Signals
%   Jounal              : IEEE - Journal of Ligthwave Technology
%   Volume - N°         : 38-9
%   Date                : 2020-05-09
%   DOI                 : 10.1109/JLT.2020.2984270
%
%   Authors             : Alix MAY - Fabien BOITIER - Elie AWWAD -
%                         Petros RAMANTANIS - Matteo LONARDI -
%                         Philippe CIBLAT
%   Title               : Receiver-Based Experimental Estimation of Power
%                         Losses in Optical Networks
%   Jounal              : IEEE - Photonics Technology Letters
%   Volume - N°         : 33-22
%   Date                : 2021-11-22
%   DOI                 : 10.1109/LPT.2021.3115627
% ---------------------------------------------

%% MAINTENANCE
clear
% where = {"server"|"local"}
where   = "server";
% what 	= {PPE|nli - quick test|data - calibration|estimation - pd|correction - help apo| - help loc|}
what    = ["PPE","quick test","estimation","pd","help apo","help loc"];

LossySpan       = 3; 	%input(' ------- Enter in what span there is the loss\n');
LossLocInSpan   = 20;	%input(' ------- Enter distance\n'); % [km]
LossValues      = 3;	%linspace(1,2,2);% [dB]

DISPs = 20;

for k = 1:length(DISPs)

    clear("Axis","las"",TX","amp","PPEparams","ppmon")
    init_step(where,what)

    %% PARAMETERS
    % -------------------------------------------------------------------------
    % Master ones
    % -------------------------------------------------------------------------
    %%% Modulation parameters
    tx.Nch                  = 1;            % []        number of WDM channels
    tx.PdBm                 = 5;            % [dBm]     power at fibre input

    %%% Axis parameters
    Axis.Nsymb              = 2^10;         % []        number of symbols
    Axis.symbrate           = 32;           % [Gbd]     symbol rate

    %%% Laser parameters
    las.n_polar             = 1;            % []        number of polarisations

    %%% Optical amplifiers
    amp.type                = 'classic';    % choose among noAWGN-ideal
    amp.Nspan               = LossySpan+1;            % []        number of spans;

    %%% Power Profile Estimator
    if is_string_in_stringArray(what,"quick test") == 1
        Nrea                = 10;
        Ntries              = 10;
    elseif is_string_in_stringArray(what,"data") == 1
        if strcmp(where,"server") == 1
            Nrea            = 20;
            Ntries          = 1000;
        else
            Nrea            = 1;
            Ntries          = 100;
        end
    end

    PPEparams.repet.Ntries  = Ntries;
    PPEparams.anomalies.get_ai = 1;

    if is_string_in_stringArray(what,'pd') == 1
        tx.pd               = 10;
    end

    %%% Anomalies detection
    if is_string_in_stringArray(what,"calibration") == 1
        PPEparams.anomalies.ai.calibration.getCal       = 1;
    elseif is_string_in_stringArray(what,"estimation") == 1
        PPEparams.anomalies.ai.calibration.getCal       = 0;
    end

    if is_string_in_stringArray(what,'help apo') == 1
        PPEparams.anomalies.ai.apodisation.helpApo     = 1;
    else
        PPEparams.anomalies.ai.apodisation.helpApo     = 0;
    end

    % -------------------------------------------------------------------------
    % Slave ones
    % -------------------------------------------------------------------------

    tx          = set_tx(tx);
    las         = set_las(tx,las);
    Axis        = set_axis(tx,Axis);

    if exist('ft','var')
        ft      = set_ft(las,ft);
    else
        ft      = set_ft(las);
    end
    ft.disp      = DISPs(k);

    SPANS       = LossySpan;
    LOSSES      = LossValues;         % [dB]
    DISTS       = LossLocInSpan;      % [km]

    PPEparams.anomalies.ai.correction.getAFC = 1;

    if is_string_in_stringArray(what,'help loc') == 1
        PPEparams.anomalies.ai.estimation.zloss.helpLocVal  = (SPANS-1)*ft.length*1e-3+DISTS;
    else
        PPEparams.anomalies.ai.estimation.zloss.helpLocVal = [];
    end

    % -------------------------------------------------------------------------
    % Simulation
    % -------------------------------------------------------------------------

    % initalisation
    inigstate(Axis.Nsamp,Axis.fs)

    for LL = 1:length(LOSSES)
        sprintf(' ========== DISP = %i [ps/nm/km]',ft.disp)
        %f(x,y) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2
        x = LL; % loss [dB]
        y = ft.length*1e-3; % Lspan [km]

        p00 =   0.003647;% (0.003526, 0.003768)
        p10 =   1.908e-05;% (7.972e-06, 3.019e-05)
        p01 =   0.0001431;% (0.00014, 0.0001462)
        p20 =  -1.336e-06;% (-2.037e-06, -6.344e-07)
        p11 =  -3.306e-07;% (-4.345e-07, -2.267e-07)
        p02 =  -6.145e-07;% (-6.35e-07, -5.941e-07)

        tmp = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2;

        if is_string_in_stringArray(what,"estimation") == 1
            if is_string_in_stringArray(what,"pd") == 1
                PPEparams.anomalies.ai.calibration.calC = 1.1105e-2/tx.Plin;
                PPEparams.anomalies.ai.calibration.calF = 1.1105e-2;
            elseif is_string_in_stringArray(what,"correction") == 1
                PPEparams.anomalies.ai.calibration.calC = tmp/tx.Plin;
                PPEparams.anomalies.ai.calibration.calF = tmp;
            end
        end

        if is_string_in_stringArray(what,'calibration') == 1
            PPEparams.anomalies.ai.calibration.lossdB   = LOSSES(LL);   % [dB]
        end

        for DD = 1:length(DISTS)
            sprintf(' ---------- DISTS = %i [km]',DISTS(DD))
            if is_string_in_stringArray(what,'calibration') == 1
                PPEparams.anomalies.ai.calibration.zloss = (SPANS-1)*ft.length*1e-3+DISTS(DD); % [km]
            end

            losses              = [SPANS,DISTS(DD),LOSSES(LL)];
            amp                 = set_topology(tx,ft,amp,losses);
            if is_string_in_stringArray(what,'correction') == 1
                PPEparams       = set_PPEparams(tx,ft,amp,PPEparams,Axis);
            else
                PPEparams       = set_PPEparams(tx,ft,amp,PPEparams);
            end

            % enabling parfoor outputs
            if strcmp(PPEparams.anomalies.method,'pd') == 1
                PPEparamsOUT    = cell(1,Nrea);
                PP_REF_OUT  = cell(1,Nrea);
                PP_S_OUT  = cell(1,Nrea);
                AI_OUT          = cell(1,Nrea);
            else
                PPEparamsOUT    = cell(1,Nrea);
                PP_REF_cSc_OUT  = cell(1,Nrea);
                PP_MON_cSc_OUT  = cell(1,Nrea);
                AI_OUT          = cell(1,Nrea);
            end

            for k = 1:Nrea
                sprintf('---------- %i-th iteration',k)

                %%% getting the reference and monitored profiles
                ppref           = get_ppr(PPEparams,tx,ft,amp,Axis,las);
                ppmon           = get_pp(PPEparams,Axis,las,tx,ft,amp);

                %%% estimation of the anomaly - with the DPD
                if strcmp(PPEparams.anomalies.method,'pd') == 1
                    % [PPEparamsOUT{k},ppref,pps,AI] = get_anomalies( ...
                    %                             PPEparams,ppref,ppmon,amp,ft);
                    [PPEparamsOUT{k},PP_REF_OUT{k},PP_S_OUT{k},AI_OUT{k}] = get_anomalies( ...
                                                PPEparams,ppref,ppmon,amp,ft);

                %%% estimation of the anomaly - with the correction
                else
                    % [PPEparamsOUT{k},ppref_cSc,ppmon_cSc,AI] = get_anomalies( ...
                    %                         PPEparams,ppref,ppmon,amp,ft,Axis,las,tx);
                    [PPEparamsOUT{k},PP_REF_cSc_OUT{k},PP_MON_cSc_OUT{k},AI_OUT{k}] = get_anomalies( ...
                                            PPEparams,ppref,ppmon,amp,ft,Axis,las,tx);
                end

            end

        end
    end
end

Ncols = zeros(1,Nrea);
for k = 1:Nrea
    Ncols(k) = PPEparamsOUT{k}.link.nsteps_tot;
end
NCOLS = min(Ncols);

pp_ref_out  = zeros(Nrea, NCOLS);
pp_mon_out  = zeros(Nrea, NCOLS);
ai_out      = zeros(Nrea, NCOLS);

for k = 1:Nrea
    if strcmp(PPEparams.anomalies.method,'pd') == 1
        pp_ref_out(k,:) = PP_REF_OUT{k}(end-NCOLS+1:end);
        pp_mon_out(k,:) = PP_S_OUT{k}(end-NCOLS+1:end);
    else
        pp_ref_out(k,:) = PP_REF_cSc_OUT{k}(end-NCOLS+1:end);
        pp_mon_out(k,:) = PP_MON_cSc_OUT{k}(end-NCOLS+1:end);
    end
    ai_out(k,:)         = AI_OUT{k}(end-NCOLS+1:end);
end

pp_ref_mean = mean(pp_ref_out);
pp_mon_mean = mean(pp_mon_out);
ai_mean         = mean(ai_out);

figure
    subplot(1,2,1)
        hold on
        plot(pp_ref_mean)
        plot(pp_mon_mean)
        legend("ref","mon")
    subplot(1,2,2)
        plot(ai_mean)

% clear all
% close all
% exit
