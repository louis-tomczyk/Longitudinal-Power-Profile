function Epmdc = pmdc(Ecdc,Axis,rx)
    
% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : PMDC - POLARISATION MODE DISPERSION COMPENSATION
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-10-03
%   Version         : 1.3
%
% ----- Main idea -----
%   Polarisation Mode Dispersion compensation using CMA algorithm.
%   This code comes, 99.99% of it, from SOPrecovery_FSE_CMAandRDE.m,
%   see the ref
%   I only take a small of part of it to put it in a new function
%
% ----- INPUTS -----
%  ECDC    (structure) of the Chromatic Dispersion Compensated field
%               - LAMBDA    (scalar)[nm]        Wavelength
%               - FIELD     (array)[sqrt(mW)]   Normalised electric fields
%   Axis    (structure) axis parameters, see SET_AXIS
%   RX      (structure) recever parameters, see SET_RX
%
% ----- OUTPUTS -----
%  EPMDC    (structure) of the Chromatic Dispersion Compensated field
%               - LAMBDA    (scalar)[nm]        Wavelength
%               - FIELD     (array)[sqrt(mW)]   Normalised electric fields
%
% ----- BIBLIOGRAPHY -----
%   Functions           : SOPRECOVERY_FSE_CMAANDDRE
%   Author              : Yves JAOUEN
%   Author contact      : yves.jaouen@telecom-paris.fr
%   Date                : 2021-02-09
%   Title of program    : NA
%   Code version        : 2021
%   Type                : Optical simulator toolbox - source code
%   Web Address         : NA
% ---------------------------------------------

    Input   = transpose(Ecdc.field);
    Ni      = Axis.Nsps/(Axis.Nsamp/size(Input,2));
    assert(Ni == 2,'---  there should be 2 samples/symbol at this step.')
    Ntaps   = rx.CMA.taps;

    Nsamples    = length(Input(1,:));
    Input(1,:)  = Input(1,:) / mean(abs(Input(1,:)));
    Input(2,:)  = Input(2,:) / mean(abs(Input(2,:)));

    D           = floor(Ntaps/2);   % Retard de restitution    

    % Calcul des coefficients du filtre MIMO
    % ***************************************
    Naffichage = 1000;

    if length(Input) > 20000
        Napprentissage = ceil((Nsamples-Ntaps)/2)-1;  % Nombre d'�chantillons dans la phase apprentissage
    else
        Napprentissage = floor(length(Input)/2) - Naffichage - Ni;
    end

    hxx         = zeros(1,Ntaps);
    hxx(1,D+1)  = 1.0 ;
    hyy         = hxx;
    hxy         = zeros(1,Ntaps);
    hyx         = hxy;

    mu          = rx.CMA.mu;

    if rx.CMA.plot == 1
        disp('Hello')
        Eout = zeros(2,Naffichage);
    
%         for jter = 1:Naffichage
%              Xin = Input(1,indice:(indice+Ntaps-1));
%              Yin = Input(2,indice:(indice+Ntaps-1));  
%     
%              Eout(1,jter) = hxx*Xin.' + hxy*Yin.';
%              Eout(2,jter) = hyx*Xin.' + hyy*Yin.';
%     
%              indice = indice+Ni;
%         end

        figure(51)
            subplot(3,2,1)
                polarplot(angle(Eout(1,:)), abs(Eout(1,:)),'.')
                title('FSE Initialisation PolX')
            subplot(3,2,2)
                polarplot(angle(Eout(2,:)), abs(Eout(2,:)),'.')
                title('FSE Initialisation PolY')
    end
    
    for iter = 1:Napprentissage
           
         isample= 1+(iter-1)*Ni;
         Xin    = Input(1,isample:(isample+Ntaps-1));
         Yin    = Input(2,isample:(isample+Ntaps-1));
              
         Xout   = hxx*Xin.' + hxy*Yin.';
         Yout   = hyx*Xin.' + hyy*Yin.';
     
         epsx   = 2*(abs(Xout)^2 - 1)*Xout;
         epsy   = 2*(abs(Yout)^2 - 1)*Yout;
          
         hxx    = hxx - mu*epsx*conj(Xin);
         hxy    = hxy - mu*epsx*conj(Yin);
         hyx    = hyx - mu*epsy*conj(Xin);
         hyy    = hyy - mu*epsy*conj(Yin);
        
         try
             if rx.CMA.plot == 1
                 % affichage de l'evolution de la constellation en temps reel
                 % sur 1000 Points, toutes les 500 iterations
                 if (mod(iter,1000) == 1)
                    ibegin = isample;
        
                    for jter = 1:Naffichage
                        Xin     = Input(1,ibegin:(ibegin+Ntaps-1));
                        Yin     = Input(2,ibegin:(ibegin+Ntaps-1));  
        
                        Eout(1,jter) = hxx*Xin.' + hxy*Yin.';
                        Eout(2,jter) = hyx*Xin.' + hyy*Yin.';
        
                        ibegin  = ibegin+Ni;
                    end 
                     
                    figure(51)
                        subplot(3,2,3)
                            polarplot(angle(Eout(1,:)), abs(Eout(1,:)),'.')
                            title('CMA + FSE')
                        subplot(3,2,4)
                            polarplot(angle(Eout(2,:)), abs(Eout(2,:)),'.')
                            title(sprintf('After %3d symbols', iter'))
                    
                        subplot(3,2,5)
                            stem(abs(hxx),'b')
                            hold on
                            stem(-abs(hxy),'k')
                            xlabel('Taps number')
                            %legend('hxx','hxy')
                            title('hxx & hxy')
                            axis([1,Ntaps,-1.2,1.2])
                            
                        subplot(3,2,6)
                            hold on
                            stem(abs(hyy),'b')
                            stem(-abs(hyx),'k')
                            xlabel('Taps number')
                            %legend('hyy','hyx')
                            title('hyy & hyx')
                            axis([1,Ntaps,-1.2,1.2])
                             
                    pause(0.1)
    %                 close(51)
                 end
            end
         catch

         end
    end 
    
    
    % Etape 2 : Séquence complète avec égalisation calculée précédemment

    Nsymbols    = floor(Nsamples/Ni) - Ntaps;
    Sout        = zeros(2, Nsamples/Ni);
    
    for iterOUT = 1:Nsymbols
         isample    = 1+(iterOUT-1)*Ni;
         Xin        = Input(1,isample:(isample+Ntaps-1));
         Yin        = Input(2,isample:(isample+Ntaps-1));   

         Sout(1,iterOUT) = hxx*Xin.' + hxy*Yin.';
         Sout(2,iterOUT) = hyx*Xin.' + hyy*Yin.'; 
    end  
    
    % Normalisation des amplitudes des sorties X & Y
    Sout(1,:) = Sout(1,:) / mean(abs(Sout(1,:)));
    Sout(2,:) = Sout(2,:) / mean(abs(Sout(2,:)));
    
    tmpx = Sout(1,:);
    tmpy = Sout(2,:);

    tmpx = tmpx(tmpx~=0);
    tmpy = tmpy(tmpy~=0);
%     mmm  = min(length(tmpx),length(tmpy));
%     size_min = [mmm,2];
%     scatterplot(tmpx)
%     scatterplot(tmpy)

    clear Sout

%     Sout = zeros(size_min);
%     Sout(:,1) = tmpx(1:mmm);
%     Sout(:,2) = tmpy(1:mmm);

    Sout(:,1) = tmpx;
    Sout(:,2) = tmpy;

    Epmdc.lambda    = Ecdc.lambda;
    Epmdc.field     = Sout;


end

