function En = fn(varargin)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : FN - FIELD NORMALISATION
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-09-09
%   Version         : 1.1
%
% ----- Main idea -----
%   Normalise the fields by wether their own \sqrt{mean power} or by a
%   given \sqrt{mean power} (OPTIONNAL).
%
% ----- INPUTS -----
%   EIN:    (structure) containing the Fields to be normlised
%               - LAMBDA    [nm]        wavelength
%               - FIELD     [sqrt(mW)]  Normalised electric fields
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------
    
    % checking number of input arguments
    assert(nargin <3,"Wrong number of input arguments." + ...
        " Should be (ENN) or (ENN,PNORM).")

    Enn = varargin{1};

    % input type flexibility
    if isstruct(Enn) == 0
        tmp.lambda = 1550;
        tmp.field  = Enn;
        clear Enn
        Enn = tmp;
    end
    
    if nargin == 2
        Pnorm = varargin{2};
    end
     
    En.lambda = Enn.lambda;
    
    % n_polar = 2;
    if size(Enn.field,2) == 2*size(Enn.lambda,1)
    
        % split the polarisations
        [Ennx,Enny]= sep_XYfields(Enn);
        Enx = Ennx;
        Eny = Enny;
            
        % normalise each polarisation
        Xnn = Ennx.field;
        Ynn = Enny.field;

        if nargin == 1
            Xn  = Xnn/sqrt(mean(empower(Xnn)));
            Yn  = Ynn/sqrt(mean(empower(Ynn)));
        else
            Xn = Xnn/sqrt(Pnorm);
            Yn = Ynn/sqrt(Pnorm);
        end

        Enx.field = Xn;
        Eny.field = Yn;

        % recombine both polarisations
        En    = merge_XYfields(Enx,Eny);
        
    % n_polar = 1;
    else
        if nargin == 1
            En.field = Enn.field/sqrt(mean(empower(Enn)));
        else
            En.field = Enn.field/sqrt(Pnorm);
        end
    end


end

