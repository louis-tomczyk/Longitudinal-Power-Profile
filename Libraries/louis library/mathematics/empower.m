function power = empower(field)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : EMPOWER
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-07-21
%   Version         : 1.0
%
% ----- Main idea -----
%   Put to the power an input signal
% 
% ----- INPUTS -----
%   FIELD   (structure/array)   The signal to put to |.|^2
%
% ----- OUTPUTS -----
%   POWER   (array)             The empowered signal
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------

    if isstruct(field) == 1
        power = abs(field.field).^2;
    else
        power = abs(field).^2;
    end
end