function Xft = FFT(X)

% ---------------------------------------------
% ----- INFORMATIONS -----
%   Function name   : GET_DISTANCE
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk.work@gmail.com
%   Date            : 2022-09-17
%   Version         : 1.0
%
% ----- MAIN IDEA -----
% ----- INPUTS -----
%   varargin{1}     INPUT
%   varargin{2}     REFERENCE
%   varargin{3}     METHOD
%
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------


Xft = fftshift(fft(fftshift(X)));