function X = iFT(Xft)

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
%   We aim to get the distance between an INPUT and a REFERENCE object in a
%   N dimensions space with a given (or not) calculation METHOD.
%   Default method is L² distance.
%
% ----- INPUTS -----
%   varargin{1}     INPUT
%   varargin{2}     REFERENCE
%   varargin{3}     METHOD
%
% ----- BIBLIOGRAPHY -----
% ---------------------------------------------


X = fftshift(ifft(fftshift(Xft)));