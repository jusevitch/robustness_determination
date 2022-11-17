function outvec = de2bifree(varargin)

% Takes an integer and returns its binary representation
% This is NOT a copy of the paid "de2bi" function. Rather, it does the same
% job but isn't hidden behind a paywall
% Inputs:
%   1. Integer to convert to binary
%   2. (optional) Desired vector length

% Reference: https://www.mathworks.com/matlabcentral/answers/89526-binary-string-to-vector

if nargin == 0
    error('de2bifree -- Needs at least one argument');
elseif nargin == 1
    outvec = dec2bin(varargin{1})-'0';
    Iinv = rot90(eye(length(outvec)));
    outvec = outvec*Iinv; % Reversed order of matrix and vector since outvec is a row vector
elseif nargin == 2
    outvec = dec2bin(varargin{1})-'0';
    Iinv = rot90(eye(length(outvec)));
    outvec = outvec*Iinv;
    if length(outvec) < varargin{2}
        outvec = [outvec zeros(1,varargin{2} - length(outvec))];
    end
else
    error('de2bifree -- Too many input arguments')
end

end