function varargout = smoothLine(varargin)

smoothAmt = 1000;

if nargin >= 2
    x = varargin{1};
    y = varargin{2};
end
if nargin >= 3
    smoothAmt = varargin{3};
end

% The type of interpolation can be:
% 'cubic' - creates smooth data points for a smooth line (best results).
% 'pchip' - same as cubic.
% 'nearest' - creates a square wave effect.
% Do NOT use interpolation:
% 'linear' - results are as though this function were not used.
% 'spline' - resulting data/line will not be an accurate representation of the data. Might as well use the SPLINE function.

xInterp = x(1):(x(2)-x(1))/smoothAmt:x(end); % Set up the mesh.
yInterp = interp1(x,y,xInterp,'cubic');

varargout{1} = xInterp;
varargout{2} = yInterp;