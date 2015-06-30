function y = fshift(x,s)
% FSHIFT Fractional circular shift
%   Syntax:
%
%       >> y = fshift(x,s)
%
%   FSHIFT circularly shifts the elements of vector x by a (possibly
%   non-integer) number of elements s. FSHIFT works by applying a linear
%   phase in the spectrum domain and is equivalent to CIRCSHIFT for integer
%   values of argument s (to machine precision).

% Author:   François Bouffard
%           fbouffard@gmail.com

needtr = 0; if size(x,1) == 1; x = x(:); needtr = 1; end;
N = size(x,1); 
r = floor(N/2)+1; f = ((1:N)-r)/(N/2); 
p = exp(-1j*s*pi*f).'; 
y = ifft(fft(x).*ifftshift(p)); if isreal(x); y = real(y); end;
if needtr; y = y.'; end;