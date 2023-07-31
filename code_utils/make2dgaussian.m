function gauss = make2dgaussian(xysize, params)

%function gauss = make2dgaussian(xysize, params)
% 
%   returns a 2d gaussian function of size X-by-Y specified by a vector PARAMS.
%
%   Input
%       xysize >>> vector of x and y size, i.e. [64 64]
%       params(1:2) >>> center coordinates (x,y): axes are normalized
%                 to 0 (lower left corner) to 1(upper right corner)
%       params(3) >>> The direction of the 2d gaussian main axis in degree (0-360).
%       params(4:5) >>> Spatial (relative axis w. respect to orientation) envelope size in standard deviation
%
% OUTPUT:
%       gauss >>> a 2d gauss function of size X-by-Y, specified by a vector PARAMS.
%
%   2017 - Paolo Papale fecit

cx = params(1);
cy = params(2);
dir = params(3);
senvX = params(4);
senvY = params(5);

dx = 0:(1/(xysize(1)-1)):1;
dy = 0:(1/(xysize(2)-1)):1;

[iy ix] = ndgrid(dx, dy);

a = cos(dir/180*pi)^2/(2*senvX^2) + sin(dir/180*pi)^2/(2*senvY^2);
b = -sin(2*dir/180*pi)/(4*senvX^2) + sin(2*dir/180*pi)/(4*senvY^2);
c = sin(dir/180*pi)^2/(2*senvX^2) + cos(dir/180*pi)^2/(2*senvY^2);

gauss = exp(-(a*(ix-cx).^2 + 2*b*(ix-cx).*(iy-cy) + c*(iy-cy).^2));
end