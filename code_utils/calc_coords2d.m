function [xx,yy] = calc_coords2d(xva,yva,xres,yres)
[xx,yy] = meshgrid(linspace(-xva,xva,xres),linspace(yva,-yva,yres));