function gb=gabor_fn(gaussiansd,orientation,wavelength,phaseoffset,ellipticity)
 
sigma_x = gaussiansd;
sigma_y = gaussiansd/ellipticity;
 
% Bounding box
nstds = 3;
xmax = max(abs(nstds*sigma_x*cos(orientation)),abs(nstds*sigma_y*sin(orientation)));
xmax = ceil(max(1,xmax));
ymax = max(abs(nstds*sigma_x*sin(orientation)),abs(nstds*sigma_y*cos(orientation)));
ymax = ceil(max(1,ymax));
xmin = -xmax; ymin = -ymax;
[x,y] = meshgrid(xmin:xmax,ymin:ymax);
 
% Rotation 
x_theta=x*cos(orientation)+y*sin(orientation);
y_theta=-x*sin(orientation)+y*cos(orientation);
 
% gb= exp(-.5*(x_theta.^2/sigma_x^2+y_theta.^2/sigma_y^2)).*cos(2*pi/wavelength*x_theta+phaseoffset);
gb= exp(-.5*(x_theta.^2/sigma_x^2+y_theta.^2/sigma_y^2)).*sin(2*pi/wavelength*x_theta+phaseoffset);


%% to then save gabor / non-gabor
% imwrite((gb+1)*intensity, 'name','bmp')
% imwrite(ones(size(gb))*intensity, 'name','bmp')