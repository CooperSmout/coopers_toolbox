%% create image matrices to be used with minimum motion experiment
% Cooper Smout, Queensland Brain Institute, University of Queensland

% phosphor to adjust
adjust = 1;

% image size
maxrad = 150; % outer radius of windmill
minrad = 100; % inner radius of windmill
w = 2*maxrad; % resolution x
h = 2*maxrad; % resolution y
x = 1:w;
y = 1:h;
ox = w/2;
oy = h/2;

% coordinates
[X,Y] = meshgrid(x-ox,y-oy); % cartesian coords in 400 x 400 matrix
[TH, R] = cart2pol(X,Y); % converts each cartesian coordinate to polar coords in 400 x 400 matrix
TH = rem(TH + 2*pi, 2*pi); % shift from -pi:pi to 0:2pi
inc = pi/120; % determines pixelation or not

% variables
cycles = 6; % number of times cycling in 360 degrees
lure = .16; % maximum modulation of luminance lure (x2)
rad = 0:inc:2*pi-inc; % radial segments in space (the windmill)
radnext = inc:inc:2*pi; % next radial segment beyond current one
testspeed = 1; % higher test speed -> quicker test run (minimum 1)
timesteps = 0:2*pi/(60/testspeed):2*pi - 2*pi/(60/testspeed);

gamma(1) = 2.4550;
gamma(2) = 2.3354;

% draw images for each rgb value
test_array =  5:5:255;
for ii = 1:length(test_array) ; % equivalent rgb value to be put into subsequent experiment
    
    
    realRGB = (test_array(ii)/255).^gamma(adjust)*255;
    adj = realRGB/255; % modulation of red phosphor expressed as a fraction of 1
    
    % cycle timesteps
    IM = zeros(h,w,3,length(timesteps));
    for t = 1:length(timesteps);

        time = timesteps(t);
        disp(time)

        % phosphor equations
        phosphor(:,adjust) = adj/2 - .5*adj*cos( cycles*rad )*cos( time ) + .5*adj*lure*sin( cycles*rad )*sin( time ) ; % colour grating component + luminance grating component
        phosphor(:,3-adjust) = .5 + .5*cos( cycles*rad )*cos( time ) + .5*lure*sin( cycles*rad )*sin( time ) ; % colour grating component + luminance grating component
        phosphor(:,3) = 0; 

        phosphor(:,1) = phosphor(:,1) .^ (1/gamma(1));
        phosphor(:,2) = phosphor(:,2) .^ (1/gamma(2));


        wm{adjust} = zeros(h,w);
        wm{3-adjust} = zeros(h,w);
        wm{3} = zeros(h,w);   
        
        for RR = 1:length(rad);
            wm{adjust}( ( TH >= rad(RR) & TH <= radnext(RR) ) & ( R > minrad & R < maxrad ) ) = phosphor(RR,1); % red phosphor intensity of windmill
            wm{3-adjust}( ( TH >= rad(RR) & TH <= radnext(RR) ) & ( R > minrad & R < maxrad ) ) = phosphor(RR,2); % green phosphor intensity of windmill   
%             wm{3}( ( TH >= rad(RR) & TH <= radnext(RR) ) & ( R > minrad & R < maxrad ) ) = phosphor(RR,3); % blue phosphor intensity of windmill
        end

        % concatenate
        im = cat(3, wm{adjust}, wm{3-adjust}, wm{3} );
        IM(:,:,:,t) = im;

    end
    IM = round(IM*10000)/10000;
    
    % save
    save( [ 'IM' num2str(test_array(ii)) '.mat'], 'IM')
    
end


    