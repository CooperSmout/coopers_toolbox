%%
clear all
close all
clc

% add convolve2 function to path
addpath('./convolve2');

% set test parameters
ppd = 40; % pixels per degree of visual angle
n = 1000; % number of test images

% set noisy Gabor parameters
cfg          = [];
cfg.patchsiz = ppd*8.0; % patch size (pix)
cfg.patchenv = ppd*1.0; % patch spatial envelope s.d. (pix)
cfg.patchlum = 0.5; % patch background luminance
cfg.gaborper = ppd*1.5; % Gabor spatial period (pix)
cfg.gaborang = pi/180*45; % Gabor orientation (rad)
cfg.gaborphi = 0; % Gabor unit phase
cfg.gaborcon = 0; % Gabor Michelson contrast
cfg.noisedim = cfg.gaborper/6; % noise dimension (pix)
cfg.noisecon = 0.5/3; % noise RMS contrast

% precompute contrast weights
cfg_x = get_patch_contrast(cfg);

% precompute signal energy filters
gaborang = pi/180*(0:5:180); % filter orientations (rad)
for j = 1:length(gaborang)
    cfg_e{j} = get_patch_energy(setfield(cfg,'gaborang',gaborang(j)));
end

% make test images
x = nan(n,1); % RMS contrast
e = nan(n,length(gaborang)); % signal energies
hbar = waitbar(0,'making test images...');
set(get(findobj(hbar,'Type','Axes'),'Title'),'FontSize',16);
for i = 1:n
    if mod(i,10) == 0, waitbar(i/n,hbar); end
    % clip patch luminance to [0,1] range
    while true
        img = make_gabor_noise(cfg);
        if all(img(:) > 0 & img(:) < 1)
            break
        end
    end
    % get weighted RMS constrast
    x(i,1) = get_patch_contrast(cfg_x,img);
    % get signal energies
    for j = 1:length(gaborang)
        e(i,j) = get_patch_energy(cfg_e{j},img);
    end
end
close(hbar);

% save results
save('test_gabor_noise.mat','n','gaborang','x','e','cfg');

%%
clear all
close all
clc

% load results
load('test_gabor_noise.mat');

% plot correlation profiles
i = find(gaborang == pi/4); % reference energy channel
re = corr(e,e(:,10)); % correlation between signal energies
rx = corr(e,x); % correlation between signal energies and RMS contrast
figure;
hold on
plot(0:5:180,re);
plot(0:5:180,rx);
hold off
xlabel('orientation (deg)');
ylabel('correlation strength');

% simulate hypothetical observer
eprime = [8,4]; % raw energy sensitivities for the two target orientations
% compute decision variable
j = find(ismember(gaborang,pi/4*[1,3]));
dv = sum(bsxfun(@times,e(:,j),eprime),2);
% simulate observer many times
for i = 1:1000
    resp = dv+sqrt(2)*randn(n,1) > mean(dv);
    % estimate energy sensitivities
    b = glmfit([e(:,j),x],resp,'binomial','link','probit');
    eprime_hat(i,:) = b(2:3)*sqrt(2);
    xprime_hat(i,1) = b(4);
end
