function decoder = train_beamformer2(cfg0, X, Y)
% [decoder] = train_beamformer(cfg, X, Y)
%    Trains a linear decoder "beamformer style" to optimally recover the latent components as 
%    prescribed in X. Several decoders may be trained indepedently, corresponding to several
%    latent components.
%
%    X           Vector or matrix of size C x N, where C is the number of components and N is
%                the number of trials, that contains the expected/prescribed component activity
%                in the training data.
%    Y           Matrix of size F x N, where F is the number of features, that contains the
%                training data.
%    cfg         Configuration struct that can possess the following fields:
%                .gamma = [scalar]                Shrinkage regularization parameter, with range [0 1]. 
%                                                 Default is to compute automatically (see covdiag.m).
%                .discardNan = 'yes' or 'no'      Whether trials with NaN in either X or Y should be
%                                                 removed prior to training. Default is 'no'.
%                .returnPattern = 'yes' or 'no'   Whether the spatial patterns of the components should
%                                                 be returned. Default = 'no';
%                .demean = 'yes' or 'no'          Whether to demean the data first (per feature, over
%                                                 trials). Default = 'yes';.
%                .demeanX = 'yes' or 'no'         Whether to demean the design matrix. Default = 'no';.
%
%    decoder     The (set of) decoder(s), that may be passed on to an appropriate decoding function,
%                e.g. decode_beamformer. It may contain a field .pattern of size C x F
%
%    See also DECODE_BEAMFORMER.

%    Created by Pim Mostert, 2016
%       edits:  added option to automatically estimate noise covariance matrix  - Sept 2017 Cooper Smout
%               demeaned design matrix (reported in Kok2017 but not implemented here) - 1/3/18 Cooper Smout

if ~isfield(cfg0, 'discardNan')
    cfg0.discardNan = 0;
end
if ~isfield(cfg0, 'returnPattern')
    cfg0.returnPattern = 'no';
end    
if ~isfield(cfg0, 'demean')
    cfg0.demean = 'yes';
end    
if ~isfield(cfg0, 'demeanX')
    cfg0.demeanX = 'no';
end    

decoder = [];

% Tranpose
Y = Y';
X = X';

if cfg0.discardNan
    X = X(~any(isnan(Y), 2), :);
    Y = Y(~any(isnan(Y), 2), :);
end

numF = size(Y, 2);
numC = size(X, 2);
numN = size(Y, 1);

% Demean data 
if strcmp(cfg0.demean, 'yes')
    mY = mean(Y, 1);
    Y = Y - repmat(mY, [numN, 1]);
    decoder.mY = mY';
end

% Demean design matrix --- edited Cooper Smout April 2018
if strcmp(cfg0.demeanX, 'yes')
    mX = mean(X, 1);
    X = X - repmat(mX, [numN, 1]);
end

% Estimate filters
if strcmp(cfg0.returnPattern, 'yes')
    decoder.pattern = zeros(numF, numC);
end
decoder.W = zeros(numC, numF);

for ic = 1:numC
    
    % Estimate leadfield for current channel    
    l = (X(:, ic)'*X(:, ic))\(X(:, ic)'*Y);
    if strcmp(cfg0.returnPattern, 'yes')
        decoder.pattern(:, ic) = l';
    end
    
    % Estimate noise
    N = Y - X(:, ic)*l; 
    if ~isfield(cfg0, 'gamma') % edited Cooper Smout Sept 2017
        S = covdiag(N); 
    else
        S = cov(N); % Estimate noise covariance
        S = (1-cfg0.gamma)*S + cfg0.gamma*eye(numF)*trace(S)/numF; % Regularize
    end
    
    % Calculate filter
    decoder.W(ic, :) = l/S;
    decoder.W(ic, :) = decoder.W(ic, :)/(decoder.W(ic, :)*l');
end

end
