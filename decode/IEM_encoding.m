


function [chan_resp_cv_coeffs, chan_resp_cv_coeffs_shift] = IEM_encoding(curr_trials,curr_reps,curr_orientationAngles,curr_trnX,orientationsExperiment,chan_center,windowSize)

error ('DEPRECATED, USE IEM_encoding2')

% 
if nargin<7
    windowSize = 1;
end


% re-order dimensions: etT -> Tet
numFeatures = size(curr_trials,1);
if numFeatures>64
    warning('More than 64 features detected')
elseif numFeatures<64
    warning('Less than 64 features detected')
end
curr_trials = permute(curr_trials,[3 1 2]);


% Cross-validated training/testing on contrast discrimination task
chan_resp_cv_coeffs = nan(size(curr_trials,1),length(chan_center),size(curr_trials,3));
tic
for ii = 1:max(curr_reps)
    fprintf('.')

    trainID = curr_reps~=ii;
    testID = curr_reps==ii;

    train_trials = curr_trials(trainID,:,:);
    test_trials = curr_trials(testID,:,:);

    % loop over timepoints
%     for tt = 1:size(train_trials,3)
    for tt = 1:size(train_trials,3)-windowSize+1
        
        % average over time window
        win = tt:tt-1+windowSize;
        
        % training
        w_coeffs = curr_trnX(trainID,:)\mean(train_trials(:,:,win),3); 
%         w_coeffs = curr_trnX(trainID,:)\train_trials(:,:,tt); 

        % testing
        chan_resp_cv_coeffs(testID,:,tt) = (inv(w_coeffs*w_coeffs.')*w_coeffs*mean(test_trials(:,:,win),3).').'; % can also be written (w_coeffs.'\thistst_tpt.').';

    end
end
toc



% coregistered reconstructions - adjust trials so that the correct orientation falls on 100 deg
targ_ori = chan_center(round(length(chan_center)/2));
targ_ori_idx = find(chan_center==targ_ori);

chan_resp_cv_coeffs_shift = nan(size(chan_resp_cv_coeffs));
for ii = 1:length(orientationsExperiment)
    thisidx = curr_orientationAngles==orientationsExperiment(ii);
    chan_resp_cv_coeffs_shift(thisidx,:,:) = circshift(chan_resp_cv_coeffs(thisidx,:,:), targ_ori_idx-find(orientationsExperiment(ii)==chan_center) , 2 );
end



% re-order dimensions: Tet -> etT
chan_resp_cv_coeffs = permute(chan_resp_cv_coeffs,[2 3 1]);
chan_resp_cv_coeffs_shift = permute(chan_resp_cv_coeffs_shift,[2 3 1]);




