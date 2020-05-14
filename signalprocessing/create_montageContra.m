% create_montageLR

% create contralateralised montage for use in Fieldtrip
% scalp channels: average-referenced
% EOG channels: bipolar referenced


% DEFINE LEFT MONTAGE
montage_biosemi64_L.labelorg = {    
    'Fp1'
    'AF7'
    'AF3'
    'F1'
    'F3'
    'F5'
    'F7'
    'FT7'
    'FC5'
    'FC3'
    'FC1'
    'C1'
    'C3'
    'C5'
    'T7'
    'TP7'
    'CP5'
    'CP3'
    'CP1'
    'P1'
    'P3'
    'P5'
    'P7'
    'P9'
    'PO7'
    'PO3'
    'O1'
    'Iz'
    'Oz'
    'POz'
    'Pz'
    'CPz'
    'Fpz'
    'Fp2'
    'AF8'
    'AF4'
    'AFz'
    'Fz'
    'F2'
    'F4'
    'F6'
    'F8'
    'FT8'
    'FC6'
    'FC4'
    'FC2'
    'FCz'
    'Cz'
    'C2'
    'C4'
    'C6'
    'T8'
    'TP8'
    'CP6'
    'CP4'
    'CP2'
    'P2'
    'P4'
    'P6'
    'P8'
    'P10'
    'PO8'
    'PO4'
    'O2'
    'LH'
    'RH'
    'LV'
    'UV'}; % Nx1 cell-array
montage_biosemi64_L.labelnew = {
    'Fp1/2i'
    'AF7/8i'
    'AF3/4i'
    'F1/2i'
    'F3/4i'
    'F5/6i'
    'F7/8i'
    'FT7/8i'
    'FC5/6i'
    'FC3/4i'
    'FC1/2i'
    'C1/2i'
    'C3/4i'
    'C5/6i'
    'T7/8i'
    'TP7/8i'
    'CP5/6i'
    'CP3/4i'
    'CP1/2i'
    'P1/2i'
    'P3/4i'
    'P5/6i'
    'P7/8i'
    'P9/10i'
    'PO7/8i'
    'PO3/4i'
    'O1/2i'
    'Iz'
    'Oz'
    'POz'
    'Pz'
    'CPz'
    'Fpz'
    'Fp1/2c'
    'AF7/8c'
    'AF3/4c'
    'AFz'
    'Fz'
    'F1/2c'
    'F3/4c'
    'F5/6c'
    'F7/8c'
    'FT7/8c'
    'FC5/6c'
    'FC3/4c'
    'FC1/2c'
    'FCz'
    'Cz'
    'C1/2c'
    'C3/4c'
    'C5/6c'
    'T7/8c'
    'TP7/8c'
    'CP5/6c'
    'CP3/4c'
    'CP1/2c'
    'P1/2c'
    'P3/4c'
    'P5/6c'
    'P7/8c'
    'P9/10c'
    'PO7/8c'
    'PO3/4c'
    'O1/2c'
    'H-EOG' 
    'V-EOG'}; % Mx1 cell-array

num_org = length(montage_biosemi64_L.labelorg);
num_new = length(montage_biosemi64_L.labelnew);

montage_biosemi64_L.tra = zeros(num_new,num_org);


% average reference scalp channels
for N = 1:num_new-2
    montage_biosemi64_L.tra(N,N)             = 1;
    montage_biosemi64_L.tra(N,1:N-1)         = -1/(num_new-3);
    montage_biosemi64_L.tra(N,N+1:end-4)       = -1/(num_new-3);
end

% horizontal EOG
montage_biosemi64_L.tra(num_new-1,num_org-3) = +1;
montage_biosemi64_L.tra(num_new-1,num_org-2) = -1;

% vertical EOG
montage_biosemi64_L.tra(num_new,num_org-1)   = +1;
montage_biosemi64_L.tra(num_new,num_org)     = -1;



% DEFINE RIGHT MONTAGE
montage_biosemi64_R = montage_biosemi64_L;
montage_biosemi64_R.tra(1:3,:) = montage_biosemi64_L.tra(34:36,:);
montage_biosemi64_R.tra(4:11,:) = montage_biosemi64_L.tra(39:46,:);
montage_biosemi64_R.tra(12:27,:) = montage_biosemi64_L.tra(49:64,:);
montage_biosemi64_R.tra(34:36,:) = montage_biosemi64_L.tra(1:3,:);
montage_biosemi64_R.tra(39:46,:) = montage_biosemi64_L.tra(4:11,:);
montage_biosemi64_R.tra(49:64,:) = montage_biosemi64_L.tra(12:27,:);




% check
montage_biosemi64_R.labelnew(1:3,:),  montage_biosemi64_R.labelorg(34:36,:)
montage_biosemi64_R.labelnew(4:11,:),  montage_biosemi64_R.labelorg(39:46,:)
montage_biosemi64_R.labelnew(12:27,:),  montage_biosemi64_R.labelorg(49:64,:)
montage_biosemi64_R.labelnew(34:36,:),  montage_biosemi64_R.labelorg(1:3,:)
montage_biosemi64_R.labelnew(39:46,:),  montage_biosemi64_R.labelorg(4:11,:)
montage_biosemi64_R.labelnew(49:64,:),  montage_biosemi64_R.labelorg(12:27,:)



save('montage_biosemi64contra','montage_biosemi64_L','montage_biosemi64_R')



