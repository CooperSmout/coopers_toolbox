function data_minusERP = subtractERP(data)
% function to subtract ERP (average across trials) from each individual trial
% used in calculating non-phase-locked power (see Cohen 2014)
%
% data must be in format (elec x time x trials)
% data_minusERP is in same format

if size(data,1)~=64
    warning('64 electrodes not detected')
end

% create erp 
erp = mean(data,3);

% subtract ERP from raw trials
data_minusERP = zeros(size(data),class(data));
for t = 1:size(data,3)
    data_minusERP(:,:,t) = data(:,:,t) - erp;
end


