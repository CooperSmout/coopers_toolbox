
function data = filterdata_cs(data,type,lowfreq,highfreq,lowpole,highpole)
% this function filters the data signal with the appropriate filter and
% returns the filtered data

switch type
    case 'Butterworth'
        if ~isempty(highpole) && ~isempty(highfreq) && ~isnan(highpole) && ~isnan(highfreq)
%             [B,A] = butter(highpole,highfreq*2/data); % makes the low pass filter
            [B,A] = butter(highpole,1/highfreq); % makes the low pass filter
            data = filtfilt(B, A, data); % filters the data
        end
        if ~isempty(lowpole) && ~isempty(lowfreq) && ~isnan(lowpole) && ~isnan(lowfreq)
%             [B,A] = butter(lowpole,lowfreq*2/data,'high'); % high pass
            [B,A] = butter(lowpole,1/lowfreq,'high'); % high pass
            data = filtfilt(B, A, data);
        end
%     case 'Elliptic'
%         Rp = 0.5; % Rp decibels of peak-to-peak ripple
%         Rs = 20; %minimum stopband attenuation of Rs decibels
%         if ~isempty(highpole) && ~isempty(highfreq) && ~isnan(highpole) && ~isnan(highfreq)
%             [B,A] = ellip(highpole,Rp,Rs,highfreq*2/data); % makes the low pass filter
%             data = filtfilt(B, A, data); % filters the data
%         end
%         if ~isempty(lowpole) && ~isempty(lowfreq) && ~isnan(lowpole) && ~isnan(lowfreq)
%             [B,A] = ellip(lowpole,Rp,Rs,lowfreq*2/data,'high'); % high pass
%             data = filtfilt(B, A, data);
%         end
end
% ---------------------------- end of function ----------------------------
