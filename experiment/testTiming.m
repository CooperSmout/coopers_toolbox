clear m id a b difs

% exp = 2;

for t = 1:8
%     trialLag(t) = timing(opt.mode).trials(t) - timing(2).trials(t-1);
    
    for f = 2:size(timing(opt.mode).frames,2)
        frameLag(t,f-1) = timing(opt.mode).frames(t,f)-timing(opt.mode).frames(t,f-1);
    end
    [a b] = sort(frameLag(t,:));
    m(t,:) = a;
    id(t,:) = b;
end

id(:,end-5:end)
m(:,end-5:end)

% max(difs(:))





% for F = 1:20
%     
% f1 = cgflip;
% wait(2)
% f2 = cgflip;
% 
% t(F) = f2-f1;
% end