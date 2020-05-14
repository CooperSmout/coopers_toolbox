function v = pickTargets(v,opt)


% determine target times for attended side
frameGap = opt.targGap + 2*opt.targRamp; % frames gap on either side of target
v.framesPerImage = opt.mon_rate./(2*v.freqs);
v.targetImages = {[] []};
v.targetFrameStart = {[] []};
v.targetFramePeak = {[] []};
v.targetFrameEnd = {[] []};

allFrames = 1:v.frames;
allFrames(1 : opt.signalRamp * opt.mon_rate) = NaN;
allFrames(v.frames-frameGap : end) = NaN;


% cycle through targets, starting with cued side
for side = [v.cue 3-v.cue]

    % determine target times for each side
    for targ = v.n_targets(side):-1:1

        availableFrames = allFrames(~isnan(allFrames));
        
        % pick frame for target to appear
        if side==v.cue && targ==2 % ensure second target of attended side appears toward end of trial
            pick_frame = round(rand*frameGap + v.frames - 2*frameGap);
%             pick_frame = randsample(availableFrames,1);
        else
            pick_frame = randsample(availableFrames,1); % select from remaining available frames
        end
        
        % locate nearest image to picked frame
        target_image = round(pick_frame*2*v.freqs(side)/opt.mon_rate) + 1;
        targetFrame = (target_image-1) * v.framesPerImage(side) + 1;
        
        % record target images and frame times
        v.targetImages{side}(targ) = target_image;
        v.targetFrameStart{side}(targ) = targetFrame;
        v.targetFramePeak{side}(targ) = v.targetFrameStart{side}(targ) + opt.targRamp;
        v.targetFrameEnd{side}(targ) = v.targetFrameStart{side}(targ) + opt.targRamp*2 - 1;
        
        % update frames register
        allFrames(  max(1,targetFrame-frameGap)  :  min(v.frames,targetFrame+frameGap)  ) = NaN;
        
%         v
%         allFrames'
%         availableFrames'
%         pause
        
    end
    
end

