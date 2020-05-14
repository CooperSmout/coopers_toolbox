function EXP = readTriggerWaituntil(EXP,t)

% read trigger when function first called
trigOld = readTrigger(EXP.trig.ioObj,EXP.trig.ioIn,EXP.trig.ioInOff); 

while time <  t

    % trigger state
    trigCurr = readTrigger(EXP.trig.ioObj,EXP.trig.ioIn,EXP.trig.ioInOff); % read trigger

    % detect trigger state change
    if trigCurr && ~trigOld % heartbeat detected (trigger off -> on)

        EXP.hbCurr(end+1) = time; % record heartbeat time
        fprintf([num2str(time-EXP.trialStartTemp) ' '])
        
    end
    
    % reset trigger for next pass
    trigOld = trigCurr;
    
%     EXP.temp(end+1) = time-EXP.trialStartTemp;
    
end   