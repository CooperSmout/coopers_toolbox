%% Script for executing Expectation x Relevance Experiment 1
% Cooper Smout, Queensland Brain Institute, University of Queensland


clear all



%% PSYCHTOOLBOX

AssertOpenGL;

% screen
Screen('Preference', 'SkipSyncTests', 0)
screenNumber=0;%max(Screen('Screens'));
[ptb.main,rect]=PsychImaging('OpenWindow',screenNumber,0);
if ~isequal(rect, [0 0 1680 1050])
   error('check resolution') 
end
HideCursor();



%%
%%
%% COLOUR

m = hsv(12);
m(13,:)=[.5 .5 .5];
m2 = rgb2hsv(m);
col = 9;

notFin = 1;

while notFin
    
    Screen('FillRect',ptb.main,255*m(col,:),rect)
    Screen('Flip',ptb.main);

    [a b c]=KbStrokeWait(-3);

    switch find(b)
        case 38
            m2(col,3)=m2(col,3)+.1;

        case 40
            m2(col,3)=m2(col,3)-.1;

        case 37
            m2(col,3)=m2(col,3)-.01;
            
        case 39
            m2(col,3)=m2(col,3)+.01;
            
        case 13
            col = col+1;
            
        case 32
            notFin=0;
            
    end
    
    
    
    m = hsv2rgb(m2);

end
    

sca



