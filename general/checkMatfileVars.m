function [answer] = checkMatfileVars(m,varargin)

% CHECKMATFILEVARS Reads matfile m, and returns true if m contains field
%
%
%   created by Cooper Smout (ORCID: 0000-0003-1144-3272)

    answer = false(1,length(varargin));
    vars = who(m);
    for v = 1:length(varargin)
        if any(strcmp(vars,varargin{v}))
            answer(v) = 1;
        end
    end
    
    