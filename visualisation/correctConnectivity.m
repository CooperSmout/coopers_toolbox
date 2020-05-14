%% Function to fix correlation matrix after swapping electrode order
% accepts up to 12 dimensions AFTER the first two electrode dimensions
% pass in second argument as 0 to replace obsolete values with zeros (rather than NaN)

function data = correctConnectivity(data,~)

% set defaults
if nargin<2
    obs = NaN; % value to overwrite top-right half of matrix
else
    obs = 0;
end
    

% check data
numElecs=size(data,1);
if numElecs~=size(data,2)
    error('matrix not square')
end


% cycle through and correct all conditions (if they exist)... there must be a better way of doing this!!!
for a = 1:size(data,3)
    for b = 1:size(data,4)
        for c = 1:size(data,5)
            for d = 1:size(data,6)
                for e = 1:size(data,7)
                    for f = 1:size(data,8)
                        for g = 1:size(data,9)
                            for h = 1:size(data,10)
                                for k = 1:size(data,11)
                                    for l = 1:size(data,12)
                                        for m = 1:size(data,13)
                                            for n = 1:size(data,14)
                                                

            % move flipped values back into bottom half of matrix
            for E1 = 1:numElecs
                for E2 = E1+1 : numElecs % only check top(bottom) right half
                    if abs(data(E1,E2,a,b,c,d,e,f,g,h,k,l,m,n))>0
                        if abs(data(E2,E1,a,b,c,d,e,f,g,h,k,l,m,n))>0
                            flag=1;
                        end
                        data(E2,E1,a,b,c,d,e,f,g,h,k,l,m,n)=data(E1,E2,a,b,c,d,e,f,g,h,k,l,m,n); % overwrite data in bottom half of matrix
                    end
                    data(E1,E2,a,b,c,d,e,f,g,h,k,l,m,n) = obs; % reset to NaN or zero
                end
            end

            % warning
            if exist('flag')
                warning('values present in top half of matrix have overwritten non-zero values in bottom half...')
            end
                                             
                                            end
                                        end
                                    end
                                end
                            end         
                        end
                    end
                end                  
            end
        end
    end                      
end
                                    
