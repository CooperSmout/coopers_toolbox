% Bayes Analysis script by Zoltan Dienes
% http://www.lifesci.sussex.ac.uk/home/Zoltan_Dienes/inference/Bayesfactor.html

% modified by Cooper Smout (ORCID: 0000-0003-1144-3272) to turn into function 30/7/2015
% can take input values for uniform distibution ONLY, or run without inputs to enter values manually


function stat  = bayes(sampleMeanDiff,SEM,lower,upper)


normaly = @(mn, variance, x) 2.718283^(- (x - mn)*(x - mn)/(2*variance))/realsqrt(2*pi*variance); 


% enter values if function passed without inputs
if nargin
    uniform = 1;
else
    sampleMeanDiff = input('What is the sample mean (difference?)? ');
    SEM = input('What is the sample standard error (press enter if only have t score)? ');

    % NOTE: SEM = s/sqrt(n)
    % where 
    %   s = sample standard deviation of sample
    %   n = size of sample 

    if isempty(SEM)
        tscore = input('What is the t score between conditions? ');
        SEM = sampleMeanDiff/tscore
    end

    uniform = input('is the distribution of p(population value|theory) uniform? 1= yes 0=no ');

    if uniform == 0
        meanoftheory = input('What is the mean of p(population value|theory)? ');
        sdtheory = input('What is the standard deviation of p(population value|theory)? ');
        omega = sdtheory*sdtheory;
        tail = input('is the distribution one-tailed or two-tailed? (1/2) ');
    end


    if uniform == 1
        lower = input('What is the lower bound? ');
        upper = input('What is the upper bound? ');
    end

end

% calculate bayes
area = 0;
if uniform == 1
    theta = lower;
else theta = meanoftheory - 5*(omega)^0.5;
end
if uniform == 1
    incr = (upper- lower)/2000;
else incr =  (omega)^0.5/200;
end

for A = -1000:1000
    theta = theta + incr;
    if uniform == 1
        dist_theta = 0;
        if and(theta >= lower, theta <= upper)
            dist_theta = 1/(upper-lower);
        end
    else %distribution is normal
        if tail == 2
            dist_theta = normaly(meanoftheory, omega, theta);
        else
            dist_theta = 0;
            if theta > 0
                dist_theta = 2*normaly(meanoftheory, omega, theta);
            end
        end
    end
    
    height = dist_theta * normaly(theta, SEM^2, sampleMeanDiff); %p(population value=theta|theory)*p(data|theta)
    area = area + height*incr; %integrating the above over theta
end

% save output
stat.data = sampleMeanDiff;
stat.SEM = SEM;
stat.lowerBound = lower;
stat.upperBound = upper;
stat.distribution = 'uniform';
stat.likelihoodTheory = area;
stat.likelihoodNull = normaly(0, SEM^2, sampleMeanDiff);
stat.bayesFactor = stat.likelihoodTheory/stat.likelihoodNull;
 
 