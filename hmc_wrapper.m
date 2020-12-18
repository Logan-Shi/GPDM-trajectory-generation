function [samples, energies, diagn] = hmc_wrapper(f, x, options, gradf, varargin)

ioptions = foptions; 
ioptions(14) = 10000;
num_samples = options(14); 
num_skip = options(15);
num_test = 25;

x = scg(f, x, ioptions, gradf, varargin{:});
ix = x; 

while (1)
    options(14) = 1;
    options(15) = 0;
    
    x = ix; 
    [samples, energies, diagn] = hmc(f, x, options, gradf, varargin{:});
    if diagn.acc(1) < 1e-4
        ratio = 0.1;
    else
        options(14) = num_test;
        options(15) = num_skip;
        
        x = ix;
        [samples, energies, diagn] = hmc(f, x, options, gradf, varargin{:});

        diagn.acc(find(diagn.acc > 1)) = 1;
        ratio = sum(diagn.acc)/size(diagn.acc, 1)
        options(18)
    end
    if (ratio < 0.6)
        options(18) = 0.7*options(18);
    end
    if (ratio > 0.9)
        options(18) = 1.5*options(18);
    end
    if (ratio > 0.6 && ratio < 0.95)
        options(15) = 0;
        options(14) = num_samples-num_test;
        
        if (options(14) > 0) 
            [samples2] = hmc(f, samples(end,:), options, gradf, varargin{:});
            samples = [samples; samples2];
        end
        break;
    end
end