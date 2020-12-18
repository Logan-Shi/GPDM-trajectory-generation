function [q, rej, qtraj] = myhmc(q0, stepsize, traj, n, f, df, varargin)
%
% [q, rej, qtraj] = hmc(q0, stepsize, traj, n, f, df)
%
% Metropolis algorithm with Gaussian proposals to sample from 
% a general distribution.
%
% q0 -- Initial state as vector.
% stepsize -- Step size of states as vector.
% traj -- Trajectory length of Markov chain.
% n -- Number of samples to obtain.
% f -- Target distribution.
% df -- Gradient of target distribution.
% varargin -- Input arguments to f and df.
% q -- Matrix of states, one on each row.  q0 is first state.
% rej -- Number of rejections.
%

  seed = 93136785; 
  
  dim = size(MakeCol(q0), 1);
  q(1, :) = q0;
  rej = 0;

  % For each sample
  for i = 2 : n + 1

    q1 = q(i - 1, :);

    % Save leapfrog trajectory position if required
    if (nargout == 3)   
	qtraj{i-1}(1,:) = q1;  
    end

    % Use the same seed for random number generator 
    %  randn ('seed', mod (seed * i, 159)  );

    % Initialize momentum variable from a unit, isotropic Gaussian
    p1 = randn(size(q1));
 
    % Initial kinetic energy 
    KE = 0.5 * sum (p1 .^ 2);
    
    % Initial potential energy 
    %    (This is the log posterior we want samples from,
    %     which we take to be Gaussian for this example)
    PE = feval(f, q1, varargin{:});
    
    % Intial Hamiltonian
    h1 = KE + PE;
    
    % Calculate gradient of energy E with respect to q
    %   (this is just the - gradient of the log posterior evaluated at q1)
    dEdq = feval(df, q1, varargin{:});
    
    % Initial half-step update of p
    p1 = p1 - 0.5 * ( stepsize .* dEdq );

    % Leapfrog trajectory
    for len = 1 : traj 
      
      % Full-step update of q
      q1 = q1 + stepsize .* p1;

      % Save leapfrog trajectory position if required
      if (nargout == 3)  
          qtraj{i-1}(len+1,:) = q1;  
      end
	
      % Recalculate gradient of q
      dEdq = feval(df, q1, varargin{:});

      if (len ~= traj) 
	    % Full-step update of p
	    p1 = p1 - stepsize .* dEdq;
      else
	    % Last half-step update of p
	    p1 = p1 - 0.5 * ( stepsize .* dEdq );
      end

    end
   
    % Final kinetic energy
    KE = 0.5 * sum (p1 .^ 2);
    
    % Final potential energy
    PE = feval(f, q1, varargin{:});
    
    % Final Hamiltonian
    h2 = KE + PE;

  %   rand ('seed', mod (seed, 13))
    % Metropolis acceptance test
    if ( h2 < h1 ) | ( rand < exp( h1 - h2 ) )
      q(i, :) = q1;
    else
      q(i, :) = q(i - 1, :);
      rej = rej + 1;
    end

  end  % All samples obtained

 
