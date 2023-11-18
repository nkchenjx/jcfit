% Copyright (c) 2021 Jixin Chen @ Ohio University

% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

% coded by Jixin Chen @ Ohio University 11/2021
% http://blog.sina.com.cn/s/blog_4a8e595e01014tvb.html
% the Levenberg–Marquardt Method using numerical method to find the
% Jacobian matrix
% http://www2.imm.dtu.dk/pubdb/edoc/imm3215.pdf

%% example testing code
%{
clear

%% generate raw data
paraTrue = [50.22, 20.3, 10.23, 258.321, 20.5, 819, 10.33, 2430, 5.2, 6456, 2.15]; 
%           A1, tau1, A2, tau2, ..., baseline

A1 = paraTrue(1);
tau1 = paraTrue(2);
A2 = paraTrue(3);
tau2 = paraTrue(4);
A3 = paraTrue(5);
tau3 = paraTrue(6);
A4 = paraTrue(7);
tau4 = paraTrue(8);
A5 = paraTrue(9);
tau5 = paraTrue(10);
bs = paraTrue(11);

exp5 = @(x) A1*exp(-x/tau1) + A2*exp(-x/tau2) + A3*exp(-x/tau3) +  A4*exp(-x/tau4) +  A5*exp(-x/tau5) + bs;
x = 0:10:10000;
y = exp5(x);

noise1 = randn(1,length(x))/10;
noise2 = randn(1,length(x), 'like', y)/10;

y1 = y + noise1;
y2 = y + noise2; % white noise with equal weight

% figure; plot(x, y); hold on; plot(x, y1); plot(x, y2); 

 %% load raw data
 
x = x; % a row vector
y = y2; % a row vector
figure; plot(x,y); title('raw data');
    
    
% x=[0.25 0.5 1 1.5 2 3 4 6 8];
% y=[19.21 18.15 15.36 14.10 12.89 9.32 7.45 5.24 3.01];
%------------END of loading data: x and y in row vectors--------


%% initialize the fitting


option.miu = 0.0001; % 0.1 to 0.0001, minimum searching steps. 
option.v = 2; % 2 to 10, times miu to give a new searching step on each failed trial
option.eps2 = 1e-20; % targeting precision when on average this decimal does not change.
option.n_maxiters = 500; % maximum iteration number


% guess model and initial parameter
mdl = @(para, x) para(1)*exp(-(x/para(2))) + para(3)*exp(-(x/para(4))) + para(5)*exp(-(x/para(6)))...
                            + para(7)*exp(-(x/para(8)))+ para(9)*exp(-(x/para(10))) + para(11);

paraGuess = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];  % A1, tau1,  A2, tau2,..., baseline

    
%%  Levenberg–Marquardt Method using numerical method to find the Jacobian matrix
[paraHist, parafinal, sigma_p, chisq, rsq] = fit_LM_num(mdl, x, y, paraGuess, option);


%% report results
y_fit = mdl(parafinal, x);
residual = y - y_fit;
figure; plot(x,y); hold on; plot(x, y_fit);
    title(['rsq = ', num2str(rsq)]);
    ax = gca;
    ax.LineWidth = 1.5;
    ax.Box = 'on';
    ax.TickLength = [0.02, 0.02];
    ax.FontName = 'Arial';
    ax.FontSize = 20;
    ax.FontWeight = 'Bold';

hold on; plot(x, residual); title(['fitting and residual, R^2= ', num2str(rsq)]);

figure; plot(paraHist(:, end)); title('history of square residual score over iterations');
disp(['Fitted final parameters = ', num2str(parafinal)]);
disp(['sigma (66%) = ', num2str(sigma_p)]);
%}

%% evenberg–Marquardt fitting function http://www2.imm.dtu.dk/pubdb/edoc/imm3215.pdf
function [paraHist, parafinal, sigma_p, chisq, rsq] = fit_LM_num(mdl, x, y, paraGuess, option)
    
    
    %step 0: set fitting options
    miu_init = option.miu; % =0.0001; % 0.1 to 0.0001
 %   miu_init = 10;
    v_init = option.v; % =2; % 2 to 10
    eps2 = option.eps2; %=1e-3; % targeting precision when on average this decimal does not change.
    n_maxiters = option.n_maxiters; %=50; % maximum iteration number
    
    % step1: set parameters
%     paraGuess = [5, 2000, 3, 8000, 3];
%     mdl = @(para, x) para(1)*exp(-(x/para(2))) + para(3)*exp(-(x/para(4))) + para(5);

    % Jacobian = @(para, x, i) [exp(-x(i)/para(2)),  para(1)/para(2)^2*x(i)*exp(-x(i)/para(2)), exp(-x(i)/para(4)), para(3)/para(4)^2*x(i)*exp(-x(i)/para(4)),  1]; 
    % Jacobian(i, m) = df(i)/dp(m), where i is ith data point, 
    % and m is mth parameter. It is the slope of the residual^2 at the
    % parameter space.
    % J(i,:) = [exp(-x(i)/para(2)),  para(1)/para(2)^2*x(i)*exp(-x(i)/para(2)), exp(-x(i)/para(4)), para(3)/para(4)^2*x(i)*exp(-x(i)/para(4)),  1];
    % mdl = @(a_buess, b_guess, x) a_guess*exp(-b_guess*x(i));
    % J(i,:) = [exp(-b_guess*x(i)),  -a_guess*x(i)*exp(-b_guess*x(i))]; % [df/da, df/db]

    N_data=length(y);
    n_para=length(paraGuess);

    %find initial error score
    y_guess = mdl(paraGuess, x);
    residual = y-y_guess;
    error = dot(residual, residual);

    miu = miu_init;
    v = v_init;

    % step2: run Levenberg–Marquardt iterations
    updateJ = 1;
    para = paraGuess';

    paraHist = zeros(n_maxiters, n_para + 1); % last one error score
    
    tic
    for it=1:n_maxiters
        paraHist(it,:) = [para', error/length(y)];
        fprintf('.');
        if rem(it, 100) == 0 % progressing indicator
            fprintf('\n');
        end

        if updateJ == 1
            % calculate the Jacobian matrix
            J = Jacobian_num(mdl, para, x, miu);
            % calculate residual
            y_guess = mdl(para, x);
            residual = y-y_guess;
            % Calculate the 
            A=J'*J;
        end

        % calculate the new point to try
        H_lm = A +(miu*eye(n_para, n_para)); % Hessian matrix
        g = J'*residual(:);
        dp = H_lm\g;
     %   dp(isnan(dp)) = 0;
        para_lm = para + dp;


        % calculate the new residual and rou value, gain ratio
        y_guess_lm = mdl(para_lm, x);
        residual_lm = y - y_guess_lm;
        error_lm = dot(residual_lm, residual_lm);
        rou = (error - error_lm)/(dp'*(miu*dp + g));

        % check if the new point is good
        if rou > 0
            para = para_lm;
            error = error_lm;
%            miu = miu*max(1/3, 1-(2*rou-1)^3);
            miu = miu*max(miu_init, 1-(2*rou-1)^3);
%            disp(miu);
            v = v_init;
            updateJ = 1; % update to the new point and recalculate the Jacobian matrix 

            %test convergence
            if dot(dp, dp) <= (eps2*(eps2 + error_lm)) % reach the target precision
                disp('convergent at the designed precision');
                break;
            end
        else
            miu = miu*v; 
            v = v_init*v;
            updateJ = 0; % recalculate dp using new miu and v
        end

    end
    toc

    %report results
    parafinal = para';
    y_fit = mdl(parafinal, x);
    residual = y - y_fit;
%     figure; plot(x,y); hold on; plot(x, y_fit);
%     figure; plot(x, residual); title('residual');
    
    
    %---------------
    % error estimation of sigma, chi^2, R^2
    J = Jacobian_num(mdl, parafinal, x, miu_init);

    weight = ones(N_data, 1); % asssuming equal weight of all data points
    
    JtWJ  = J' * ( J .* ( weight * ones(1,n_para) ) );  
  
    covar_p = inv(JtWJ);
    sigma_p = sqrt(diag(covar_p));
    sigma_p = sigma_p';  % 1sgima, 66%
    
%     H = J'*J +(miu*eye(n_para, n_para)); % Hessian matrix
%     g = J'*residual(:);
%     dp = H\g;
%        
%     sigma = abs(std(y-y_fit)^2*(N_data - n_para)*dp +1/2*H*dp.^2);
% 
%     sigma = abs(4*error*dp);
    chisq = 0;
    rsq = 0;
 
    a = (residual).^2./y_fit;
    a(isinf(a)) = 0;
    a(isnan(a)) = 0;
    chisq = sum(a);
 
    meany = mean(y);
    sumy =  sum((y-meany).^2);
    sumr = sum(residual(:).^2);
    rsq = 1- sumr/sumy;
    
end


%% numerically calculate the Jacobian matrix
function J = Jacobian_num(mdl, para, x, miu) % numerically calculate the Jacobian matrix
  lp  = length(para);
  lx = length(x);
  J = zeros(lx, lp);
  
  for p = 1:lp
      parap = para;
      parap(p) = para(p)+miu;
      p1 =  mdl(para,x)';
      p2 =  mdl(parap,x)';
      Jp = (p2-p1)/miu;
      J(:, p) = Jp;
  end 
end



