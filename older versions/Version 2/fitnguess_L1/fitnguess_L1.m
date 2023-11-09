% All rights reserved by Jixin Chen
% First coding: 2016/04/05 by Jixin Chen @ Department of Chemistry and Biochemistry, Ohio University
% 20170110 Jixin Chen modified it to a function
% 20180609 Jixin Chen simplified it into single curve fitting
% 20200617 Jixin Chen add fitting error analysis
% 20211116 Jixin Chen cleaned the code and change a few variable names for
%                     easy to read
% 

% Copyright (c) 2018 Jixin Chen @ Ohio University
% 
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

%% an example for triple exponential decay fitting of a single curve.
%{

clear

%% generate raw data
paraTrue = [100.22, 105.3, 6.23, 958.321, 15.33, 8230 2.15]; 
%           A1, tau1, A2, tau2, baseline

A1 = paraTrue(1);
tau1 = paraTrue(2);
A2 = paraTrue(3);
tau2 = paraTrue(4);
A3 = paraTrue(5);
tau3 = paraTrue(6);
bs = paraTrue(7);

exp2 = @(x) A1*exp(-x/tau1) + A2*exp(-x/tau2) + A3*exp(-x/tau3) + bs;
x = 0:100000;
y = exp2(x);

noise1 = randn(1,length(x))/10;
noise2 = randn(1,length(x), 'like', y)/10;

y1 = y + noise1;
y2 = y + noise2; % white noise with equal weight

% figure; plot(x, y); hold on; plot(x, y1); plot(x, y2); 

    
    x = x; % a row vector
    y = y2; % a row vector
    figure; plot(x,y); title('raw data');
    
    
% x=[0.25 0.5 1 1.5 2 3 4 6 8];
% y=[19.21 18.15 15.36 14.10 12.89 9.32 7.45 5.24 3.01];
%------------END of loading data: x and y in row vectors--------

%% ----Setting up fitting model and parameters-------------
    %           the rest can be a double expnential function, any custom function 
    %           or a group of functions in a separated Matlab
    %           function. Just pass the function handle to the fitting
    %           funciton, e.g.
    %           function [yfit1, yfit2, yfit3,...] = yournamedoubleExp(para, x1, x2, x3,...)
    %                 
    %           All functions are numerical solved with no efforts to fit
    %           analytically with x and y data.
    %-----------------------------------------

    % set fitting options
    option.maxiteration = 50;  % number of iteration fixed, the fitting will stop either this iteration or convergence reached first 
    option.precision = 1E-3;  % best searching precision, recommend 1 decimal better than desired. e.g want 0.01, set to 0.001.
    option.convgtest = 1e-100; % difference between two iterations on the square difference between fitting and data.

    % ----------------Attn: change below for different fitting equations-----------------
    % set the fitting equation to double exponential decay with a base line
    mdl = @(para, x) para(1)*exp(-(x/para(2))) + para(3)*exp(-(x/para(4))) + para(5)*exp(-(x/para(6))) + para(7);

%    mdl = @(para, x) para(1)*exp(-(x/para(2))); 
    % equation grammar: modle name 'mdl' use as function y = mdl(para, x), the rest is the equation.
    % you can also generate a function to replace the above equation: 
    % function newy = mdl(para, x)
    % which will allow combination of equations and global fitting using
    % different equations for different pieces of data that share some
    % common parameters.
    
    % initial guess
    paraGuess = [7, 10, 4, 300, 8, 1000, 3];  % A1, tau1,  A2, tau2, baseline
%     paraGuess = [1, 33];
    % boundarys
    bounds = [0, 0, 0, 0, 0, 0, 0;   % lower boundary
              1E3, 1E5, 1E3, 1E5, 1E3, 1E5, 10]; % upper boundary
%     bounds = [0, 0,;   % lower boundary
%               100, 100]; % upper boundary

    %-------------------------------

    d1 = paraGuess-bounds(1,:);
    d2 = bounds(2,:)-paraGuess;
    if prod(d1)*prod(d2)<=0
        display('WARNING: initial guess out of boundary');
    end
    %--------------END of fitting option setting, equation, initial guess,
    %              and guessed parameter boundaries.

    
    %------------------and start fitting:------------------------
     
    tic
         [paraHist, parafinal, paraBounds_95, chisq, rsq] = fitnguess_L1(mdl, x, y, paraGuess, bounds, option);
    toc
%     fprintf(['\n rsq = ', num2str(rsq), '\n']);
    % parafinal is the fitted results; yfit is the fittered curve; 
    % use residual = y-yfit; to get the residual
    % rsq: root mean sqare value best will be close to 1
    
    %--------- plot results -----------------
    yfit = mdl(parafinal, x);
    residual = y - yfit;
    figure; plot(x,y,'linewidth',1.5); hold on; plot(x,yfit,'linewidth',1.5); plot(x, residual,'linewidth',1.5);
    title(['rsq = ', num2str(rsq)]);
    ax = gca;
    ax.LineWidth = 1.5;
    ax.Box = 'on';
    ax.TickLength = [0.02, 0.02];
    ax.FontName = 'Arial';
    ax.FontSize = 20;
    ax.FontWeight = 'Bold';
    
    %-------------------------------------
    % End. by Jixin Chen @ Ohio University

 
%}


%% main function with data error bar not fitted
% -------if mdl is a complicated function, remove it from here and add it
% in the end as a separated function.
function [paraHist, parafinal, paraBounds_95, chisq, rsq] = fitnguess_L1(mdl, x, y, paraGuess, bounds, option)
    % load options
    if isfield(option, 'maxiteration')
        maxiteration = option.maxiteration;
    else
        maxiteration = 50;  % number of iteration fixed
    end
    
    if isfield(option, 'precision')
        precision =option.precision;
    else
        precision = 0.0001;    % best searching precision default 0.0
    end
    
    if isfield(option, 'convgtest')
        convgtest = option.convgtest;
    else
        convgtest = 1e-10; % difference between two iterations on the square difference between fitting and data.
    end
    if isfield(option,'step')
        step = option.step;
    else
        step = 0.5; % spacing fraction. Spacing is 2^(step*i)*precision, i is the number of guessing point away from the initial guess.
    end

   
    
    y_guess = mdl(paraGuess, x);
    residual = y-y_guess;
%    minerrorlast = dot(residual, residual); %L2, least square
    minerrorlast = sum(abs(residual(:))); %L1, least absolute various
  
    
    
    para = paraGuess;
    paraHist = zeros(maxiteration, length(paraGuess) + 1); % last one error score
    
     tic
     for iteration = 1:maxiteration % fixed number of iteration.
         paraHist(iteration,:) = [para, minerrorlast]; 
         
         fprintf('.');
         if rem(iteration, 100) == 0 % progressing indicator
             fprintf('\n');
         end
         for i = 1:length(paraGuess) % scan each parameter
             %set the scanning scale withing the boundary.
            p = para(i);
%             if abs(p) > precision
%                 precision = abs(p)*precision;
%             end        
            lb = bounds(1, i);
            ub = bounds(2, i);
            ll = p-lb;
            nl = floor(log2(ll/precision+1)):-step:1;
            ul = ub-p;
            nu = 1:step:floor(log2(ul/precision+1));
            ps = [p, lb, p-2.^nl*precision, p+2.^nu*precision, ub];
            error = zeros(length(ps),1);
             % scan the parameter across the scale
             for j = 1: length(ps)
                para(i) = ps(j);
                residual = y - mdl(para,x);
                %error(j) = sum((mdl(para, x)-y).^2); %---the key equation: sum square of residual
                % error(j) = dot(residual, residual); % square residual, dot() is ~5% faster than sum()
                error(j) = sum(abs(residual(:))); % L1 norm
             end
             % find the best
             [minerror, ind] = min(error); % find the least square.
             para(i) = ps(ind(1)); % if there are multiple minima, take the first, could be p.
         end
         %test if converged
         if abs(minerror-minerrorlast)< convgtest; %convergence test positive
             fprintf(['\n converged at iteration = ', num2str(iteration)]);
             break
         end
         minerrorlast = minerror;
     end
    fprintf('\n');
    toc


     parafinal = para;

    [paraBounds_95, chisq, rsq] = finderror(mdl, x, y, parafinal, bounds, precision);
    
%-------------Note out the figure plotting or cut it to main function if
%-------------you don't like it here:
%     % plot figures 
%     yfit = mdl(parafinal, x);
%     figure; plot(x,y,'linewidth',1.5); hold on; plot(x,yfit,'linewidth',1.5); plot(x, residual,'linewidth',1.5);
%     title(['rsq = ', num2str(rsq)]);
%     ax = gca;
%     ax.LineWidth = 1.5;
%     ax.Box = 'on';
%     ax.TickLength = [0.02, 0.02];
%     ax.FontName = 'Arial';
%     ax.FontSize = 20;
%     ax.FontWeight = 'Bold';
%--------end plotting figure    
end

%% find errors of the fitting
function [paraBounds_95, chisq, rsq] = finderror(mdl, x, y, para, bounds, precision)
% parafinal structure: parameter, lower boundary, upper bounday at 95
% confidence 2 sigma.
% error structure: upper error and lower error. Take the average as the
% final for 1 sigma.
%precision = precision*10;

num_para = length(para);

 yfit =  mdl(para, x);
 
 residual = y - yfit;
 a = (y - yfit).^2./yfit;
 a(isinf(a)) = 0;
 chisq = sum(a);
 meany = mean(y);
 sumy =  sum((y-meany).^2);
 N =  length(x);

 
sumr = sum(residual(:).^2);
rsq = 1- sumr/sumy;

sigma = sqrt(sumr/N);
% meanrsd = mean(residual(:));

% upper bounds
% usigma = zeros(num_para);
upbounds = zeros(num_para, 1);
lowbounds = zeros(num_para, 1);
error = zeros(num_para, 2);

for i = 1: num_para
    parau = para;
    p = parau(i);
    ub = bounds(2,i);
    ul = ub-p;
    nu = 1:0.2:floor(log2(ul/precision+1));
    psu = [p, p+2.^nu*precision, ub];
    sigmau = zeros(length(psu),1);
    for f = 1: length(psu)
       parau(i) = psu(f);
     %   esqu(i) = 0;

      yfitu = mdl(parau, x);
   
      %%% -------------------- use y + y_error as upper limit if y_error is
      %%% avaailable (this idea needs reference and justification):    
      residualu =  y - yfitu;
 %     residualu = y+y_error - yfitu; %------------------------
       %    esqu(i) = esqu(i) + sum((yfitu{j} -y).^2);
      sigmau(f) = sqrt(sum(residualu(:).^2)/N);
    end
    a = abs(sigmau - sigma);
    b1 = abs(a - 1*sigma/sqrt(N-num_para)); % 66% confidence level 1 sigma
    b2 = abs(a - 2*sigma/sqrt(N-num_para)); % 95% confidence level 2 sigma
    [~, ind1] = min(b1);
    [~, ind2] = min(b2);
    upbounds(i) = psu(ind2(1));         % 2 sigma 95 confidence upper bounds
    error(i,1) = psu(ind1(1))-para(i); % 1 sigma error upper bound 66% confidence
  
 end

%usigma = sigma*precision./abs(sigmau-sigma);%/sqrt(N-num_para);


for i = 1: num_para
    paral = para;
    p = paral(i);
    lb = bounds(1,i);
    ll = p-lb;
    nl = floor(log2(ll/precision+1)):-0.2:1;
    psl = [lb, p-2.^nl*precision, p];
    sigmal = zeros(length(psl),1);
    for f = 1: length(psl)
       paral(i) = psl(f);
     %   esqu(i) = 0;
    
      yfitl =  mdl(paral, x);
       %%% -------------------- use y - y_error as upper limit if y_error is
      %%% avaailable to calculate resitual (this idea needs reference and justification): 
      residuall =  y - yfitl;
  %     residualu = y - y_error - yfitu; %----------------------------
      sigmal(f) = sqrt(sum(residuall(:).^2)/N);
    end
    a = abs(sigmal - sigma);
    b1 = abs(a - 1*sigma/sqrt(N-num_para)); % 66% confidence level 1 sigma
    b2 = abs(a - 2*sigma/sqrt(N-num_para)); % 95% confidence level 2 sigma
    [~, ind1] = min(b1);
    [~, ind2] = min(b2);
    lowbounds(i) = psl(ind2(1));          % 2 sigma 95 confidence lower bounds
    error(i,2) = para(i) - psl(ind1(1));  % 1 sigma  error lower bounds 66% confidence
end

%----
% lower and upper bounds of the parameters at 95% confidence or 2 sigma of the final noise.
paraBounds_95 = [lowbounds(:)'; upbounds(:)']; 
    % warning: the parameter 95% confidence lower and upper bounds are based on estimation of the local minimum,
    % not considering global minimum and other local minima.
%------
% parafinal2 = [para; lowbounds; upbounds];
end

% End. by Jixin Chen @ Ohio University
