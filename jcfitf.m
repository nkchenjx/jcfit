% by Jixin Chen
% First coding: 2016/04/05 by Jixin Chen @ Department of Chemistry and Biochemistry, Ohio University
% 20170110 Jixin Chen modified it to a function
% 20180609 Jixin Chen simplified it into single curve fitting
% 20190426 Jixin Chen package it onto GitHub

%-----------------------------------
% For publications please cite: 
% https://pubs.acs.org/doi/abs/10.1021/acs.jpcb.6b05697
% Jixin Chen, Joseph R Pyle, Kurt Waldo Sy Piecco, Anatoly B Kolomeisky,
% Christy F Landes, A Two-Step Method for smFRET Data Analysis, 
% J. Phys. Chem. B, 2016, 120 (29), pp 7128–7132


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

%% an example for double exponential decay fitting of a single curve.
%{
  %  load('jcfitExampleData.mat');

    x1 = exampledata_x; % a row vector
    y1 = exampledata_y; % a row vector
      
    
    %% -------% the rest can be a double expnential function which only need feed
    % with x and y data.---------------------------
    
    %function [parafinal, rsq] = yournamedoubleExp(x1,y1)
    
    figure; plot(x1,y1); title('raw data');
    % this is a double exponential decay curve that have a dark part
    % find time 0 and remove the data before it starts, and normalize the curve 
    [m,ind] = max(y1);
    l = length(y1);
    x = x1-x1(ind);
    x = x(ind:l);
    y = y1(ind:l);
    y = y/m; 
    figure; plot(x,y); title('treated data'); 

    %------------END of loading data: x and y in row vectors--------
    
    % set fitting options
    option.maxiteration = 10;  % number of iteration fixed, the fitting will stop either this iteration or convergence reached first 
    option.accuracy = 0.0001;  % best searching accuracy
    option.convgtest = 1e-10; % difference between two iterations on the square difference between fitting and data.

    % ----------------Attn: change below for different fitting equations-----------------
    % set the fitting equation to double exponential decay with a base line
    mdl = @(para, x) para(1)*exp(-(x/para(2))) + para(3)*exp(-(x/para(4))) + para(5);
    % equation grammar: modle name 'mdl' use as function y = mdl(para, x), the rest is the equation.
    % you can also generate a function to replace the above equation: 
    % function newy = mdl(para, x)
    
    % initial guess
    paraGuess = [0.5, 10, 0.5, 100, 0.1];  % A1, tau1,  A2, tau2, baseline
    % boundarys
    bounds = [0, 0.1, 0, 0.1, -0.1;   % lower boundary
              1, 2000, 1, 20000, 0.2]; % upper boundary
    %-------------------------------

    d1 = paraGuess-bounds(1,:);
    d2 = bounds(2,:)-paraGuess;
    if prod(d1)*prod(d2)<=0
        display('WARNING: initial guess out of boundary');
    end
    %--------------END of fitting option setting, equation, initial guess, and 

    %------------------and start fitting:------------------------
    [parafinal, yfit, chisq, rsq] = jcfit(mdl, x, y, paraGuess, bounds, option);
    fprintf(['\n rsq = ', num2str(rsq), '\n']);
    % parafinal is the fitted results; yfit is the fittered curve; 
    % use residual = y-yfit; to get the residual
    % rsq: root mean sqare value best will be close to 1
%}



%% main function with data error bar not fitted
function [parafinal2, xfit, yfit, residual, chisq, rsq] = jcfitf(x, y, paraGuess, bounds, option)
    % load options
    if ~isfield(option,'maxiteration')
     option.maxiteration = 50;  % number of iteration fixed
     disp('set default interation to 50');
    end
    if ~isfield(option,'accuracy')
     option.accuracy = 0.0001;    % best searching accuracy
     disp('set default decimal to 1E-4');
    end
    if ~isfield(option,'convgtest')
     option.convgtest = 1e-10; % difference between two iterations on the square difference between fitting and data.
     disp('set default convergent minimum to 1E-10');
    end
    if ~isfield(option,'numdatafit')
     option.numdatafit = 100; % How many points to show in the fitted data plot
  %   disp('set default fitted data point to 100');
    end
    %------------SearchExpScale is very important, the smaller the more
    %exhaustive searching with more time consuming------------------
    if ~isfield(option,'searchExpScale')
     option.searchExpScale = 0.5; % default guessing spacing is n*exp(0.5)*accuracy times from the initial guess, where n is number of guesses away from the initial guess.
     disp('set default search space to be exp(0.5)');
    end
    maxiteration = option.maxiteration;
    accuracy = option.accuracy;
    convgtest = option.convgtest;
    numdatafit = option.numdatafit;
    searchExpScale = option.searchExpScale;
    minerrorlast = 0;
       
    para = paraGuess;
    
    
     tic
     for iteration = 1:maxiteration % fixed number of iteration.
         fprintf('.');
         if rem(iteration, 100) == 0 % progressing indicator
             fprintf('\n');
         end
         for i = 1:length(paraGuess) % scan each parameter
             %set the scanning scale withing the boundary.
            p = para(i);
%             accuracy = option.accuracy(i); % if accuracy is a vector for each parameter
%             if abs(p) > option.accuracy
%                 accuracy = abs(p)*option.accuracy;
%             else
%                 accuracy = option.accuracy;
%             end
            lb = bounds(1, i);
            ub = bounds(2, i);
            ll = p-lb;
            nl = floor(log2(ll/accuracy+1)):-searchExpScale:1;
            ul = ub-p;
            nu = 1:searchExpScale:floor(log2(ul/accuracy+1));
            ps = [lb, p-2.^nl*accuracy, p+2.^nu*accuracy, ub];
            error = [];
             % scan the parameter across the scale
             for j = 1: length(ps)
                para(i) = ps(j);
                error(j) = sum((mdl(para, x)-y).^2); %---the key equation: sum square of residual
             end
             % find the best
             [minerror, ind] = min(error); % find the least square.
             para(i) = ps(ind);
         end
         %test if converged
         if abs(minerror-minerrorlast)< convgtest; %convergence test positive
             fprintf('\n converged at iteration = ', num2str(iteration));
             break
         end
         minerrorlast = minerror;
     end
    fprintf('\n');
    toc
    %options=optimset('Maxiter',100000,'TolFun',1e10);
    % [para,r,J,COVB,mse] = nlinmultifit(x_cell, y_cell, mdl_cell, para0);

    % Calculate results

    %[ypred1, delta1] = nlpredci(mdl1, x, parareal, r, 'covar', mse);
    % figure

    parafinal = para;
    yfit = mdl(para, x);
    
    residual = y - yfit;
    a = (y - yfit).^2./yfit;
    a(isinf(a)) = 0;
    chisq = sum(a);

    sumy = sum(y.^2);
    N = length(x);

    sumr = sum(residual.^2);
    rsq = 1- sumr/sumy;
    
    if length(x) < numdatafit
         xfit = [x(1):(x(end)-x(1))/numdatafit : x(end)];
         xfit = xfit';
         yfit = mdl(para, xfit);
    else
        xfit = x;
    end
    
    % plot figures
    figure; plot(x,y,'linewidth',1.5); hold on; plot(xfit,yfit,'linewidth',1.5); plot(x, residual,'linewidth',1.5);
    title(['rsq = ', num2str(rsq)]);
    ax = gca;
    ax.LineWidth = 1.5;
    ax.Box = 'on';
    ax.TickLength = [0.02, 0.02];
    ax.FontName = 'Arial';
    ax.FontSize = 20;
    ax.FontWeight = 'Bold';

    [parafinal2, error, chisq, rsq] = finderror(x, y, para, bounds, accuracy);
    
end


function [parafinal2,error, chisq, rsq] = finderror(x, y, para, bounds, accuracy)
% parafinal structure: parameter, lower boundary, upper bounday at 95
% confidence 2 sigma.
% error structure: upper error and lower error. Take the average as the
% final for 1 sigma.
%accuracy = accuracy*10;

num_para = length(para);

chisq = 0;
rsq = 0;


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
for i = 1: num_para
    parau = para;
    p = parau(i);
    ub = bounds(2,i);
    ul = ub-p;
    nu = 1:0.2:floor(log2(ul/accuracy+1));
    ps = [p, p+2.^nu*accuracy, ub];
    sigmau = [];
    for f = 1: length(ps)
       parau(i) = ps(f);
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
    [c1, ind1] = min(b1);
    [c2, ind2] = min(b2);
    upbounds(i) = ps(ind2(1));
    error(i,1) = ps(ind1(1))-para(i);
  
 end

%usigma = sigma*accuracy./abs(sigmau-sigma);%/sqrt(N-num_para);

% lower bounds
a = [];
b1 = [];b2 = [];
c1 = []; c2 = [];
ind1 = []; ind2 = [];
ps = [];
for i = 1: num_para
    paral = para;
    p = paral(i);
    lb = bounds(1,i);
    ll = p-lb;
    nl = floor(log2(ll/accuracy+1)):-0.2:1;
    ps = [lb, p-2.^nl*accuracy, p];
    sigmal = [];
    for f = 1: length(ps)
       paral(i) = ps(f);
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
    [c1, ind1] = min(b1);
    [c2, ind2] = min(b2);
    lowbounds(i) = ps(ind2(1));
    error(i,2) = para(i) - ps(ind1(1));
 end

parafinal2 = [para; lowbounds; upbounds];
end



%% custum your own model here for any functions including global fitting functions or multiple peak functions
%------------------ Fitting model -----------------------------------------
function [Yguess] = mdl(para, x) 
% Yguess = para(1)*exp(-x/para(2));
mdl2 = @(para, x) para(1)*exp(-(x/para(2))) + para(3)*exp(-(x/para(4))) + para(5);

Yguess = mdl2(para, x);
end

% by Jixin Chen