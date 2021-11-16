% All rights reserved by Jixin Chen
% First coding: 2016/04/05 by Jixin Chen @ Department of Chemistry and Biochemistry, Ohio University
% 20170110 Jixin Chen modified it to a function
% 20180609 Jixin Chen simplified it into single curve fitting

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

 %   load('jcfit2DExampleData.mat');
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
      
% x=[0.25 0.5 1 1.5 2 3 4 6 8];
% y=[19.21 18.15 15.36 14.10 12.89 9.32 7.45 5.24 3.01];

    %% -------% the rest can be a double expnential function which only need feed
    % with x and y data.---------------------------
    
    %function [parafinal, rsq] = yournamedoubleExp(x1,y1)
    
    figure; plot(x,y); title('raw data');


    %------------END of loading data: x and y in row vectors--------
    
    % set fitting options
    option.maxiteration = 50;  % number of iteration fixed, the fitting will stop either this iteration or convergence reached first 
    option.accuracy = 1E-3;  % best searching accuracy
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
    %--------------END of fitting option setting, equation, initial guess, and 

    %------------------and start fitting:------------------------
    
 
    tic
     [paraHist, parafinal, paraBounds_95, chisq, rsq] = jcfit(mdl, x, y, paraGuess, bounds, option);
    toc
%    fprintf(['\n rsq = ', num2str(rsq), '\n']);
    % parafinal is the fitted results; yfit is the fittered curve; 
    % use residual = y-yfit; to get the residual
    % rsq: root mean sqare value best will be close to 1
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
    
    % End. by Jixin Chen @ Ohio University