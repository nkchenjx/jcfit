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


%% load and treat data for the proposed model.
    load('jcfitExampleData.mat');

    x1 = exampledata_x'; % a cln vector
    y1 = exampledata_y'; % a cln vector
      
    
    %-------% the rest can be a double expnential function which only need feed
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

    %------------END of loading data: x and y in clumn vectors--------
    
    % set fitting options
    option.maxiteration = 10;  % number of iteration fixed, the fitting will stop either this iteration or convergence reached first 
    option.accuracy = 0.0001;  % best searching accuracy
    option.convgtest = 1e-10; % difference between two iterations on the square difference between fitting and data.

    % ----------------Attn: change the fitting equations in the file jcftf.m -----------------
    % -------------------in function mdl ----------------
    
    % initial guess
    paraGuess = [0.5, 10, 0.5, 100, 0.1];  % A1, tau1,  A2, tau2, baseline
    % boundarys
    bounds = [0, 0.1, 0, 0.1, -0.1;   % lower boundary
              1, 2000, 1, 20000, 0.2]; % upper boundary
    %-------------------------------
    if option.accuracy > min(abs(paraGuess))/1E4
        option.accuracy = min(abs(paraGuess))/1E4;
    end % to accomodate for the acuracy requirment of the smallest parameter. Set the search mininum 1E-4* the smallest parameter.
    

    d1 = paraGuess-bounds(1,:);
    d2 = bounds(2,:)-paraGuess;
    if prod(d1)*prod(d2)<=0
        display('WARNING: initial guess out of boundary');
    end
    %--------------END of fitting option setting, equation, initial guess, and 

    %------------------and start fitting:------------------------
    [parafinal, xfit, yfit, residual, chisq, rsq] = jcfitf(x, y, paraGuess, bounds, option);
    
    display('--------------plot fitting results (x,y), (xfit, yfit), and residual -----------');
    display(parafinal);
    fprintf(['\n r square = ', num2str(rsq), '\n']);
    fprintf(['\n chi square = ', num2str(chisq), '\n']);

    % parafinal is the fitted results; yfit is the fittered curve; 
    % use residual = y-yfit; to get the residual
    % rsq: root mean sqare value best will be close to 1
%


% by Jixin Chen
