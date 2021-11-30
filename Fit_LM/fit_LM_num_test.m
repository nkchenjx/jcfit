% coded by Jixin Chen @ Ohio University 11/2021
% http://www2.imm.dtu.dk/pubdb/edoc/imm3215.pdf
% https://people.duke.edu/~hpgavin/ce281/lm.pdf

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


