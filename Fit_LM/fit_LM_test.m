% coded by Jixin Chen @ Ohio University 11/2021
% http://www2.imm.dtu.dk/pubdb/edoc/imm3215.pdf
% https://people.duke.edu/~hpgavin/ce281/lm.pdf

clear

%% generate raw data
paratrue = [1.22, 2206.3, 6.23, 6958.321, 2.15]; 
%           A1, tau1, A2, tau2, baseline

A1 = paratrue(1);
tau1 = paratrue(2);
A2 = paratrue(3);
tau2 = paratrue(4);
bs = paratrue(5);

exp2 = @(x) A1*exp(-x/tau1) + A2*exp(-x/tau2) + bs;
x = 0:100000;
y = exp2(x);

noise1 = randn(1,length(x))/10;
noise2 = randn(1,length(x), 'like', y)/10;

y1 = y + noise1;
y2 = y + noise2;

figure; plot(x, y); hold on; plot(x, y1); plot(x, y2);

% figure; plot(x, y); hold on; plot(x, y1); plot(x, y2); 

 %% load raw data
 
x = x; % a row vector
y = y2; % a row vector
figure; plot(x,y); title('raw data');
    
    
% x=[0.25 0.5 1 1.5 2 3 4 6 8];
% y=[19.21 18.15 15.36 14.10 12.89 9.32 7.45 5.24 3.01];
%------------END of loading data: x and y in row vectors--------


%% initialize the fitting

x = x(:); % column
y = y2(:);% column


% step 0: set fitting options
option.miu_init = 0.0001; % initial searching step
option.v_init = 2; % searching step scaler
option.eps2 = 1e-3; % % targeting precision when on average this decimal does not change.
option.n_maxiters = 50; % maximum iteration number

% step1: set parameters
paraGuess = [5, 2000, 3, 8000, 3];
mdl = @(para, x) para(1)*exp(-(x/para(2))) + para(3)*exp(-(x/para(4))) + para(5);
Jacobian = @(para, x) [exp(-x/para(2)),  para(1)/para(2)^2*x.*exp(-x/para(2)), exp(-x/para(4)), para(3)/para(4)^2*x.*exp(-x/para(4)),  ones(length(x), 1)]; 
% Jacobian(i, m) = df(i)/dp(m), where i is ith data point, 
% and m is mth parameter. It is the slope of the residual^2 at the
% parameter space.
% J(i,:) = [exp(-x(i)/para(2)),  para(1)/para(2)^2*x(i)*exp(-x(i)/para(2)), exp(-x(i)/para(4)), para(3)/para(4)^2*x(i)*exp(-x(i)/para(4)),  1];
% mdl = @(a_buess, b_guess, x) a_guess*exp(-b_guess*x(i));
% J(i,:) = [exp(-b_guess*x(i)),  -a_guess*x(i)*exp(-b_guess*x(i))]; % [df/da, df/db]

    
%%  Levenberg–Marquardt Method using numerical method to find the Jacobian matrix
[paraHist, parafinal, sigma_p, chisq, rsq] = fit_LM(mdl, Jacobian, x, y, paraGuess, option);


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



