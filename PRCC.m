% It calculates PRCCs and their significances
% (uncorrected p-value and Bonferroni correction)
% LHSmatrix: LHS matrix (N x k) 
% Y: output matrix (time x N) 
% s: time points to test (row vector) 
% PRCC_var: cell array of strings {'p1','p2',...,'pk'}
%            to label all the parameters varied in the LHS
% by Simeone Marino, May 29 2007 
% modified by Marissa Renardy, October 8 2020
% N.B.: this version uses ONE output at a time

function [prcc,sign,sign_label]=PRCC(LHSmatrix,Y,s,PRCC_var,alpha)
% Y=Y(s,:)';% Define the output. Comment out if the Y is already 
%           % a subset of all the time points and it already comprises
%           % ONLY the s rows of interest
Y = Y'; % N x time_points vector (N x 1 here)
[~,num_parameters]=size(LHSmatrix); % Define the size of LHS matrix, k = number of parameters
[~,num_time_points]=size(Y); % number of times
prcc = zeros(num_parameters, num_time_points); % k x N matrix
prcc_significance = zeros(num_parameters, num_time_points);
for i=1:num_parameters  % Loop for the whole submatrices
    % Get all 
    Z = [LHSmatrix(:, 1:i-1), LHSmatrix(:, i+1:end)];
    vals = [LHSmatrix(:,i), Y]; % N x time_points+1 matrix ( Here: N x 2 ).
    % Compute PRCCs, controlling for the other parameters in Z
    [rho,pvals]=partialcorr(vals, Z, 'type', 'Spearman', 'rows', 'complete');
    for j=1:num_time_points
    % Only look at first row because you are only interested in how the
    % parameter correlates with the time points, not the time points to
    % each other.
        prcc(i,j) = rho(1,j+1); % prcc = (num_param) x (num_timepoints)
        prcc_significance(i,j) = pvals(1,j+1); % prcc_sign = (num_param) x (num_timepoints)
    end
end
PRCCs=prcc'; % PRCCs = (num_timepoints) x (num_param)
uncorrected_significance=prcc_significance'; % uncorrected_significance = (num_timepoints) x (num_param)
prcc=PRCCs;
sign=uncorrected_significance;

%% Multiple tests correction: Bonferroni
%tests=length(s)*num_parameters; % # of tests performed
%correction_factor=tests;
%Bonf_sign=uncorrected_sign*tests;
%uncorrected_sign; % uncorrected p-value
%Bonf_sign;  % Bonferroni correction

sign_label=struct;
sign_label.uncorrected_significance=uncorrected_significance;

% figure;
for r=1:length(s)
    % Find all significant variables from PRCC calculation
    a=find(uncorrected_significance(r,:)<alpha);
%     PRCC_var(a);
%     prcc(r,a);
    b=num2cell(prcc(r,a));
    sign_label.index{r}=a;
    sign_label.label{r}=PRCC_var(a);
    sign_label.value{r}=b;
    %% Plots of PRCCs and their p-values
%     subplot(1,2,1),bar(PRCCs(r,:)),title(['PRCCs at time = ' num2str(s(r))]);%,set(gca,'XTickLabel',PRCC_var,'XTick',[1:num_parameters]),Title('PRCCs');
%     subplot(1,2,2),bar(uncorrected_significance(r,:)),...
%        set(gca,'YLim',[0,.1]),title('P values');%'XTickLabel',PRCC_var,'XTick',[1:num_parameters],
end % for r
end %function