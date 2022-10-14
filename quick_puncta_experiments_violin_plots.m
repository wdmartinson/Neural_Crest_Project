function quick_puncta_experiments_violin_plots(D_indices, D_labels, Scatter_Titles)
    D_max = length(D_indices);
    Lt = length(Scatter_Titles);
    Data = zeros(200, D_max, Lt);
    for kk = 1:D_max
    J = D_indices(kk);
    system(['tar -xf ',strcat('/Volumes/easystore/NCC_puncta_model_output_data_experiments/',num2str(J, '%04.f'),'.tar.gz'), ' summary_statistics.mat']);
    load('summary_statistics.mat', 'Summary_Statistics');
    Data(:, kk, 1) = Summary_Statistics.Max_X_Length_of_Stream;
    Data(:, kk, 2) = Summary_Statistics.Max_Y_Length_of_Stream;
    Data(:, kk, 3) = Summary_Statistics.Gap_Statistic_Standard_1_Error;
    Data(:, kk, 4) = Summary_Statistics.Persistence_Homology_Largest_Radius_Before_All_Connected;
    Data(:, kk, 5) = Summary_Statistics.Average_Min_Distance_per_Realization;
%     Data(:, kk, 5) = Summary_Statistics.Average_Min_Distance_per_Realization./15;
    Data(:, kk, 6) = Summary_Statistics.Average_Number_of_Cells_within_50_um_Ball_per_Realization;
    Data(:, kk, 7) = Summary_Statistics.Average_Speed;
    Data(:, kk, 8) = Summary_Statistics.Order_Parameter_Leaders;
    idx = find(all(diff(Data(:,:,8)) < 1e-16) == 1);
    Data(:, idx, 8) = zeros(size(Data, 1), length(idx));
    Data(:, kk, 9) = Summary_Statistics.Order_Parameter_Followers;
    Data(:, kk, 10) = Summary_Statistics.Percent_Target_Area_Covered;
    Data(:, kk, 11) = Summary_Statistics.Percent_Cell_Area_Outside_Target_Sites;
    if Lt > 11 && isfield(Summary_Statistics, 'Average_Leader_Speed_per_Realization')
        Data(:, kk, 12) = Summary_Statistics.Average_Leader_Speed_per_Realization;
    end
    if Lt > 11 && isfield(Summary_Statistics, 'Average_Follower_Speed_per_Realization')
        Data(:, kk, 13) = Summary_Statistics.Average_Follower_Speed_per_Realization;
    end
    if Lt > 11 && isfield(Summary_Statistics, 'Average_Min_Distance_from_Follower_to_Leader_per_Realization')
        Data(:, kk, 14) = Summary_Statistics.Average_Min_Distance_from_Follower_to_Leader_per_Realization;
%         Data(:, kk, 14) = Summary_Statistics.Average_Min_Distance_from_Follower_to_Leader_per_Realization./15;
    end
    if Lt > 11 && isfield(Summary_Statistics, 'Number_of_Outlier_Cells')
        Data(:, kk, 15) = Summary_Statistics.Number_of_Outlier_Cells;
    end
    if Lt > 11 && isfield(Summary_Statistics, 'Average_Leader_Speed_per_Realization')
        Data(:, kk, 16) = Summary_Statistics.Average_Leader_Speed_per_Realization - Summary_Statistics.Average_Follower_Speed_per_Realization;
    end
    if Lt > 11 && isfield(Summary_Statistics, 'Average_Number_of_Cells_within_50_um_Ball_per_Realization')
        Data(:, kk, 17) = Summary_Statistics.Average_Number_of_Cells_within_50_um_Ball_per_Realization./(Summary_Statistics.Number_of_Leaders + Summary_Statistics.Number_of_Followers);
    end
    if Lt > 11 && isfield(Summary_Statistics, 'Porosity_Measure')
        Data(:, kk, 18) = Summary_Statistics.Porosity_Measure;
    end
    system('rm -rf summary_statistics.mat');
    end % for kk
    
    % Compute Mann-Whitney U tests and KS tests for all pairwise 
    % comparisons, and determine signficance with the bonferroni correction
    % lambda = (#Pairwise Comparisons)*(#Output variables)
    alpha = 0.01; % Statistical significance cutoff
    u_vals = zeros(D_max, Lt);
    ks_vals = zeros(D_max, Lt);
    bonferroni_utest_pvals = zeros(D_max, Lt);
    bonferroni_kstest_pvals = zeros(D_max, Lt);
    for ll = 2:D_max
        for qq = 1:Lt
            [p,~,stats] = ranksum(Data(:,1,qq), Data(:, ll, qq), 'alpha', alpha);
            % Store statistic, raw p-value, and corrected p-value:
            u_vals(ll,qq) = stats.ranksum - length(Data(:,1,qq))*(length(Data(:,1,qq))+1)./2;
            bonferroni_utest_pvals(ll, qq) = (D_max-1)*Lt*p;
            bonf_alpha = alpha;
            [~,p,ks2stat] = kstest2(Data(:,1,qq), Data(:, ll, qq), 'alpha', bonf_alpha);
            ks_vals(ll,qq) = ks2stat;
            bonferroni_kstest_pvals(ll, qq) = p*Lt*(D_max-1);
        end % for qq
    end % for ll
    
    % Violin plot
    for i = 1:Lt
        figure('units', 'inches', 'position', [0,0,14,10]);
        Data(isinf(Data)) = 100;
        violinplot(Data(:, :, i));
        for ll = 2:D_max
            if (bonferroni_utest_pvals(ll, i) < alpha) && (bonferroni_kstest_pvals(ll, i) < alpha)
                ptr = plot(ll, 1.1*max(Data(:, ll, i)), '*k', 'markersize', 25);
%                 ptr = plot(ll, 20 + max(Data(:, ll, i)), '*k', 'markersize', 25);
            else
                ptr = plot(NaN, NaN, '*k');
            end
%             if bonferroni_utest_pvals(ll, i) < alpha  
%                 ptr = plot(ll+0.05, 1.1*max(Data(:, ll, i)), '*k', 'markersize', 15);
%             else
%                 ptr = plot(NaN, NaN, '*k');
%             end
%             if bonferroni_kstest_pvals(ll, i) < alpha
%                 ptr2 = plot(ll-0.05, 1.1*max(Data(:, ll, i)), '*r', 'markersize', 15);
%             else
%                 ptr2 = plot(NaN, NaN, '*r');
%             end
        end % for ll
        title(Scatter_Titles{i});
%         legend([ptr(1), ptr2(1)], 'Mann-Whitney Significant w/Bonferroni ($p<0.01$)', 'KS Signficant w/Bonferroni ($p<0.01$)', 'Location', 'best', 'interpreter', 'latex');
        legend(ptr(1), 'Mann-Whitney and KS Significant w/Bonferroni ($p<0.01$)', 'Location', 'best', 'interpreter', 'latex');
        set(gca, 'XTick', 1:D_max, 'XTickLabel', D_labels,'TickLabelInterpreter','latex','fontsize',30);
        box on;
    end
end