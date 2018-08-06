clear all;
close all;
pwd; 
currentFolder = pwd;
addpath(currentFolder, genpath('D:/Program Files/MATLAB/nearestSPD'));
import SLEP.*;
import nearestSPD.*;
import Spd_Mat.*;
import fixcorrmatrix.*;
run mexC.m;
A = importdata('SAS saved as spss.csv');
%save IC_SP.mat A;
%load('IC_SP.mat');
data_double = A.data;
data_name = A.textdata(1,5:end);
%var = char(data_name); 

% to pre-process whole data set
% to get missing ratios of all variables
for i = 1:size(data_double,2)
    zeros(1,i) = sum(isnan(data_double(:,i)))/size(data_double,1);
end

% to remove variables with 75% or more missing data and get new data set data_double_new with 1688 variables
k = 1;
for i = 1:size(data_double, 2)
        if zeros(1,i) < 0.75
            data_name_new(1, k) = data_name(1,i); % to get the subset of variable names without respondents name
            data_double_new(:, k) = data_double(:,i);
            %data_double_new(k).name = data_name(1, i);
            %data_double_new(k).value = data_double(:,i);
            k = k + 1;
        end
end

% to impute missing data using the nearest neighbor method
%knnimpute() is to find nearest column to impute, but here we want to use nearest row, then transpose first 
temp = transpose(data_double_new);
temp_imputed = knnimpute(temp);
data_imputed = transpose(temp_imputed);

%Feb-16-2017 new added
%generate bootstrap sampling data sets for data_imputed
tic;
cnt = 10;
[stat, index]=bootstrp(cnt,[],data_imputed);
bootstrp_data_imputed = data_imputed(index,:);

%then for each bootstrap sample, calculate optimal SICE matrix
for t = 1: cnt
    bootstrap_sample = bootstrp_data_imputed(1+(t-1)*197:197*t,:);
    % to standardize the variables
    data_s = zscore(bootstrap_sample);

    % to drop the last 6 columns(unknown meaning)
    data_s(:, 1683:end) = [];
    data_imputed(:, 1683:end) = []; % un-standardized data
    data_name_new(:, 1683:end) = [];

    % to subset data set by constructs-standardized
    for i = 1: length(data_name_new)
        construct(i). name = data_name_new{1,i}(4:end);
        construct(i). value = data_s(:,i);
    end

    % to subset data set by constructs-unstandardized
    for i = 1: length(data_name_new)
        construct_us(i). name = data_name_new{1,i}(4:end);
        construct_us(i). value = data_imputed(:,i);
    end

    % to sort the construct by names(mainly for form 125 sub-constructs for PCA)
    [unused, order] = sortrows(strvcat(construct(:).name));
    construct_sorted = construct(order); 

    k = 1;
    i = 1;
    j = 2;
    %subconstruct{1} = construct(1).value;
    %subconstruct_name{1} = construct(1).name;
    % for i = 1: length(data_name_new)
    %     for j = 1: length(data_name_new)
    %         if(strncmp(construct(i).name, construct(j).name, 5) == 1)
    %             subconstruct(k).name = construct(j).name;
    %             subconstruct(k).value = construct(j).value;
    %          end
    %         j = j + 1;
    %     end
    %     i = i +1;
    % end
    subconstruct(1).name = construct_sorted(1).name(1:5);
    subconstruct(1).value = construct_sorted(1).value;
    % to group data set into constructs by 4th~7th characters in variable names 
    while(i <= length(data_name_new) && j <= length(data_name_new))
        if(strncmp(construct_sorted(i).name, construct_sorted(j).name, 4) == 1)
            subconstruct(k).name = construct_sorted(j).name(1:4);
            subconstruct(k).value = [subconstruct(k).value, construct_sorted(j).value];
        else
            i = j;
            k = k + 1;
            subconstruct(k).name = construct_sorted(i).name(1:4);
            subconstruct(k).value = construct_sorted(i).value;
        end
        j = j + 1;
    end

    %********************************
    % to conduct PCA using default arguments
    %********************************
    for i = 1: length(subconstruct)
         [coeff{i},score{i},latent{i},tsquared{i},explained{i},mu{i}] = pca(subconstruct(i).value);
        % to get name of each construct
         pcadata_name(i,:) =  subconstruct(i).name;
    end
    % to get all the first PCs of 125 constructs and form a new dataset pcadata
    for i = 1: length(score)
        pcadata(:,i) = score{1,i}(:,1);
    end

    %*******************************************
    %group high and low sub-groups using criteron variables
    %*******************************************
    % for innovation dimension
    % to find column NRPSN10 and plot histogram
    %(here un-standardized data are used, same to other categories)
    for i = 1: length(construct_us)
        if(construct_us(i).name == 'nrpsn10')
            innovation = construct_us(i).value;
        end
    end
    % to group all records to high and low innovation categories (median, 
    %using unstandardized data,same to other categories)
    % at the same time group pcadata into 10 subsets
    m = 1;
    n = 1;
    for i = 1: length(innovation)
        if (innovation(i) >= median(innovation))
            %innovation_high(1,i) = construct_us(i).name;
            for j = 1:length(construct)
               innovation_high(m,j) = construct(j).value(i);
            end
             % at the same time group pcadata to high-innovation group
             %group observations to two pcadata groups(responding to high and low performance)
               pca_innovation_high(m,:) = pcadata(i, :);
               m = m + 1;
        else
            for j = 1:length(construct)
                 innovation_low(n,j) = construct(j).value(i);
            end
            pca_innovation_low(n,:) = pcadata(i, :);
            n = n + 1;
        end
    end


    % for quality dimension 
    % to find column GRCPN02 and plot histogram
    for i = 1:length(construct_us)
        if(construct_us(i).name == 'grcpn02')
            quality = construct_us(i).value;
        end
    end
    % to group all records to high and low quality categories (median)
    m = 1;
    n = 1;
    for i = 1: length(quality)
        if (quality(i) >= median(quality))
           for j = 1:length(construct)
               quality_high(m,j) = construct(j).value(i);
           end
            pca_quality_high(m,:) = pcadata(i, :);
            m = m + 1;
        else
            for j = 1:length(construct)
                 quality_low(n,j) = construct(j).value(i);
            end
             pca_quality_low(n,:) = pcadata(i, :);
            n = n + 1;
        end
    end      

    % for flexibility dimension
    % to find column GRCPN05 and plot histogram
    for i = 1:length(construct_us)
        if(construct_us(i).name == 'grcpn05')
            flexibility = construct_us(i).value;
        end
    end
    % to group all records to high and low flexibility categories (median)
    m = 1;
    n = 1;
    for i = 1: length(quality)
        if (flexibility(i) >= median(flexibility))
           for j = 1:length(construct)
               flexibility_high(m,j) = construct(j).value(i);
           end
            pca_flexibility_high(m,:) = pcadata(i, :);
            m = m + 1;
        else
            for j = 1:length(construct)
                 flexibility_low(n,j) = construct(j).value(i);
            end
            pca_flexibility_low(n,:) = pcadata(i, :);
            n = n + 1;
        end
    end       

    % for speed dimension
    % to find column GRCPN04  and plot histogram
    for i = 1:length(construct_us)
        if(construct_us(i).name == 'grcpn04')
            speed = construct_us(i).value;
        end
    end
    % to group all records to high and low speed categories (median)
    m = 1;
    n = 1;
    for i = 1: length(speed)
        if (speed(i) >= median(speed))
           for j = 1:length(construct)
               speed_high(m,j) = construct(j).value(i);
           end
            pca_speed_high(m,:) = pcadata(i, :);
            m = m + 1;
        else
            for j = 1:length(construct)
                 speed_low(n,j) = construct(j).value(i);
            end
            pca_speed_low(n,:) = pcadata(i, :);
            n = n + 1;
        end
    end       

    % for cost dimension
    % to find column GRADN00 and GRADN03   and plot histogram of GRADN00 divided by GRADN03 
    for i = 1:length(construct_us)
        if(construct_us(i).name == 'gradn00')
            cost_numerator = construct_us(i).value;
        end
    end
    for i = 1:length(construct_us)
        if(construct_us(i).name == 'gradn03')
            cost_denominator  = construct_us(i).value;
        end
    end

    for i = 1: length(cost_numerator)
        cost(i) = cost_numerator(i)/cost_denominator(i);
    end
    % to exclude the outlier(s)
    k = 1;
    for i = 1:length(cost_numerator)
        if(cost(i) <=20)
        cost_noutlier(k) = cost(i);
        k = k +1;
        end
    end
    % to group all records to high and low cost categories (median)
    m = 1;
    n = 1;
    for i = 1: length(cost_noutlier)
        if (cost_noutlier(i) >= median( cost_noutlier))
             for j = 1:length(construct)
                    cost_noutlier_high(m,j) =  construct(j).value(i);
             end
             pca_cost_high(m,:) = pcadata(i, :);
            m = m + 1;
        else
            for j = 1:length(construct)
                 cost_noutlier_low(n,j) =  construct(j).value(i);
            end
            pca_cost_low(n,:) = pcadata(i, :);
            n = n + 1;
        end
    end      
    %*******************************************
    %end of group high and low sub-groups using criteron variables
    %*******************************************


    %*********************************************************
    % Method 1: using BIC/AIC to determine optimized penalty paramters
    % Method 2: to keep arcs at 700, 400, 200 and compare the matrices to see
    % if there are statistical significance--at last we chose arc-200 to
    % demonstrate networks and central nodes
    %******************************************************************
    % SICE analysis based on principal components (first PC of each construct)-grouped
    %******************************************************************
    % to get covariance matrix of 10 subsets
    cov_pca_innovation_high = cov(pca_innovation_high);
    cov_pca_innovation_low = cov(pca_innovation_low);

    cov_pca_quality_high = cov(pca_quality_high);
    cov_pca_quality_low = cov(pca_quality_low);

    cov_pca_flexibility_high = cov(pca_flexibility_high);
    cov_pca_flexibility_low = cov(pca_flexibility_low);

    cov_pca_speed_high = cov(pca_speed_high);
    cov_pca_speed_low = cov(pca_speed_low);

    cov_pca_cost_high = cov(pca_cost_high);
    cov_pca_cost_low = cov(pca_cost_low);


    % Method 1: using BIC/AIC to determine optimized penalty paramters 0.1-10
    lamb = 0;
    for i = 1: 100
        lamb = i * 0.1;
        SICE_pca_innovation_high{i}  = sparseInverseCovariance(cov_pca_innovation_high, lamb, []);
        SICE_pca_innovation_low{i}  = sparseInverseCovariance(cov_pca_innovation_low, lamb, []);

        SICE_pca_quality_high{i}  = sparseInverseCovariance(cov_pca_quality_high, lamb, []);
        SICE_pca_quality_low{i}  = sparseInverseCovariance(cov_pca_quality_low, lamb, []);

        SICE_pca_flexibility_high{i}  = sparseInverseCovariance(cov_pca_flexibility_high, lamb, []);
        SICE_pca_flexibility_low{i}  = sparseInverseCovariance(cov_pca_flexibility_low, lamb, []);

        SICE_pca_speed_high{i}  = sparseInverseCovariance(cov_pca_speed_high, lamb, []);
        SICE_pca_speed_low{i}  = sparseInverseCovariance(cov_pca_speed_low, lamb, []);

        SICE_pca_cost_high{i}  = sparseInverseCovariance(cov_pca_cost_high, lamb, []);
        SICE_pca_cost_low{i}  = sparseInverseCovariance(cov_pca_cost_low, lamb, []);    
    end


    %to calculate BIC and AIC scores
    % to calculate Multivariate normal negative log-likelihood function
    sample_innovation_high = size(innovation_high,1);
    sample_innovation_low = size(innovation_low,1);

    sample_quality_high = size(quality_high,1);
    sample_quality_low = size(quality_low,1);

    sample_flexibility_high = size(flexibility_high,1);
    sample_flexibility_low = size(flexibility_low,1);

    sample_speed_high = size(speed_high,1);
    sample_speed_low = size(speed_low,1);

    sample_cost_high = size(cost_noutlier_high,1);
    sample_cost_low = size(cost_noutlier_low,1);

    for i = 1:100
        temp_matrix = SICE_pca_innovation_high{i};
        free_param_innovation_high(i) = (nnz(temp_matrix) - 125)/2; 
        mean_innovation_high = mean(pca_innovation_high);

      sum_rows = 0;
      for j = 1: sample_innovation_high
          sum_rows =  sum_rows + (pca_innovation_high(j,:) - mean_innovation_high) * temp_matrix * transpose(pca_innovation_high(j,:) - mean_innovation_high);
      end

       innovation_high_objective(i) = 0.5*sample_innovation_high* log(det(temp_matrix))  - 0.5 * sum_rows;

       if free_param_innovation_high(i) > 0
           [aic_innovation_high(i),bic_innovation_high(i)] = aicbic(innovation_high_objective(i),free_param_innovation_high(i),sample_innovation_high);
       else
           aic_innovation_high(i) = NaN;
           bic_innovation_high(i) = NaN;
       end   
    end
      
    %Feb-16-2017 new added
    %to get the smallest BIC score and optimal SICE matrix for
    %innovation_high subgroup
     [~, I] = min(bic_innovation_high);
    opt_SICE_innovation_high{t} = SICE_pca_innovation_high{I}; 


    for i = 1:100
        temp_matrix = SICE_pca_innovation_low{i};
         free_param_innovation_low(i)  = (nnz(temp_matrix) - 125)/2; 
         mean_innovation_low = mean(pca_innovation_low);

     sum_rows = 0;
      for j = 1: sample_innovation_low
          sum_rows =  sum_rows + (pca_innovation_low(j,:) - mean_innovation_low) * temp_matrix * transpose(pca_innovation_low(j,:) - mean_innovation_low);
      end

       innovation_low_objective(i) = 0.5*sample_innovation_low * log(det(temp_matrix))  - 0.5 * sum_rows; 
       if free_param_innovation_low(i) > 0
             [aic_innovation_low(i),bic_innovation_low(i)] = aicbic(innovation_low_objective(i),free_param_innovation_low(i),sample_innovation_low);
       else
           aic_innovation_low(i) = NaN;
           bic_innovation_low(i) = NaN;
       end

    end
     %Feb-16-2017 new added
    %to get the smallest BIC score and optimal SICE matrix for
    %innovation_low subgroup
     [~, I] = min(bic_innovation_low);
     opt_SICE_innovation_low{t} = SICE_pca_innovation_low{I}; 

    for i = 1:100
        temp_matrix = SICE_pca_quality_high{i};
         free_param_quality_high(i) = (nnz(temp_matrix) - 125)/2; 
        mean_quality_high = mean(pca_quality_high);

     sum_rows = 0;
      for j = 1: sample_quality_high
           sum_rows =  sum_rows + (pca_quality_high(j,:) - mean_quality_high) * temp_matrix * transpose(pca_quality_high(j,:) - mean_quality_high);
      end

       quality_high_objective(i) = 0.5*sample_quality_high * log(det(temp_matrix))  - 0.5 * sum_rows; 
       if free_param_quality_high(i) > 0
           [aic_quality_high(i),bic_quality_high(i)] = aicbic(quality_high_objective(i),free_param_quality_high(i),sample_quality_high);
       else
           aic_quality_high(i) = NaN;
           bic_quality_high(i) = NaN;
       end

    end

     %Feb-17-2017 new added
    %to get the smallest BIC score and optimal SICE matrix for
    %quality_high subgroup
     [~, I] = min(bic_quality_high);
     opt_SICE_quality_high{t} = SICE_pca_quality_high{I}; 


    for i = 1:100
        temp_matrix = SICE_pca_quality_low{i};
          free_param_quality_low(i) = (nnz(temp_matrix) - 125)/2; 
         mean_quality_low = mean(pca_quality_low);

      sum_rows = 0;
      for j = 1: sample_quality_low
          sum_rows =  sum_rows + (pca_quality_low(j,:) - mean_quality_low) * temp_matrix * transpose(pca_quality_low(j,:) - mean_quality_low);
      end

       quality_low_objective(i) = 0.5*sample_quality_low * log(det(temp_matrix))  - 0.5 * sum_rows;   
       if free_param_quality_low(i) > 0
            [aic_quality_low(i),bic_quality_low(i)] = aicbic(quality_low_objective(i),free_param_quality_low(i),sample_quality_low);
       else
           aic_quality_low(i) = NaN;
           bic_quality_low(i) = NaN;
       end

       end

    %Feb-17-2017 new added
    %to get the smallest BIC score and optimal SICE matrix for
    %quality_low subgroup
     [~, I] = min(bic_quality_low);
     opt_SICE_quality_low{t} = SICE_pca_quality_low{I}; 

    for i = 1:100
        temp_matrix = SICE_pca_flexibility_high{i};
          free_param_flexibility_high(i) =  (nnz(temp_matrix) - 125)/2; 
          mean_flexibility_high = mean(pca_flexibility_high);

    sum_rows = 0;
      for j = 1: sample_flexibility_high
         sum_rows =  sum_rows + (pca_flexibility_high(j,:) - mean_flexibility_high) * temp_matrix * transpose(pca_flexibility_high(j,:) - mean_flexibility_high);
      end

      flexibility_high_objective(i) = 0.5*sample_flexibility_high * log(det(temp_matrix))  - 0.5 * sum_rows;   
      if free_param_flexibility_high(i) > 0
          [aic_flexibility_high(i),bic_flexibility_high(i)] = aicbic(flexibility_high_objective(i),free_param_flexibility_high(i),sample_flexibility_high);
      else
           aic_flexibility_high(i) = NaN;
           bic_flexibility_high(i) = NaN;
      end

    end
    
    %Feb-27-2017 new added
    %to get the smallest BIC score and optimal SICE matrix for
    %flexibility_high subgroup
     [~, I] = min(bic_flexibility_high);
     opt_SICE_flexibility_high{t} = SICE_pca_flexibility_high{I}; 
     
    for i = 1:100
        temp_matrix = SICE_pca_flexibility_low{i};
         free_param_flexibility_low(i) = (nnz(temp_matrix) - 125)/2; 
        mean_flexibility_low = mean(pca_flexibility_low);

     sum_rows = 0;
      for j = 1: sample_flexibility_low
           sum_rows =  sum_rows + (pca_flexibility_low(j,:) - mean_flexibility_low) * temp_matrix * transpose(pca_flexibility_low(j,:) - mean_flexibility_low);
      end

      flexibility_low_objective(i) = 0.5*sample_flexibility_low * log(det(temp_matrix))  - 0.5 * sum_rows;     
      if free_param_flexibility_low(i) > 0
           [aic_flexibility_low(i),bic_flexibility_low(i)] = aicbic(flexibility_low_objective(i),free_param_flexibility_low(i),sample_flexibility_low);
       else
           aic_flexibility_low(i) = NaN;
           bic_flexibility_low(i) = NaN;
      end
    end
    
    %Feb-27-2017 new added
    %to get the smallest BIC score and optimal SICE matrix for
    %flexibility_low subgroup
     [~, I] = min(bic_flexibility_low);
     opt_SICE_flexibility_low{t} = SICE_pca_flexibility_low{I}; 


    for i = 1:100
        temp_matrix = SICE_pca_speed_high{i};
          free_param_speed_high(i) = (nnz(temp_matrix) - 125)/2; 
       mean_speed_high = mean(pca_speed_high);

    sum_rows = 0;
      for j = 1: sample_speed_high
           sum_rows =  sum_rows + (pca_speed_high(j,:) - mean_speed_high) * temp_matrix * transpose(pca_speed_high(j,:) - mean_speed_high);
      end

      speed_high_objective(i) = 0.5*sample_speed_high * log(det(temp_matrix))  - 0.5 * sum_rows;     
      if free_param_speed_high(i) > 0
          [aic_speed_high(i),bic_speed_high(i)] = aicbic(speed_high_objective(i),free_param_speed_high(i),sample_speed_high);
       else
           aic_speed_high(i) = NaN;
           bic_speed_high(i) = NaN;
      end
    end
    %Feb-27-2017 new added
    %to get the smallest BIC score and optimal SICE matrix for
    %speed_high subgroup
     [~, I] = min(bic_speed_high);
     opt_SICE_speed_high{t} = SICE_pca_speed_high{I}; 

    for i = 1:100
        temp_matrix = SICE_pca_speed_low{i};
         free_param_speed_low(i) =(nnz(temp_matrix) - 125)/2; 
        mean_speed_low = mean(pca_speed_low);

        sum_rows = 0;
      for j = 1: sample_speed_low
           sum_rows =  sum_rows + (pca_speed_low(j,:) - mean_speed_low) * temp_matrix * transpose(pca_speed_low(j,:) - mean_speed_low);
      end

      speed_low_objective(i) = 0.5*sample_speed_low * log(det(temp_matrix))  - 0.5 * sum_rows;    
      if free_param_speed_low(i) > 0
           [aic_speed_low(i),bic_speed_low(i)] = aicbic(speed_low_objective(i),free_param_speed_low(i),sample_speed_low);
      else
           aic_speed_low(i) = NaN;
           bic_speed_low(i) = NaN;
      end

    end
    
    %Feb-27-2017 new added
    %to get the smallest BIC score and optimal SICE matrix for
    %speed_low subgroup
     [~, I] = min(bic_speed_low);
     opt_SICE_speed_low{t} = SICE_pca_speed_low{I}; 



    for i = 1:100
        temp_matrix = SICE_pca_cost_high{i};
         free_param_cost_high(i) = (nnz(temp_matrix) - 125)/2; 
        mean_cost_high = mean(pca_cost_high);

      sum_rows = 0;
      for j = 1: sample_cost_high
         sum_rows =  sum_rows + (pca_cost_high(j,:) - mean_cost_high) * temp_matrix * transpose(pca_cost_high(j,:) - mean_cost_high);
      end

      cost_high_objective(i) = 0.5*sample_cost_high * log(det(temp_matrix))  - 0.5 * sum_rows;     
      if free_param_cost_high(i) > 0
          [aic_cost_high(i),bic_cost_high(i)] = aicbic(cost_high_objective(i),free_param_cost_high(i),sample_cost_high);
      else
           aic_cost_high(i) = NaN;
           bic_cost_high(i) = NaN;
      end
    end

    %Feb-27-2017 new added
    %to get the smallest BIC score and optimal SICE matrix for
    %cost_high subgroup
     [~, I] = min(bic_cost_high);
     opt_SICE_cost_high{t} = SICE_pca_cost_high{I}; 

    for i = 1:100
        temp_matrix = SICE_pca_cost_low{i};
         free_param_cost_low(i) = (nnz(temp_matrix) - 125)/2; 
         mean_cost_low = mean(pca_cost_low);

       sum_rows = 0;
      for j = 1: sample_cost_low
          sum_rows =  sum_rows + (pca_cost_low(j,:) - mean_cost_low) * temp_matrix * transpose(pca_cost_low(j,:) - mean_cost_low);
      end

      cost_low_objective(i) = 0.5*sample_cost_low * log(det(temp_matrix))  - 0.5 * sum_rows;      
      if free_param_cost_low(i) > 0
        [aic_cost_low(i),bic_cost_low(i)] = aicbic(cost_low_objective(i),free_param_cost_low(i),sample_cost_low);
     else
           aic_cost_low(i) = NaN;
           bic_cost_low(i) = NaN;
      end
    end
    
    %Feb-27-2017 new added
    %to get the smallest BIC score and optimal SICE matrix for
    %cost_low subgroup
     [~, I] = min(bic_cost_low);
     opt_SICE_cost_low{t} = SICE_pca_cost_low{I}; 

end
toc;

%Feb-16-2017 new added
%calculate mean and variance of innovation_high
mean_innovation_high = opt_SICE_innovation_high{1};
for i = 2:cnt
	mean_innovation_high = mean_innovation_high + opt_SICE_innovation_high{i};
end
mean_innovation_high = mean_innovation_high/cnt;

var_innovation_high = (opt_SICE_innovation_high{1} - mean_innovation_high).^2;
for i = 2:cnt
	var_innovation_high = var_innovation_high + (opt_SICE_innovation_high{i} - mean_innovation_high).^2;
end
var_innovation_high = var_innovation_high/cnt;

%cakculate the CI with 5% confidence level of mean
%kiwer limit and upper limit
LCL_innvation_high = mean_innovation_high - 1.96*sqrt(var_innovation_high)/sqrt(cnt);
UCL_innvation_high = mean_innovation_high + 1.96*sqrt(var_innovation_high)/sqrt(cnt);
%calculate multiplication of UCL_innvation_high and LCL_innvation_high
%element-wise, then negative elements repond to CI containing 0
indicator_innovation_high = LCL_innvation_high.*UCL_innvation_high;

fiiltered_mean_innovation_high = mean_innovation_high;
for i = 1:125
    for j = 1:125
        if indicator_innovation_high(i, j) < 0
            fiiltered_mean_innovation_high(i,j) = 0;
        end
    end
end

%calculate mean and variance of innovation_low
mean_innovation_low = opt_SICE_innovation_low{1};
for i = 2:cnt
	mean_innovation_low = mean_innovation_low + opt_SICE_innovation_low{i};
end
mean_innovation_low = mean_innovation_low/cnt;

var_innovation_low = (opt_SICE_innovation_low{1} - mean_innovation_low).^2;
for i = 2:cnt
	var_innovation_low = var_innovation_low + (opt_SICE_innovation_low{i} - mean_innovation_low).^2;
end
var_innovation_low = var_innovation_low/cnt;

%cakculate the CI with 5% confidence level of mean
%kiwer limit and upper limit
LCL_innovation_low = mean_innovation_low - 1.96*sqrt(var_innovation_low)/sqrt(cnt);
UCL_innovation_low = mean_innovation_low + 1.96*sqrt(var_innovation_low)/sqrt(cnt);
%calculate multiplication of UCL_innovation_low and LCL_innovation_low
%element-wise, then negative elements repond to CI containing 0
indicator_innovation_low = LCL_innovation_low.*UCL_innovation_low;

fiiltered_mean_innovation_low = mean_innovation_low;
for i = 1:125
    for j = 1:125
        if indicator_innovation_low(i, j) < 0
            fiiltered_mean_innovation_low(i,j) = 0;
        end
    end
end

%Feb-17-2017 new added
%calculate mean and variance of quality_high
mean_quality_high = opt_SICE_quality_high{1};
for i = 2:cnt
	mean_quality_high = mean_quality_high + opt_SICE_quality_high{i};
end
mean_quality_high = mean_quality_high/cnt;

var_quality_high = (opt_SICE_quality_high{1} - mean_quality_high).^2;
for i = 2:cnt
	var_quality_high = var_quality_high + (opt_SICE_quality_high{i} - mean_quality_high).^2;
end
var_quality_high = var_quality_high/cnt;

%cakculate the CI with 5% confidence level of mean
%kiwer limit and upper limit
LCL_quality_high = mean_quality_high - 1.96*sqrt(var_quality_high)/sqrt(cnt);
UCL_quality_high = mean_quality_high + 1.96*sqrt(var_quality_high)/sqrt(cnt);
%calculate multiplication of UCL_quality_high and LCL_quality_high
%element-wise, then negative elements repond to CI containing 0
indicator_quality_high = LCL_quality_high.*UCL_quality_high;

fiiltered_mean_quality_high = mean_quality_high;
for i = 1:125
    for j = 1:125
        if indicator_quality_high(i, j) < 0
            fiiltered_mean_quality_high(i,j) = 0;
        end
    end
end

%calculate mean and variance of quality_low
mean_quality_low = opt_SICE_quality_low{1};
for i = 2:cnt
	mean_quality_low = mean_quality_low + opt_SICE_quality_low{i};
end
mean_quality_low = mean_quality_low/cnt;

var_quality_low = (opt_SICE_quality_low{1} - mean_quality_low).^2;
for i = 2:cnt
	var_quality_low = var_quality_low + (opt_SICE_quality_low{i} - mean_quality_low).^2;
end
var_quality_low = var_quality_low/cnt;

%cakculate the CI with 5% confidence level of mean
%kiwer limit and upper limit
LCL_quality_low = mean_quality_low - 1.96*sqrt(var_quality_low)/sqrt(cnt);
UCL_quality_low = mean_quality_low + 1.96*sqrt(var_quality_low)/sqrt(cnt);
%calculate multiplication of UCL_quality_low and LCL_quality_low
%element-wise, then negative elements repond to CI containing 0
indicator_quality_low = LCL_quality_low.*UCL_quality_low;

fiiltered_mean_quality_low = mean_quality_low;
for i = 1:125
    for j = 1:125
        if indicator_quality_low(i, j) < 0
            fiiltered_mean_quality_low(i,j) = 0;
        end
    end
end

%Feb-27-2017 new added
%calculate mean and variance of flexibility_high
mean_flexibility_high = opt_SICE_flexibility_high{1};
for i = 2:cnt
	mean_flexibility_high = mean_flexibility_high + opt_SICE_flexibility_high{i};
end
mean_flexibility_high = mean_flexibility_high/cnt;

var_flexibility_high = (opt_SICE_flexibility_high{1} - mean_flexibility_high).^2;
for i = 2:cnt
	var_flexibility_high = var_flexibility_high + (opt_SICE_flexibility_high{i} - mean_flexibility_high).^2;
end
var_flexibility_high = var_flexibility_high/cnt;

%cakculate the CI with 5% confidence level of mean
%kiwer limit and upper limit
LCL_flexibility_high = mean_flexibility_high - 1.96*sqrt(var_flexibility_high)/sqrt(cnt);
UCL_flexibility_high = mean_flexibility_high + 1.96*sqrt(var_flexibility_high)/sqrt(cnt);
%calculate multiplication of UCL_flexibility_high and LCL_flexibility_high
%element-wise, then negative elements repond to CI containing 0
indicator_flexibility_high = LCL_flexibility_high.*UCL_flexibility_high;

fiiltered_mean_flexibility_high = mean_flexibility_high;
for i = 1:125
    for j = 1:125
        if indicator_flexibility_high(i, j) < 0
            fiiltered_mean_flexibility_high(i,j) = 0;
        end
    end
end

%calculate mean and variance of flexibility_low
mean_flexibility_low = opt_SICE_flexibility_low{1};
for i = 2:cnt
	mean_flexibility_low = mean_flexibility_low + opt_SICE_flexibility_low{i};
end
mean_flexibility_low = mean_flexibility_low/cnt;

var_flexibility_low = (opt_SICE_flexibility_low{1} - mean_flexibility_low).^2;
for i = 2:cnt
	var_flexibility_low = var_flexibility_low + (opt_SICE_flexibility_low{i} - mean_flexibility_low).^2;
end
var_flexibility_low = var_flexibility_low/cnt;

%cakculate the CI with 5% confidence level of mean
%kiwer limit and upper limit
LCL_flexibility_low = mean_flexibility_low - 1.96*sqrt(var_flexibility_low)/sqrt(cnt);
UCL_flexibility_low = mean_flexibility_low + 1.96*sqrt(var_flexibility_low)/sqrt(cnt);
%calculate multiplication of UCL_flexibility_low and LCL_flexibility_low
%element-wise, then negative elements repond to CI containing 0
indicator_flexibility_low = LCL_flexibility_low.*UCL_flexibility_low;

fiiltered_mean_flexibility_low = mean_flexibility_low;
for i = 1:125
    for j = 1:125
        if indicator_flexibility_low(i, j) < 0
            fiiltered_mean_flexibility_low(i,j) = 0;
        end
    end
end

%Feb-27-2017 new added
%calculate mean and variance of speed_high
mean_speed_high = opt_SICE_speed_high{1};
for i = 2:cnt
	mean_speed_high = mean_speed_high + opt_SICE_speed_high{i};
end
mean_speed_high = mean_speed_high/cnt;

var_speed_high = (opt_SICE_speed_high{1} - mean_speed_high).^2;
for i = 2:cnt
	var_speed_high = var_speed_high + (opt_SICE_speed_high{i} - mean_speed_high).^2;
end
var_speed_high = var_speed_high/cnt;

%cakculate the CI with 5% confidence level of mean
%kiwer limit and upper limit
LCL_speed_high = mean_speed_high - 1.96*sqrt(var_speed_high)/sqrt(cnt);
UCL_speed_high = mean_speed_high + 1.96*sqrt(var_speed_high)/sqrt(cnt);
%calculate multiplication of UCL_speed_high and LCL_speed_high
%element-wise, then negative elements repond to CI containing 0
indicator_speed_high = LCL_speed_high.*UCL_speed_high;

fiiltered_mean_speed_high = mean_speed_high;
for i = 1:125
    for j = 1:125
        if indicator_speed_high(i, j) < 0
            fiiltered_mean_speed_high(i,j) = 0;
        end
    end
end

%calculate mean and variance of speed_low
mean_speed_low = opt_SICE_speed_low{1};
for i = 2:cnt
	mean_speed_low = mean_speed_low + opt_SICE_speed_low{i};
end
mean_speed_low = mean_speed_low/cnt;

var_speed_low = (opt_SICE_speed_low{1} - mean_speed_low).^2;
for i = 2:cnt
	var_speed_low = var_speed_low + (opt_SICE_speed_low{i} - mean_speed_low).^2;
end
var_speed_low = var_speed_low/cnt;

%cakculate the CI with 5% confidence level of mean
%kiwer limit and upper limit
LCL_speed_low = mean_speed_low - 1.96*sqrt(var_speed_low)/sqrt(cnt);
UCL_speed_low = mean_speed_low + 1.96*sqrt(var_speed_low)/sqrt(cnt);
%calculate multiplication of UCL_speed_low and LCL_speed_low
%element-wise, then negative elements repond to CI containing 0
indicator_speed_low = LCL_speed_low.*UCL_speed_low;

fiiltered_mean_speed_low = mean_speed_low;
for i = 1:125
    for j = 1:125
        if indicator_speed_low(i, j) < 0
            fiiltered_mean_speed_low(i,j) = 0;
        end
    end
end

%Feb-27-2017 new added
%calculate mean and variance of cost_high
mean_cost_high = opt_SICE_cost_high{1};
for i = 2:cnt
	mean_cost_high = mean_cost_high + opt_SICE_cost_high{i};
end
mean_cost_high = mean_cost_high/cnt;

var_cost_high = (opt_SICE_cost_high{1} - mean_cost_high).^2;
for i = 2:cnt
	var_cost_high = var_cost_high + (opt_SICE_cost_high{i} - mean_cost_high).^2;
end
var_cost_high = var_cost_high/cnt;

%cakculate the CI with 5% confidence level of mean
%kiwer limit and upper limit
LCL_cost_high = mean_cost_high - 1.96*sqrt(var_cost_high)/sqrt(cnt);
UCL_cost_high = mean_cost_high + 1.96*sqrt(var_cost_high)/sqrt(cnt);
%calculate multiplication of UCL_cost_high and LCL_cost_high
%element-wise, then negative elements repond to CI containing 0
indicator_cost_high = LCL_cost_high.*UCL_cost_high;

fiiltered_mean_cost_high = mean_cost_high;
for i = 1:125
    for j = 1:125
        if indicator_cost_high(i, j) < 0
            fiiltered_mean_cost_high(i,j) = 0;
        end
    end
end

%calculate mean and variance of cost_low
mean_cost_low = opt_SICE_cost_low{1};
for i = 2:cnt
	mean_cost_low = mean_cost_low + opt_SICE_cost_low{i};
end
mean_cost_low = mean_cost_low/cnt;

var_cost_low = (opt_SICE_cost_low{1} - mean_cost_low).^2;
for i = 2:cnt
	var_cost_low = var_cost_low + (opt_SICE_cost_low{i} - mean_cost_low).^2;
end
var_cost_low = var_cost_low/cnt;

%cakculate the CI with 5% confidence level of mean
%kiwer limit and upper limit
LCL_cost_low = mean_cost_low - 1.96*sqrt(var_cost_low)/sqrt(cnt);
UCL_cost_low = mean_cost_low + 1.96*sqrt(var_cost_low)/sqrt(cnt);
%calculate multiplication of UCL_cost_low and LCL_cost_low
%element-wise, then negative elements repond to CI containing 0
indicator_cost_low = LCL_cost_low.*UCL_cost_low;

fiiltered_mean_cost_low = mean_cost_low;
for i = 1:125
    for j = 1:125
        if indicator_cost_low(i, j) < 0
            fiiltered_mean_cost_low(i,j) = 0;
        end
    end
end

% Jan-03-2018 new added
% to plot the networks from second methods
%BGobj_innovation_high = biograph(fiiltered_mean_innovation_high);
%view(BGobj_innovation_high);


%Feb-16-2017 new added
% plot the SICE matrix of innovation
figure;
subplot(1,2,1);
spy(fiiltered_mean_innovation_high);
title('Optimal IC Matrix of High-innovation');
subplot(1,2,2);
spy(fiiltered_mean_innovation_low);
title('Optimal IC Matrix of Low-innovation');

%{
%to plot the histograms of optimal IC Matrix (nonzero values)--innovation
%first exclude diagonal elements
no_diagonal_mean_innovation_high = fiiltered_mean_innovation_high;
for i = 1:125
    for j = 1:125
        if i == j
            no_diagonal_mean_innovation_high(i,j) = 0;
        end
    end
end

no_diagonal_mean_innovation_low = fiiltered_mean_innovation_low;
for i = 1:125
    for j = 1:125
        if i == j
            no_diagonal_mean_innovation_low(i,j) = 0;
        end
    end
end
%to find nonzero values 
[row_innovation_high_opt,col_innovation_high_opt,v_innovation_high_opt] = find(no_diagonal_mean_innovation_high);
[row_innovation_low_opt,col_innovation_low_opt,v_innovation_low_opt] = find(no_diagonal_mean_innovation_low);
%plot histograms
figure;
subplot(1,2,1);
hist(v_innovation_high_opt);
title('Nonzero values distribution of high-innovation IC Matrix');
subplot(1,2,2);
hist(v_innovation_low_opt);
title('Nonzero values distribution of low-innovation IC Matrix');
%}

%Feb-17-2017 new added
% plot the SICE matrix of quality
figure;
subplot(1,2,1);
spy(fiiltered_mean_quality_high);
title('Optimal IC Matrix of High-quality');
subplot(1,2,2);
spy(fiiltered_mean_quality_low);
title('Optimal IC Matrix of Low-quality');

%{
%to plot the histograms of optimal IC Matrix (nonzero values)--quality
%first exclude diagonal elements
no_diagonal_mean_quality_high = fiiltered_mean_quality_high;
for i = 1:125
    for j = 1:125
        if i == j
            no_diagonal_mean_quality_high(i,j) = 0;
        end
    end
end

no_diagonal_mean_quality_low = fiiltered_mean_quality_low;
for i = 1:125
    for j = 1:125
        if i == j
            no_diagonal_mean_quality_low(i,j) = 0;
        end
    end
end
%to find nonzero values 
[row_quality_high_opt,col_quality_high_opt,v_quality_high_opt] = find(no_diagonal_mean_quality_high);
[row_quality_low_opt,col_quality_low_opt,v_quality_low_opt] = find(no_diagonal_mean_quality_low);
%plot histograms
figure;
subplot(1,2,1);
hist(v_quality_high_opt);
title('Nonzero values distribution of high-quality IC Matrix');
subplot(1,2,2);
hist(v_quality_low_opt);
title('Nonzero values distribution of low-quality IC Matrix');
%}

%Feb-27-2017 new added
% plot the SICE matrix of flexibility
figure;
subplot(1,2,1);
spy(fiiltered_mean_flexibility_high);
title('Optimal IC Matrix of High-flexibility');
subplot(1,2,2);
spy(fiiltered_mean_flexibility_low);
title('Optimal IC Matrix of Low-flexibility');

%{
%to plot the histograms of optimal IC Matrix (nonzero values)--flexibility
%first exclude diagonal elements
no_diagonal_mean_flexibility_high = fiiltered_mean_flexibility_high;
for i = 1:125
    for j = 1:125
        if i == j
            no_diagonal_mean_flexibility_high(i,j) = 0;
        end
    end
end

no_diagonal_mean_flexibility_low = fiiltered_mean_flexibility_low;
for i = 1:125
    for j = 1:125
        if i == j
            no_diagonal_mean_flexibility_low(i,j) = 0;
        end
    end
end
%to find nonzero values 
[row_flexibility_high_opt,col_flexibility_high_opt,v_flexibility_high_opt] = find(no_diagonal_mean_flexibility_high);
[row_flexibility_low_opt,col_flexibility_low_opt,v_flexibility_low_opt] = find(no_diagonal_mean_flexibility_low);
%plot histograms
figure;
subplot(1,2,1);
hist(v_flexibility_high_opt);
title('Nonzero values distribution of high-flexibility IC Matrix');
subplot(1,2,2);
hist(v_flexibility_low_opt);
title('Nonzero values distribution of low-flexibility IC Matrix');
%}

%Feb-27-2017 new added
% plot the SICE matrix of speed
figure;
subplot(1,2,1);
spy(fiiltered_mean_speed_high);
title('Optimal IC Matrix of High-speed');
subplot(1,2,2);
spy(fiiltered_mean_speed_low);
title('Optimal IC Matrix of Low-speed');

%{
%to plot the histograms of optimal IC Matrix (nonzero values)--speed
%first exclude diagonal elements
no_diagonal_mean_speed_high = fiiltered_mean_speed_high;
for i = 1:125
    for j = 1:125
        if i == j
            no_diagonal_mean_speed_high(i,j) = 0;
        end
    end
end

no_diagonal_mean_speed_low = fiiltered_mean_speed_low;
for i = 1:125
    for j = 1:125
        if i == j
            no_diagonal_mean_speed_low(i,j) = 0;
        end
    end
end
%to find nonzero values 
[row_speed_high_opt,col_speed_high_opt,v_speed_high_opt] = find(no_diagonal_mean_speed_high);
[row_speed_low_opt,col_speed_low_opt,v_speed_low_opt] = find(no_diagonal_mean_speed_low);
%plot histograms
figure;
subplot(1,2,1);
hist(v_speed_high_opt);
title('Nonzero values distribution of high-speed IC Matrix');
subplot(1,2,2);
hist(v_speed_low_opt);
title('Nonzero values distribution of low-speed IC Matrix');
%}

%Feb-27-2017 new added
% plot the SICE matrix of cost
figure;
subplot(1,2,1);
spy(fiiltered_mean_cost_high);
title('Optimal IC Matrix of High-cost');
subplot(1,2,2);
spy(fiiltered_mean_cost_low);
title('Optimal IC Matrix of Low-cost');

toc;

%{
%to plot the histograms of optimal IC Matrix (nonzero values)--cost
%first exclude diagonal elements
no_diagonal_mean_cost_high = fiiltered_mean_cost_high;
for i = 1:125
    for j = 1:125
        if i == j
            no_diagonal_mean_cost_high(i,j) = 0;
        end
    end
end

no_diagonal_mean_cost_low = fiiltered_mean_cost_low;
for i = 1:125
    for j = 1:125
        if i == j
            no_diagonal_mean_cost_low(i,j) = 0;
        end
    end
end
%to find nonzero values 
[row_cost_high_opt,col_cost_high_opt,v_cost_high_opt] = find(no_diagonal_mean_cost_high);
[row_cost_low_opt,col_cost_low_opt,v_cost_low_opt] = find(no_diagonal_mean_cost_low);
%plot histograms
figure;
subplot(1,2,1);
hist(v_cost_high_opt);
title('Nonzero values distribution of high-cost IC Matrix');
subplot(1,2,2);
hist(v_cost_low_opt);
title('Nonzero values distribution of low-cost IC Matrix');
%}

%{

%****************************************************************
% Method2: to keep arcs at 700.400,200, run IC analysis for lambda varying from 1 to
% 10 with the interval of 0.01
%*****************************************************************
% lambda = 0 % the penalty parameter
% for i = 101:1000
%     lambda = i * 0.01;
%     pp(i) = i * 0.01;
%     SICE_pca_innovation_high{i}  = sparseInverseCovariance(cov_pca_innovation_high, lambda, []);
%     SICE_pca_innovation_high_nz(i) = nnz(SICE_pca_innovation_high{i}); % to get amounts of nonzero values
%    
%      SICE_pca_innovation_low{i}  = sparseInverseCovariance(cov_pca_innovation_low, lambda, []);
%     SICE_pca_innovation_low_nz(i) = nnz(SICE_pca_innovation_low{i});
%     % to find IC matrix with arcs of 700,400,200
%     if (SICE_pca_innovation_high_nz(i) < 705 & SICE_pca_innovation_high_nz(i) > 695)
%         SICE_pca_innovation_high_700 = SICE_pca_innovation_high{i};
%     end
%     if (SICE_pca_innovation_high_nz(i) < 405 & SICE_pca_innovation_high_nz(i) > 395)
%         SICE_pca_innovation_high_400 = SICE_pca_innovation_high{i};
%     end
%     if (SICE_pca_innovation_high_nz(i) < 205 & SICE_pca_innovation_high_nz(i) > 195)
%         SICE_pca_innovation_high_200 = SICE_pca_innovation_high{i};
%     end
%    
%     if (SICE_pca_innovation_low_nz(i) < 705 & SICE_pca_innovation_low_nz(i) > 695)
%         SICE_pca_innovation_low_700 = SICE_pca_innovation_low{i};
%     end
%     if (SICE_pca_innovation_low_nz(i) < 405 & SICE_pca_innovation_low_nz(i) > 395)
%         SICE_pca_innovation_low_400 = SICE_pca_innovation_low{i};
%     end
%     if (SICE_pca_innovation_low_nz(i) < 205 & SICE_pca_innovation_low_nz(i) > 195)
%         SICE_pca_innovation_low_200 = SICE_pca_innovation_low{i};
%     end
%    
%     SICE_pca_quality_high{i}  = sparseInverseCovariance(cov_pca_quality_high, lambda, []);
%     SICE_pca_quality_high_nz(i) = nnz(SICE_pca_quality_high{i});
%     SICE_pca_quality_low{i}  = sparseInverseCovariance(cov_pca_quality_low, lambda, []);
%     SICE_pca_quality_low_nz(i) = nnz(SICE_pca_quality_low{i});
%     % to find IC matrix with arcs of 700,400,200
%     if (SICE_pca_quality_high_nz(i) < 705 & SICE_pca_quality_high_nz(i) > 695)
%         SICE_pca_quality_high_700 = SICE_pca_quality_high{i};
%     end
%     if (SICE_pca_quality_high_nz(i) < 405 & SICE_pca_quality_high_nz(i) > 395)
%         SICE_pca_quality_high_400 = SICE_pca_quality_high{i};
%     end
%     if (SICE_pca_quality_high_nz(i) < 205 & SICE_pca_quality_high_nz(i) > 195)
%         SICE_pca_quality_high_200 = SICE_pca_quality_high{i};
%     end
%    
%     if (SICE_pca_quality_low_nz(i) < 705 & SICE_pca_quality_low_nz(i) > 695)
%         SICE_pca_quality_low_700 = SICE_pca_quality_low{i};
%     end
%     if (SICE_pca_quality_low_nz(i) < 405 & SICE_pca_quality_low_nz(i) > 395)
%         SICE_pca_quality_low_400 = SICE_pca_quality_low{i};
%     end
%     if (SICE_pca_quality_low_nz(i) < 205 & SICE_pca_quality_low_nz(i) > 195)
%         SICE_pca_quality_low_200 = SICE_pca_quality_low{i};
%     end
% 
%     SICE_pca_flexibility_high{i}  = sparseInverseCovariance(cov_pca_flexibility_high, lambda, []);
%     SICE_pca_flexibility_high_nz(i) = nnz(SICE_pca_flexibility_high{i});
%     SICE_pca_flexibility_low{i}  = sparseInverseCovariance(cov_pca_flexibility_low, lambda, []);
%     SICE_pca_flexibility_low_nz(i) = nnz(SICE_pca_flexibility_low{i});
%    
%     % to find IC matrix with arcs of 700,400,200
%     if (SICE_pca_flexibility_high_nz(i) < 705 & SICE_pca_flexibility_high_nz(i) > 695)
%         SICE_pca_flexibility_high_700 = SICE_pca_flexibility_high{i};
%     end
%     if (SICE_pca_flexibility_high_nz(i) < 405 & SICE_pca_flexibility_high_nz(i) > 395)
%         SICE_pca_flexibility_high_400 = SICE_pca_flexibility_high{i};
%     end
%     if (SICE_pca_flexibility_high_nz(i) < 205 & SICE_pca_flexibility_high_nz(i) > 195)
%         SICE_pca_flexibility_high_200 = SICE_pca_flexibility_high{i};
%     end
%    
%     if (SICE_pca_flexibility_low_nz(i) < 705 & SICE_pca_flexibility_low_nz(i) > 695)
%         SICE_pca_flexibility_low_700 = SICE_pca_flexibility_low{i};
%     end
%     if (SICE_pca_flexibility_low_nz(i) < 405 & SICE_pca_flexibility_low_nz(i) > 395)
%         SICE_pca_flexibility_low_400 = SICE_pca_flexibility_low{i};
%     end
%     if (SICE_pca_flexibility_low_nz(i) < 205 & SICE_pca_flexibility_low_nz(i) > 195)
%         SICE_pca_flexibility_low_200 = SICE_pca_flexibility_low{i};
%     end
% 
%     SICE_pca_speed_high{i}  = sparseInverseCovariance(cov_pca_speed_high, lambda, []);
%     SICE_pca_speed_high_nz(i) = nnz(SICE_pca_speed_high{i});
%     SICE_pca_speed_low{i}  = sparseInverseCovariance(cov_pca_speed_low, lambda, []);
%     SICE_pca_speed_low_nz(i) = nnz(SICE_pca_speed_low{i});
%    
%     % to find IC matrix with arcs of 700,400,200
%     if (SICE_pca_speed_high_nz(i) < 705 & SICE_pca_speed_high_nz(i) > 695)
%         SICE_pca_speed_high_700 = SICE_pca_speed_high{i};
%     end
%     if (SICE_pca_speed_high_nz(i) < 405 & SICE_pca_speed_high_nz(i) > 395)
%         SICE_pca_speed_high_400 = SICE_pca_speed_high{i};
%     end
%     if (SICE_pca_speed_high_nz(i) < 205 & SICE_pca_speed_high_nz(i) > 195)
%         SICE_pca_speed_high_200 = SICE_pca_speed_high{i};
%     end
%    
%     if (SICE_pca_speed_low_nz(i) < 705 & SICE_pca_speed_low_nz(i) > 695)
%         SICE_pca_speed_low_700 = SICE_pca_speed_low{i};
%     end
%     if (SICE_pca_speed_low_nz(i) < 405 & SICE_pca_speed_low_nz(i) > 395)
%         SICE_pca_speed_low_400 = SICE_pca_speed_low{i};
%     end
%     if (SICE_pca_speed_low_nz(i) < 205 & SICE_pca_speed_low_nz(i) > 195)
%         SICE_pca_speed_low_200 = SICE_pca_speed_low{i};
%     end
% 
%     SICE_pca_cost_high{i}  = sparseInverseCovariance(cov_pca_cost_high, lambda, []);
%     SICE_pca_cost_high_nz(i) = nnz(SICE_pca_cost_high{i});
%     SICE_pca_cost_low{i}  = sparseInverseCovariance(cov_pca_cost_low, lambda, []);
%      SICE_pca_cost_low_nz(i) = nnz(SICE_pca_cost_low{i});
%     
%      % to find IC matrix with arcs of 700,400,200
%     if (SICE_pca_cost_high_nz(i) < 705 & SICE_pca_cost_high_nz(i) > 695)
%         SICE_pca_cost_high_700 = SICE_pca_cost_high{i};
%     end
%     if (SICE_pca_cost_high_nz(i) < 405 & SICE_pca_cost_high_nz(i) > 395)
%         SICE_pca_cost_high_400 = SICE_pca_cost_high{i};
%     end
%     if (SICE_pca_cost_high_nz(i) < 205 & SICE_pca_cost_high_nz(i) > 195)
%         SICE_pca_cost_high_200 = SICE_pca_cost_high{i};
%     end
%    
%     if (SICE_pca_cost_low_nz(i) < 705 & SICE_pca_cost_low_nz(i) > 695 )
%         SICE_pca_cost_low_700 = SICE_pca_cost_low{i};
%     end
%     if (SICE_pca_cost_low_nz(i) < 405 & SICE_pca_cost_low_nz(i) > 395)
%         SICE_pca_cost_low_400 = SICE_pca_cost_low{i};
%     end
%     if (SICE_pca_cost_low_nz(i) < 205 & SICE_pca_cost_low_nz(i) > 195)
%         SICE_pca_cost_low_200 = SICE_pca_cost_low{i};
%     end
% end



%to plot the histograms of 700-,400-,200-arc IC Matrix (nonzero values)
%to find nonzero values
% [row_innovation_high_700,col_innovation_high_700,v_innovation_high_700] = find(SICE_pca_innovation_high_700);
% [row_innovation_low_700,col_innovation_low_700,v_innovation_low_700] = find(SICE_pca_innovation_low_700);
% [row_innovation_high_400,col_innovation_high_400,v_innovation_high_400] = find(SICE_pca_innovation_high_400);
% [row_innovation_low_400,col_innovation_low_400,v_innovation_low_400] = find(SICE_pca_innovation_low_400);
% [row_innovation_high_200,col_innovation_high_200,v_innovation_high_200] = find(SICE_pca_innovation_high_200);
% [row_innovation_low_200,col_innovation_low_200,v_innovation_low_200] = find(SICE_pca_innovation_low_200);
%
% [row_quality_high_700,col_quality_high_700,v_quality_high_700] = find(SICE_pca_quality_high_700);
% [row_quality_low_700,col_quality_low_700,v_quality_low_700] = find(SICE_pca_quality_low_700);
% [row_quality_high_400,col_quality_high_400,v_quality_high_400] = find(SICE_pca_quality_high_400);
% [row_quality_low_400,col_quality_low_400,v_quality_low_400] = find(SICE_pca_quality_low_400);
% [row_quality_high_200,col_quality_high_200,v_quality_high_200] = find(SICE_pca_quality_high_200);
% [row_quality_low_200,col_quality_low_200,v_quality_low_200] = find(SICE_pca_quality_low_200);
%
% [row_flexibility_high_700,col_flexibility_high_700,v_flexibility_high_700] = find(SICE_pca_flexibility_high_700);
% [row_flexibility_low_700,col_flexibility_low_700,v_flexibility_low_700] = find(SICE_pca_flexibility_low_700);
% [row_flexibility_high_400,col_flexibility_high_400,v_flexibility_high_400] = find(SICE_pca_flexibility_high_400);
% [row_flexibility_low_400,col_flexibility_low_400,v_flexibility_low_400] = find(SICE_pca_flexibility_low_400);
% [row_flexibility_high_200,col_flexibility_high_200,v_flexibility_high_200] = find(SICE_pca_flexibility_high_200);
% [row_flexibility_low_200,col_flexibility_low_200,v_flexibility_low_200] = find(SICE_pca_flexibility_low_200);
%
% [row_speed_high_700,col_speed_high_700,v_speed_high_700] = find(SICE_pca_speed_high_700);
% [row_speed_low_700,col_speed_low_700,v_speed_low_700] = find(SICE_pca_speed_low_700);
% [row_speed_high_400,col_speed_high_400,v_speed_high_400] = find(SICE_pca_speed_high_400);
% [row_speed_low_400,col_speed_low_400,v_speed_low_400] = find(SICE_pca_speed_low_400);
% [row_speed_high_200,col_speed_high_200,v_speed_high_200] = find(SICE_pca_speed_high_200);
% [row_speed_low_200,col_speed_low_200,v_speed_low_200] = find(SICE_pca_speed_low_200);
%
% [row_cost_high_700,col_cost_high_700,v_cost_high_700] = find(SICE_pca_cost_high_700);
% [row_cost_low_700,col_cost_low_700,v_cost_low_700] = find(SICE_pca_cost_low_700);
% [row_cost_high_400,col_cost_high_400,v_cost_high_400] = find(SICE_pca_cost_high_400);
% [row_cost_low_400,col_cost_low_400,v_cost_low_400] = find(SICE_pca_cost_low_400);
% [row_cost_high_200,col_cost_high_200,v_cost_high_200] = find(SICE_pca_cost_high_200);
% [row_cost_low_200,col_cost_low_200,v_cost_low_200] = find(SICE_pca_cost_low_200);

% %plot histograms
% figure;
% subplot(1,2,1);
% hist(v_innovation_high_700);
% title('Nonzero values distribution of high-flexibility IC outcome-700 ');
% subplot(1,2,2);
% hist(v_innovation_low_700);
% title('Nonzero values distribution of low-flexibility IC outcome-700');


% IC analysis for specific penalty parameters
% SICE_pca_innovation_high  = sparseInverseCovariance(cov_pca_innovation_high, 2.197, []);
% SICE_pca_innovation_low  = sparseInverseCovariance(cov_pca_innovation_low, 2.197, []);
%
% SICE_pca_quality_high  = sparseInverseCovariance(cov_pca_quality_high, 2.197, []);
% SICE_pca_quality_low  = sparseInverseCovariance(cov_pca_quality_low, 2.197, []);
%
% SICE_pca_flexibility_high  = sparseInverseCovariance(cov_pca_flexibility_high, 2.197, []);
% SICE_pca_flexibility_low  = sparseInverseCovariance(cov_pca_flexibility_low, 2.197, []);
%
% SICE_pca_speed_high  = sparseInverseCovariance(cov_pca_speed_high, 2.197, []);
% SICE_pca_speed_low  = sparseInverseCovariance(cov_pca_speed_low, 2.197, []);
%
% SICE_pca_cost_high  = sparseInverseCovariance(cov_pca_cost_high, 2.197, []);
% SICE_pca_cost_low  = sparseInverseCovariance(cov_pca_cost_low, 2.197, []);

%illustrate of PCA-based SICE outcome
% figure;
% subplot(1,2,1);
% spy(SICE_pca_innovation_high);
% title('PCA-based IC Matrix of High-Innovation Records');
% subplot(1,2,2);
% spy(SICE_pca_innovation_low);
% title('PCA-based IC Matrix of Low-Innovation Records');
%
% figure;
% subplot(1,2,1);
% spy(SICE_pca_quality_high);
% title('PCA-based IC Matrix of High-Quality Records');
% subplot(1,2,2);
% spy(SICE_pca_quality_low);
% title('PCA-based IC Matrix of Low-Quality Records');
%
% figure;
% subplot(1,2,1);
% spy(SICE_pca_flexibility_high);
% title('PCA-based IC Matrix of High-Flexibility Records');
% subplot(1,2,2);
% spy(SICE_pca_flexibility_low);
% title('PCA-based IC Matrix of Low-Flexibility Records');
%
% figure;
% subplot(1,2,1);
% spy(SICE_pca_speed_high);
% title('PCA-based IC Matrix of High-Speed Records');
% subplot(1,2,2);
% spy(SICE_pca_speed_low);
% title('PCA-based IC Matrix of Low-Speed Records');
%


% figure;
% subplot(1,2,1);
% spy(SICE_pca_cost_high_700);
% title('PCA-based IC Matrix of High-Cost-700 Records');
% subplot(1,2,2);
% spy(SICE_pca_cost_low_700);
% title('PCA-based IC Matrix of Low-Cost-700 Records');
%
% figure;
% subplot(1,2,1);
% spy(SICE_pca_cost_high_400);
% title('PCA-based IC Matrix of High-Cost-400 Records');
% subplot(1,2,2);
% spy(SICE_pca_cost_low_400);
% title('PCA-based IC Matrix of Low-Cost-400 Records');

figure;
subplot(1,2,1);
spy(SICE_pca_innovation_high_200);
title('PCA-based IC Matrix of High-Innovation-200 Records');
subplot(1,2,2);
spy(SICE_pca_innovation_low_200);
title('PCA-based IC Matrix of Low-Innovation-200 Records');

figure;
subplot(1,2,1);
spy(SICE_pca_quality_high_200);
title('PCA-based IC Matrix of High-Quality-200 Records');
subplot(1,2,2);
spy(SICE_pca_quality_low_200);
title('PCA-based IC Matrix of Low-Quality-200 Records');

figure;
subplot(1,2,1);
spy(SICE_pca_flexibility_high_200);
title('PCA-based IC Matrix of High-flexibility-200 Records');
subplot(1,2,2);
spy(SICE_pca_flexibility_low_200);
title('PCA-based IC Matrix of Low-flexibility-200 Records');

figure;
subplot(1,2,1);
spy(SICE_pca_speed_high_200);
title('PCA-based IC Matrix of High-Speed-200 Records');
subplot(1,2,2);
spy(SICE_pca_speed_low_200);
title('PCA-based IC Matrix of Low-Speed-200 Records');

figure;
subplot(1,2,1);
spy(SICE_pca_cost_high_200);
title('PCA-based IC Matrix of High-Cost-200 Records');
subplot(1,2,2);
spy(SICE_pca_cost_low_200);
title('PCA-based IC Matrix of Low-Cost-200 Records');

BGobj_pca_innovation_high_200 = biograph(SICE_pca_innovation_high_200);
BGobj_pca_innovation_low_200 = biograph(SICE_pca_innovation_low_200);
BGobj_pca_quality_high_200 = biograph(SICE_pca_quality_high_200);
BGobj_pca_quality_low_200 = biograph(SICE_pca_quality_low_200);
BGobj_pca_flexibility_high_200 = biograph(SICE_pca_flexibility_high_200);
BGobj_pca_flexibility_low_200 = biograph(SICE_pca_flexibility_low_200);
BGobj_pca_speed_high_200 = biograph(SICE_pca_speed_high_200);
BGobj_pca_speed_low_200 = biograph(SICE_pca_speed_low_200);
BGobj_pca_cost_high_200 = biograph(SICE_pca_cost_high_200);
BGobj_pca_cost_low_200 = biograph(SICE_pca_cost_low_200);

% to view network of SICE outcome-PCA-based,grouped
view(BGobj_pca_innovation_high_200);
view(BGobj_pca_innovation_low_200);
view(BGobj_pca_quality_high_200);
view(BGobj_pca_quality_low_200);
view(BGobj_pca_flexibility_high_200);
view(BGobj_pca_flexibility_low_200);
view(BGobj_pca_speed_high_200);
view(BGobj_pca_speed_low_200);

% to extract construct pairs' names of SICE outcome-nonzero
% k = 1;
% for i = 1: length(SICE_pca_cost_high_700)
%     for j =  (i + 1): length(SICE_pca_cost_high_700)
%         if(i ~= j)  %exclude values on diagonal
%            if(abs(SICE_pca_cost_high_700 (i,j)) ~=0)
%                arc_pca_cost_high_700(k).var1 = pcadata_name(i,:);
%                arc_pca_cost_high_700(k).var2 = pcadata_name(j,:);
%                arc_pca_cost_high_700(k).value = SICE_pca_cost_high_700 (i,j);
%                k = k + 1;
%         end
%     end
%   end
% end
%
% k = 1;
% for i = 1: length(SICE_pca_cost_low_700)
%     for j =  (i + 1): length(SICE_pca_cost_low_700)
%         if(i ~= j)  %exclude values on diagonal
%            if(abs(SICE_pca_cost_low_700 (i,j)) ~=0)
%                arc_pca_cost_low_700(k).var1 = pcadata_name(i,:);
%                arc_pca_cost_low_700(k).var2 = pcadata_name(j,:);
%                arc_pca_cost_low_700(k).value = SICE_pca_cost_low_700 (i,j);
%                k = k + 1;
%         end
%     end
%   end
% end
%
% k = 1;
% for i = 1: length(SICE_pca_cost_high_400)
%     for j =  (i + 1): length(SICE_pca_cost_high_400)
%         if(i ~= j)  %exclude values on diagonal
%            if(abs(SICE_pca_cost_high_400 (i,j)) ~=0)
%                arc_pca_cost_high_400(k).var1 = pcadata_name(i,:);
%                arc_pca_cost_high_400(k).var2 = pcadata_name(j,:);
%                arc_pca_cost_high_400(k).value = SICE_pca_cost_high_400 (i,j);
%                k = k + 1;
%         end
%     end
%   end
% end
%
% k = 1;
% for i = 1: length(SICE_pca_cost_low_400)
%     for j =  (i + 1): length(SICE_pca_cost_low_400)
%         if(i ~= j)  %exclude values on diagonal
%            if(abs(SICE_pca_cost_low_400 (i,j)) ~=0)
%                arc_pca_cost_low_400(k).var1 = pcadata_name(i,:);
%                arc_pca_cost_low_400(k).var2 = pcadata_name(j,:);
%                arc_pca_cost_low_400(k).value = SICE_pca_cost_low_400 (i,j);
%                k = k + 1;
%         end
%     end
%   end
% end

k = 1;
for i = 1: length(SICE_pca_cost_high_200)
    for j =  (i + 1): length(SICE_pca_cost_high_200)
        if(i ~= j)  %exclude values on diagonal
           if(abs(SICE_pca_cost_high_200 (i,j)) ~=0)
               arc_pca_cost_high_200(k).var1 = pcadata_name(i,:);
               arc_pca_cost_high_200(k).var2 = pcadata_name(j,:);
               arc_pca_cost_high_200(k).value = SICE_pca_cost_high_200 (i,j);
               k = k + 1;
        end
    end
  end
end

k = 1;
for i = 1: length(SICE_pca_cost_low_200)
    for j =  (i + 1): length(SICE_pca_cost_low_200)
        if(i ~= j)  %exclude values on diagonal
           if(abs(SICE_pca_cost_low_200 (i,j)) ~=0)
               arc_pca_cost_low_200(k).var1 = pcadata_name(i,:);
               arc_pca_cost_low_200(k).var2 = pcadata_name(j,:);
               arc_pca_cost_low_200(k).value = SICE_pca_cost_low_200 (i,j);
               k = k + 1;
        end
    end
  end
end

% to extract construct pairs' names of SICE outcome-Innovation-200
k = 1;
for i = 1: length(SICE_pca_innovation_high_200)
    for j =  (i + 1): length(SICE_pca_innovation_high_200)
        if(i ~= j)  %exclude values on diagonal
           if(abs(SICE_pca_innovation_high_200 (i,j)) ~=0)
               arc_pca_innovation_high_200(k).var1 = pcadata_name(i,:);
               arc_pca_innovation_high_200(k).var2 = pcadata_name(j,:);
               arc_pca_innovation_high_200(k).value = SICE_pca_innovation_high_200 (i,j);
               k = k + 1;
        end
    end
  end
end

k = 1;
for i = 1: length(SICE_pca_innovation_low_200)
    for j =  (i + 1): length(SICE_pca_innovation_low_200)
        if(i ~= j)  %exclude values on diagonal
           if(abs(SICE_pca_innovation_low_200 (i,j) )~=0)
               arc_pca_innovation_low_200(k).var1 = pcadata_name(i,:);
               arc_pca_innovation_low_200(k).var2 = pcadata_name(j,:);
               arc_pca_innovation_low_200(k).value = SICE_pca_innovation_low_200 (i,j);
               k = k + 1;
        end
    end
  end
end

% to extract construct pairs' names of SICE outcome-quality-200
k = 1;
for i = 1: length(SICE_pca_quality_high_200)
    for j =  (i + 1): length(SICE_pca_quality_high_200)
        if(i ~= j)  %exclude values on diagonal
           if(abs(SICE_pca_quality_high_200 (i,j)) ~=0)
               arc_pca_quality_high_200(k).var1 = pcadata_name(i,:);
               arc_pca_quality_high_200(k).var2 = pcadata_name(j,:);
               arc_pca_quality_high_200(k).value = SICE_pca_quality_high_200 (i,j);
               k = k + 1;
        end
    end
  end
end

k = 1;
for i = 1: length(SICE_pca_quality_low_200)
    for j =  (i + 1): length(SICE_pca_quality_low_200)
        if(i ~= j)  %exclude values on diagonal
           if(abs(SICE_pca_quality_low_200 (i,j) )~=0)
               arc_pca_quality_low_200(k).var1 = pcadata_name(i,:);
               arc_pca_quality_low_200(k).var2 = pcadata_name(j,:);
               arc_pca_quality_low_200(k).value = SICE_pca_quality_low_200 (i,j);
               k = k + 1;
        end
    end
  end
end

% to extract construct pairs' names of SICE outcome-flexibility-200
k = 1;
for i = 1: length(SICE_pca_flexibility_high_200)
    for j =  (i + 1): length(SICE_pca_flexibility_high_200)
        if(i ~= j)  %exclude values on diagonal
           if(abs(SICE_pca_flexibility_high_200 (i,j)) ~=0)
               arc_pca_flexibility_high_200(k).var1 = pcadata_name(i,:);
               arc_pca_flexibility_high_200(k).var2 = pcadata_name(j,:);
               arc_pca_flexibility_high_200(k).value = SICE_pca_flexibility_high_200 (i,j);
               k = k + 1;
        end
    end
  end
end

k = 1;
for i = 1: length(SICE_pca_flexibility_low_200)
    for j =  (i + 1): length(SICE_pca_flexibility_low_200)
        if(i ~= j)  %exclude values on diagonal
           if(abs(SICE_pca_flexibility_low_200 (i,j) )~=0)
               arc_pca_flexibility_low_200(k).var1 = pcadata_name(i,:);
               arc_pca_flexibility_low_200(k).var2 = pcadata_name(j,:);
               arc_pca_flexibility_low_200(k).value = SICE_pca_flexibility_low_200 (i,j);
               k = k + 1;
        end
    end
  end
end

% to extract construct pairs' names of SICE outcome-speed-200
k = 1;
for i = 1: length(SICE_pca_speed_high_200)
    for j =  (i + 1): length(SICE_pca_speed_high_200)
        if(i ~= j)  %exclude values on diagonal
           if(abs(SICE_pca_speed_high_200 (i,j)) ~=0)
               arc_pca_speed_high_200(k).var1 = pcadata_name(i,:);
               arc_pca_speed_high_200(k).var2 = pcadata_name(j,:);
               arc_pca_speed_high_200(k).value = SICE_pca_speed_high_200 (i,j);
               k = k + 1;
        end
    end
  end
end

k = 1;
for i = 1: length(SICE_pca_speed_low_200)
    for j =  (i + 1): length(SICE_pca_speed_low_200)
        if(i ~= j)  %exclude values on diagonal
           if(abs(SICE_pca_speed_low_200 (i,j) )~=0)
               arc_pca_speed_low_200(k).var1 = pcadata_name(i,:);
               arc_pca_speed_low_200(k).var2 = pcadata_name(j,:);
               arc_pca_speed_low_200(k).value = SICE_pca_speed_low_200 (i,j);
               k = k + 1;
        end
    end
  end
end




%  %to form vector of  lambda
% for i = 1: 100
%      lambda(i) = i * 0.01;
% end
% 
% figure;
% subplot(1,2,1);
% plot(lambda, bic_innovation_high);
% title('BIC of innovation-high ');
% subplot(1,2,2);
% plot(lambda, aic_innovation_high);;
% title('AIC of innovation-high');
% 
% figure;
% subplot(1,2,1);
% plot(lambda, bic_innovation_low);
% title('BIC of innovation-low ');
% subplot(1,2,2);
% plot(lambda, aic_innovation_low);;
% title('AIC of innovation-low');
% 
% figure;
% subplot(1,2,1);
% plot(lambda, bic_quality_high);
% title('BIC of quality-high ');
% subplot(1,2,2);
% plot(lambda, aic_quality_high);;
% title('AIC of quality-high');
% 
% figure;
% subplot(1,2,1);
% plot(lambda, bic_quality_low);
% title('BIC of quality-low ');
% subplot(1,2,2);
% plot(lambda, aic_quality_low);;
% title('AIC of quality-low');
% 
% figure;
% subplot(1,2,1);
% plot(lambda, bic_flexibility_high);
% title('BIC of flexibility-high ');
% subplot(1,2,2);
% plot(lambda, aic_flexibility_high);;
% title('AIC of flexibility-high');
% 
% figure;
% subplot(1,2,1);
% plot(lambda, bic_flexibility_low);
% title('BIC of flexibility-low ');
% subplot(1,2,2);
% plot(lambda, aic_flexibility_low);;
% title('AIC of flexibility-low');
% 
% figure;
% subplot(1,2,1);
% plot(lambda, bic_speed_high);
% title('BIC of speed-high ');
% subplot(1,2,2);
% plot(lambda, aic_speed_high);;
% title('AIC of speed-high');
% 
% figure;
% subplot(1,2,1);
% plot(lambda, bic_speed_low);
% title('BIC of speed-low ');
% subplot(1,2,2);
% plot(lambda, aic_speed_low);;
% title('AIC of speed-low');
% 
% figure;
% subplot(1,2,1);
% plot(lambda, bic_cost_high);
% title('BIC of cost-high ');
% subplot(1,2,2);
% plot(lambda, aic_cost_high);;
% title('AIC of cost-high');
% 
% figure;
% subplot(1,2,1);
% plot(lambda, bic_cost_low);
% title('BIC of cost-low ');
% subplot(1,2,2);
% plot(lambda, aic_cost_low);;
% title('AIC of cost-low');

% to get strongest connections in each IC matrix according to optimal/adjusted lambda
% k = 1;
% for i = 1: 125
%     for j =  (i + 1): 125
%         if(i ~= j)  %exclude values on diagonal
%            if(abs(SICE_pca_innovation_high{5} (i,j)) >0.05)
%                best_innovation_high(k).var1 = pcadata_name(i,:);
%                best_innovation_high(k).var2 = pcadata_name(j,:);
%                best_innovation_high(k).value = SICE_pca_innovation_high{5} (i,j);
%                k = k + 1;
%             end
%         end
%     end
% end
% 
% k = 1;
% for i = 1: 125
%     for j =  (i + 1): 125
%         if(i ~= j)  %exclude values on diagonal
%            if(abs(SICE_pca_innovation_low{6} (i,j)) >0.05)
%                best_innovation_low(k).var1 = pcadata_name(i,:);
%                best_innovation_low(k).var2 = pcadata_name(j,:);
%                best_innovation_low(k).value = SICE_pca_innovation_low{6} (i,j);
%                k = k + 1;
%             end
%         end
%     end
% end
% 
% k = 1;
% for i = 1: 125
%     for j =  (i + 1): 125
%         if(i ~= j)  %exclude values on diagonal
%            if(abs(SICE_pca_quality_high{5} (i,j)) >0.05)
%                best_quality_high(k).var1 = pcadata_name(i,:);
%                best_quality_high(k).var2 = pcadata_name(j,:);
%                best_quality_high(k).value = SICE_pca_quality_high{5} (i,j);
%                k = k + 1;
%             end
%         end
%     end
% end
% 
% k = 1;
% for i = 1: 125
%     for j =  (i + 1): 125
%         if(i ~= j)  %exclude values on diagonal
%            if(abs(SICE_pca_quality_low{7} (i,j)) >0.05)
%                best_quality_low(k).var1 = pcadata_name(i,:);
%                best_quality_low(k).var2 = pcadata_name(j,:);
%                best_quality_low(k).value = SICE_pca_quality_low{7} (i,j);
%                k = k + 1;
%             end
%         end
%     end
% end
% 
% k = 1;
% for i = 1: 125
%     for j =  (i + 1): 125
%         if(i ~= j)  %exclude values on diagonal
%            if(abs(SICE_pca_flexibility_high{5} (i,j)) >0.05)
%                best_flexibility_high(k).var1 = pcadata_name(i,:);
%                best_flexibility_high(k).var2 = pcadata_name(j,:);
%                best_flexibility_high(k).value = SICE_pca_flexibility_high{5} (i,j);
%                k = k + 1;
%             end
%         end
%     end
% end
% 
% k = 1;
% for i = 1: 125
%     for j =  (i + 1): 125
%         if(i ~= j)  %exclude values on diagonal
%            if(abs(SICE_pca_flexibility_low{7} (i,j)) >0.05)
%                best_flexibility_low(k).var1 = pcadata_name(i,:);
%                best_flexibility_low(k).var2 = pcadata_name(j,:);
%                best_flexibility_low(k).value = SICE_pca_flexibility_low{7} (i,j);
%                k = k + 1;
%             end
%         end
%     end
% end
% 
% k = 1;
% for i = 1: 125
%     for j =  (i + 1): 125
%         if(i ~= j)  %exclude values on diagonal
%            if(abs(SICE_pca_speed_high{5} (i,j)) >0.05)
%                best_speed_high(k).var1 = pcadata_name(i,:);
%                best_speed_high(k).var2 = pcadata_name(j,:);
%                best_speed_high(k).value = SICE_pca_speed_high{5} (i,j);
%                k = k + 1;
%             end
%         end
%     end
% end
% 
% k = 1;
% for i = 1: 125
%     for j =  (i + 1): 125
%         if(i ~= j)  %exclude values on diagonal
%            if(abs(SICE_pca_speed_low{7} (i,j)) >0.05)
%                best_speed_low(k).var1 = pcadata_name(i,:);
%                best_speed_low(k).var2 = pcadata_name(j,:);
%                best_speed_low(k).value = SICE_pca_speed_low{7} (i,j);
%                k = k + 1;
%             end
%         end
%     end
% end
% 
% k = 1;
% for i = 1: 125
%     for j =  (i + 1): 125
%         if(i ~= j)  %exclude values on diagonal
%            if(abs(SICE_pca_cost_high{6} (i,j)) >0.05)
%                best_cost_high(k).var1 = pcadata_name(i,:);
%                best_cost_high(k).var2 = pcadata_name(j,:);
%                best_cost_high(k).value = SICE_pca_cost_high{6} (i,j);
%                k = k + 1;
%             end
%         end
%     end
% end
% 
% k = 1;
% for i = 1: 125
%     for j =  (i + 1): 125
%         if(i ~= j)  %exclude values on diagonal
%            if(abs(SICE_pca_cost_low{6} (i,j)) >0.05)
%                best_cost_low(k).var1 = pcadata_name(i,:);
%                best_cost_low(k).var2 = pcadata_name(j,:);
%                best_cost_low(k).value = SICE_pca_cost_low{6} (i,j);
%                k = k + 1;
%             end
%         end
%     end
% end
% 
% plot the IC matrix
% figure;
% subplot(1,2,1);
% spy(SICE_pca_innovation_high{5});
% title('Optimal IC Matrix of High-innovation');
% subplot(1,2,2);
% spy(SICE_pca_innovation_low{6});
% title('Optimal IC Matrix of Low-innovation');
% 
% figure;
% subplot(1,2,1);
% spy(SICE_pca_quality_high{5});
% title('Optimal IC Matrix of High-quality');
% subplot(1,2,2);
% spy(SICE_pca_quality_low{7});
% title('Optimal IC Matrix of Low-quality');
% 
% figure;
% subplot(1,2,1);
% spy(SICE_pca_flexibility_high{5});
% title('Optimal IC Matrix of High-flexibility');
% subplot(1,2,2);
% spy(SICE_pca_flexibility_low{7});
% title('Optimal IC Matrix of Low-flexibility');
% 
% figure;
% subplot(1,2,1);
% spy(SICE_pca_speed_high{5});
% title('Optimal IC Matrix of High-speed');
% subplot(1,2,2);
% spy(SICE_pca_speed_low{7});
% title('Optimal IC Matrix of Low-speed');
% 
% figure;
% subplot(1,2,1);
% spy(SICE_pca_cost_high{6});
% title('Optimal IC Matrix of High-cost');
% subplot(1,2,2);
% spy(SICE_pca_cost_low{6});
% title('Optimal IC Matrix of Low-cost');

% not necessary to figure out the networks in BIC/AIC method 
% BGobj_best_innovation_high = biograph(SICE_pca_innovation_high{5});
% BGobj_best_innovation_low = biograph(SICE_pca_innovation_low{6});
% view(BGobj_best_innovation_high);
% view(BGobj_best_innovation_low);
% 
% BGobj_best_quality_high = biograph(SICE_pca_quality_high{5});
% BGobj_best_quality_low = biograph(SICE_pca_quality_low{7});
% view(BGobj_best_quality_high);
% view(BGobj_best_quality_low);
% 
% BGobj_best_flexibility_high = biograph(SICE_pca_flexibility_high{5});
% BGobj_best_flexibility_low = biograph(SICE_pca_flexibility_low{7});
% view(BGobj_best_flexibility_high);
% view(BGobj_best_flexibility_low);
% 
% BGobj_best_speed_high = biograph(SICE_pca_speed_high{5});
% BGobj_best_speed_low = biograph(SICE_pca_speed_low{7});
% view(BGobj_best_speed_high);
% view(BGobj_best_speed_low);
% 
% BGobj_best_cost_high = biograph(SICE_pca_cost_high{6});
% BGobj_best_cost_low = biograph(SICE_pca_cost_low{6});
% view(BGobj_best_cost_high);
% view(BGobj_best_cost_low);

%}
