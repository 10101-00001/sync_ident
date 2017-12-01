%% written and developed by Uwe Altmann
%% please cite: Altmann, U. (2013). Synchronisation nonverbalen Verhaltens. Wiesbaden: VS Springer. ISBN 978-3-531-19815-6

%% ********************************************************
%% ********************************************************
function [R2] = ...
            compute_WCLR(data, bandwidth, step, max_lag, add_noise)

    %If smooth=4 was selected: function moving_median.m required
    
    % the variablen must have 6 columns: patient_body, therapist_body,
    % background roi1+2, patient_head, therapist_head

    % *** set default values, check input
    if nargin < 5,
        add_noise = true;
        if nargin < 4,
            max_lag = 125;
            if nargin < 3,
                step = 1;
                 if nargin < 2,
                    bandwidth = 75;
                                if nargin <1,
                                     error('More input is needed.')
                                end
                 end    
            end
        end
    end
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
    % *** compute WCLR for person 1
        disp('WCLR for person 1 is computed.');
    [R2_1] = ...
            compute_WCLR_for_one(data, ...
                bandwidth, step, max_lag, add_noise);
    
    % *** compute WCLR for person 2
    % *** we have to change the columns of the data matrix and compute WCLR
         disp('WCLR for person 2 is computed.');
    [R2_2] = ...
            compute_WCLR_for_one(data(:, [1 3 2]), ...
                bandwidth, step, max_lag, add_noise);
    
    % *** concat  R2 of person 1 with R2 of person 2
    R2 = [R2_1(:,end:-1:2) R2_2];
    
    
%% ********************************************************
%% *** windowed cross-lagged regresstion (WCLR) ***
function [R2] = ...
            compute_WCLR_for_one(data, bandwidth, step, max_lag, add_noise)
                 
    % *** model selection
            var_Y = [2];
            var_X_with_cross_correl = [2 3];
            var_X_without_cross_correl = [2];
    
    % ** in the case of artifical data, it should noise added
    % ** otherwise the regression could fail because the regressors 
    % ** have no variance
    if add_noise == true,
        data = data + random('Normal', 0, .1, size(data));
    end

    % *** the empirical value should be bigger
    % *** choose alpha and compute 1 - alpha,  e.g. 1 - 0.001 = 0.999
    r_square_test_critical_value = finv(0.999, 1, bandwidth - 3);  

    
    % *** compute the start points of the intervals which later under study
    n_time_points = length(data(:, 1));
    all_positions = (1:step:(n_time_points - bandwidth +1))';
    all_time_lags = 0:step:max_lag;

    % *** WCLR
    for pos_x = 1:length(all_positions),

        % get the indenpendent variable
        cases_independent = (all_positions(pos_x):(all_positions(pos_x) + bandwidth-1))';

        % **  start settings for the loop 
        do_next_loop = true;    lag_n = 1;  

        % compute the relationship between Y and X for every time lag
        while do_next_loop == true,

            %  pos_y = pos_x + time_lag;
            cases_dependent = cases_independent + all_time_lags( lag_n );

            % get the variance of Y
            var_y = var( data(cases_dependent, var_Y) );
            var_x_without = var( data(cases_independent, var_X_without_cross_correl) );
            var_x_with= var( data(cases_independent, var_X_with_cross_correl) );

            % if the variance of Y = 0, the computation of the regression 
            % is not necessary. R2 could set on 0.
            if var_y > 0,
                % model without cross-correlation
                [c, dev, stats] = glmfit( ...
                    data(cases_independent, var_X_without_cross_correl), ...
                    data(cases_dependent, var_Y) );
                R2_model_1 = 1 - var(stats.resid)/var_y;

                % model with cross-correlation
                [c, dev, stats] = glmfit( ...
                    data(cases_independent, var_X_with_cross_correl), ...
                    data(cases_dependent, var_Y) );
                R2_model_2 = 1 - var(stats.resid)/var_y;
                
                if R2_model_2 == 1, % this aviods problems with r square test
                    R2_model_2 = 0.9999999;
                end
                

                % *** r square difference test  and  save r square
                r2_diff = R2_model_2 - R2_model_1;

                if r2_diff * (bandwidth-3) / (1 - R2_model_2) < r_square_test_critical_value,
                    R2(pos_x, lag_n) = 0;
                else
                    R2(pos_x, lag_n) = r2_diff;
                end
            
            else
                R2(pos_x, lag_n) = 0;
                %disp(['Variance( Y ) < 0.1 at pos_x = ', ...
                %      num2str(pos_x),' and lag_n = ', num2str(lag_n) ]);
            end

            % settings for the next step: 
            % 1) last time lag reached?  2) end of time series reached?
            lag_n = lag_n + 1;
            if  lag_n > length(all_time_lags) || ...
                all_positions(pos_x) + bandwidth-1 + all_time_lags( lag_n ) > n_time_points ,
                do_next_loop = false;
            end

        end

    end


%% ********************************************************