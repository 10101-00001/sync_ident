%% written and developed by Uwe Altmann
%% please cite: Altmann, U. (2013). Synchronisation nonverbalen Verhaltens. Wiesbaden: VS Springer. ISBN 978-3-531-19815-6

%% ********************************************************
%% ***  compute the windowed cross-lagged correlation (WCLC)
function [R2] = ...
            compute_WCLC(data, bandwidth, step, max_lag, add_noise)

    %If smooth=4 was selected: function moving_median.m required
    
    % the variablen must have 3 columns: frame number; patient_METS, therapist_METS
    
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
     
    % *** compute WCLC for person 1
    disp('WCLC for person 1 is computed.');
    [R2_1] = ...
            compute_WCLC_for_one(data, ...
                bandwidth, step, max_lag, add_noise);

     disp('WCLC for person 2 is computed.');
    % *** compute WCLC for person 2
    [R2_2] = ...
            compute_WCLC_for_one(data(:, [1 3 2]), ...
                bandwidth, step, max_lag, add_noise);
            
    % *** concat  R2 of person 1 with R2 of person 2
    R2 = [R2_1(:,end:-1:2) R2_2];
   % correl_coef = [correl_coef1(:,end:-1:2) correl_coef2];
   
    
  %  if nargout > 1,   % output only if requested
  %      correl_coef_out = correl_coef;
  %  end
 

    
    
%% ********************************************************
%% *** windowed cross-lagged correlation (WCLC) ***
function [R2] = ...
            compute_WCLC_for_one(data, bandwidth, step, max_lag, add_noise)
               
        
    % *** pre set 
    minimum_p_value = .001;  % for tests of H_0 :  r = 0
    
    % *** model selection
    % WCLC  
            var_Y = [2];
            var_X = [3];
    
    % ** in the case of artifical data, it should noise added
    % ** otherwise the regression could fail because the regressors 
    % ** have no variance
    if add_noise == true,  
        data(:, var_Y) = ...
            data(:, var_Y) + ...
            random('Normal', 0, .1, size(data(:, var_Y)));
        
        data(:, var_X) = ...
            data(:, var_X) + ...
            random('Normal', 0, .1, size(data(:, var_X)));
    end

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

            % get the variance of Y, if it is zero, then break
            std_y = std( data(cases_dependent,   var_Y) );
            std_x = std( data(cases_independent, var_X) );

            % if the variance of Y = 0, the computation of the correlation 
            % is not necessary. R2 could set on 0.
            if std_y > 0 && std_x > 0,

                [r, p] = corr(data(cases_dependent,var_Y), ...
                              data(cases_independent, var_X));
                
                correl_coef(pos_x, lag_n) = r;
                
                if p < minimum_p_value,
                    R2(pos_x, lag_n) = r * r;
                else
                    R2(pos_x, lag_n) = 0;
                end
                              
            else
                R2(pos_x, lag_n) = 0;
            end

            % settings for the next step
            lag_n = lag_n + 1;
            if  lag_n > length(all_time_lags) || ...
                all_positions(pos_x) + bandwidth-1 + all_time_lags( lag_n ) > n_time_points ,
                do_next_loop = false;
            end

        end

    end


%% ********************************************************