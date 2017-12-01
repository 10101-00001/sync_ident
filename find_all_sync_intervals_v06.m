%% written and developed by Uwe Altmann
%% please cite: Altmann, U. (2013). Synchronisation nonverbalen Verhaltens. Wiesbaden: VS Springer. ISBN 978-3-531-19815-6


% ***********************************************************************
% ***  find all sync intervals, using the R^2 matrix

% ***  input  : the R2 matrix (output of compute_WCLC or compute_WCLR)
% ***           X_axis_time - output of axis_frame2time
% ***           Y_axis_time - output of axis_frame2time
% ***           time_lag_tolerance -> integer how much changes of time lag
%               allowed
% ***           minimum_length -> minimum length of an sync interval

% ***  output : losi  = list of sync interval, 1col time lag, 2col begin,
%                       3col end, 4col mean(R^2) which sync interval
% ***           loR2p = list of R2 peaks, for the plots
% ***           lop   = list of peaks, peak at   1col = time lag,  2col = time point

function [losi, loR2p, lop] = ...
    find_all_sync_intervals_v06(R2, X_axis_time, Y_axis_time, ...
             time_lag_tolerance, minimum_length)

    % proof the input
    if nargin < 5,
        minimum_length = 5;          % set a default value
        if nargin < 4,
            time_lag_tolerance = 1;  % set a default value
            if nargin < 3,
                error('I need R2, X_axis_time, and Y_axis_time as input.');
            end
        end
    end
    
    if minimum_length < 0 || minimum_length ~= round(minimum_length) || ...
            time_lag_tolerance < 0 || ...
            time_lag_tolerance ~= round(time_lag_tolerance),
        error(['time_lag_tolerance and minimum_length must be ', ...
               'positive integer.'])
    end
    
    % proof X_axis_time 
    v_s = size(X_axis_time);
    
    if v_s(1) < 2 || v_s(2) ~= 1,
        error(['X_axis_time must be column-vector. '...
               'Currently X_axis_time has ', num2str(v_s(1)), ...
               'columns and ', num2str(v_s(2)), ' rows'])
    end
    
    % proof Y_axis_time
    v_s = size(Y_axis_time);
    
    if v_s(1) ~= 1 || v_s(2) < 2,
        error(['Y_axis_time must be row-vector. '...
               'Currently Y_axis_time has ', num2str(v_s(1)), ...
               ' columns and ', num2str(v_s(2)), ' rows'])
    end
    
    
    % *********************************************************************
    % peak picking of R2  
    % list of peaks 1col = time lag, 2col = time point
    disp('Peak Picking starts.');
    [lop] = peak_picking_of_R2( R2 );
    clear R2;
    
    % connect peaks given a time lag tolerance and a sync minimum length
    % loR2p: list of R2 peaks 
    disp('Connecting peaks...');
        [loR2p] = connect_peaks( lop, time_lag_tolerance, minimum_length );
    
    % convert from frame in time
    disp('Convert from frame to time');
    [losi] = sync_intervals_frame2time( loR2p, X_axis_time, Y_axis_time);
    
    % overlapping sync intervals will be deleted
    disp('Delete overlapping Sync-intervals.');
    [losi] = clear_overlapping_sync_intervals( losi );
    
end % function
    
    
    
% ***********************************************************************
% ***  connect peaks given a time lag tolerance and a sync minimum length
% ***  input  : list_of_peaks is the output of peak_picking_of_R2
% ***  output : losi = list of sync intervals
% ***  losi: 1col = time lag,  2col = begin,  3col = end,  4col = mean(R2)
function [losi] = connect_peaks( ...
                     list_of_peaks, time_lag_tolerance, minimum_length)

    % proof the input ************************************
    if nargin < 3,
        minimum_length = 5;   % set default
        disp('For the connection of peak intervals the default value is used.')
        disp(['minimum_length = ', num2str(minimum_length)])
        if nargin < 2,
            time_lag_tolerance = 1;  % set default
            disp('For the connection of peak intervals the default value is used.')
            disp(['time_lag_tolerance = ', num2str(time_lag_tolerance)])
            if nargin < 1,
                error('list_of_peaks as input is necessary.');
            end
        end
    end

    if minimum_length < 1 || minimum_length ~= round(minimum_length),
        error('minimum_length must be an interger and larger than 0.');
    end
    
    if time_lag_tolerance < 0 || ...
             time_lag_tolerance ~= round(time_lag_tolerance),
        error(['time_lag_tolerance must be an interger ', ...
               'and larger than 0 or equal 0.']);
    end   

    % pre set *********************************************
    losi = [];
    
    % now, try to find the neighbour of a peak
    if ~isempty(list_of_peaks),
        
        % how many peaks?
        peaks_count = length( list_of_peaks(:,1));
        
        % at first, add a column with an case index
        list_of_peaks = [list_of_peaks (1:peaks_count)'];
   
        % now, find for every single peak the possible neighbours
         M = false( peaks_count, peaks_count);
    
        
        
        for ni = 1:peaks_count,
            
            % pre set
            list_of_neighbours = [];
            
            % give the reference peak the first place in the list
            lop = [list_of_peaks(ni,:); ...
                   list_of_peaks( list_of_peaks(:,end) ~= ni,:)];
            
            % find neighbours looking foreward and look backward
            list_of_neighbours = [list_of_neighbours ; ...
                find_neighbours_moving_foreward( lop, time_lag_tolerance) ; ...
                find_neighbours_moving_backward( lop, time_lag_tolerance) ];
            
            % delete the double element
            list_of_neighbours(end, :) = []; 

            % save the result
             M( list_of_neighbours(:,end), ni) = true;
             R2_mean(1, ni) = mean( list_of_neighbours(:,3) );
             
        end
        
        % clear unused variables
        clear list_of_neighbours i j lop
        
        % filter too short sync intervals -> minimum_length
        indexM = sum( M ) < minimum_length;
        M(:, indexM ) = [];
        R2_mean(:, indexM ) = [];
        
        % sort the intervals regarding to their length and (subsequent) R2 mean
        [~, indexM] =  sortrows( [sum( M ); R2_mean]' , [-1 -2] );
        M = M(:, indexM );
        R2_mean = R2_mean(:, indexM ); 
        clear indexM
        
        % delete sync intervals with double used elements
        new_M = [];
        new_R2_mean = [];
        while ~isempty(M),
            
            % save the largest sequence
            new_M = [new_M M(:,1)];
            M(:,1) = [];
            
            new_R2_mean = [new_R2_mean R2_mean(:,1)];
            R2_mean(:,1) = [];
            
            % find shorter sequences with the same elements
            index_row = new_M(:, end) == true;
            index_col = sum(M(index_row, :) == 1) > 0;
            M(:,index_col)=[];
            R2_mean(:,index_col) = [];

        end
        
        % save results as list of sync intervals (losi)
        % 1col = time lag,  2col = begin,  3col = end,  4col = mean(R2)
        losi = zeros( length(new_R2_mean), 4 );
        
        for ni = 1:length(new_R2_mean),
            
            losi(ni,1) = round( mean( list_of_peaks(logical(new_M(:,ni)), 2) )); % mean time lag
            losi(ni,2) = min( list_of_peaks(logical(new_M(:,ni)), 1) ); % begin of sequence
            losi(ni,3) = max( list_of_peaks(logical(new_M(:,ni)), 1) ); % end of sequence
            losi(ni,4) = new_R2_mean(1,ni);  % mean( R2)
            
        end
      
    end

    
end


% ***********************************************************************
% ***  peak picking of R2
% ***  input  : R2 matrix, it is the output of compute_WCLC or compute_WCLR
% ***  output : lop = list of peaks
% ***  lop    : 1col = time point, 2col = time lag, 3col = peak height
function [lop] = peak_picking_of_R2( R2 )

    % proof input parameters
    if nargin < 1,
        error('R2 as input is necessary.');
    end
    
    % set variabels for work
    lop = [];
    
    % now, look for peaks for every t separately
    for t = 1:length( R2(:,1) ),
       
       pks = []; locs = [];
       [pks, locs] = findpeaks( R2( t, :));
       
       if ~isempty(pks),  % peaks found
           for k = 1:length(pks),
               lop = [ lop ; ...
                           [t locs(k) pks(k)] ];

           end
       end
    end

    
end


% ***********************************************************************
% ***  for the first run, use only the variables 
% ***          list_of_peaks, time_lag_tolerance, minimum_length  !!!
% ***  
function [list_of_neighbours, list_of_peaks_out] = ...
    find_neighbours_moving_foreward( ...
        list_of_peaks, time_lag_tolerance, list_of_neighbours, lag)

    % ***  proof input
    if nargin < 4,
        lag = 0;
        if nargin < 3,
            list_of_neighbours = [];
            if nargin < 2,
                error('Input is not complete.');
            end
        end
    end


    % ***  now, find neighbours 
    if ~isempty(list_of_peaks),
        
        % ***  first run
        if isempty(list_of_neighbours),
            
            % set the first element as orgin
            list_of_neighbours = list_of_peaks(1,:);
            
            % begin the next step
            [list_of_neighbours, list_of_peaks] = ...
                        find_neighbours_moving_foreward( list_of_peaks, ...
                        time_lag_tolerance, list_of_neighbours, lag);

        else % find the neighbours of the late entry in the list_of_neighbours
        
            % clear all peaks which are t or before t, which could not a neighbour
            list_of_peaks( ...
                    list_of_peaks(:,1) <= list_of_neighbours(end,1), :) = [];

            % which time lag changings allowed?
            [allowed_lags] = which_lags_are_allowed(lag, time_lag_tolerance);
                
            % find possible neighbours
            neighbour_index = [];
            for i = 1:length(allowed_lags),
                
                neighbour_index = [neighbour_index ; ...
                    find( ...
                        list_of_neighbours(end,1)+1 == list_of_peaks(:,1) & ...  
                        list_of_neighbours(end,2)+allowed_lags(i) == list_of_peaks(:,2) ) ];
            end

            if ~isempty(neighbour_index),  % if neighbours found ...

                % if neighbours, choose the one with the largest R2
                candidates = list_of_peaks(neighbour_index,:);
                candidates = candidates( ...
                    max(candidates( :, 3)) == candidates( :, 3),:);
                    
                % the impossible case of 2 maxima
                if length( candidates( :, 3) ) > 1,
                    candidates = candidates(1,:);
                end    


                % increase resp. decrease the "lag"
                lag = lag + candidates(1,2) - list_of_neighbours(end,2) ;

                % save new neighbour
                list_of_neighbours = [list_of_neighbours ; candidates ];
                
                % delete the neighbour in the peak list
                list_of_peaks( ...
                    list_of_peaks(:,1) == candidates(1,1) & ...
                    list_of_peaks(:,2) == candidates(1,2) & ...
                    list_of_peaks(:,3) == candidates(1,3) , : ) = [];

                % now, look for the next neighbour
                [list_of_neighbours, list_of_peaks] = ...
                    find_neighbours_moving_foreward( list_of_peaks, ...
                    time_lag_tolerance, list_of_neighbours, lag);
            end
        end
    end
    
    % prepare output
    if nargout == 2,
        list_of_peaks_out = list_of_peaks;
    end
    
end


% ***********************************************************************
% ***  for the first run, use only the variables 
% ***          list_of_peaks, time_lag_tolerance, minimum_length  !!!
% ***  
function [list_of_neighbours, list_of_peaks_out] = ...
    find_neighbours_moving_backward( ...
        list_of_peaks, time_lag_tolerance, list_of_neighbours, lag)

    % ***  proof input
    if nargin < 4,
        lag = 0;
        if nargin < 3,
            list_of_neighbours = [];
            if nargin < 2,
                error('Input is not complete.');
            end
        end
    end


    % ***  now, find neighbours 
    if ~isempty(list_of_peaks),
        
        % ***  first run
        if isempty(list_of_neighbours),
            
            % set the first element as orgin
            list_of_neighbours = list_of_peaks(1,:);
            
            % begin the next step
            [list_of_neighbours, list_of_peaks] = ...
                        find_neighbours_moving_backward( list_of_peaks, ...
                        time_lag_tolerance, list_of_neighbours, lag);

        else % find the neighbours of the late entry in the list_of_neighbours
        
            % clear all peaks which are t or later, which could not a neighbour
            list_of_peaks( ...
                    list_of_peaks(:,1) >= list_of_neighbours(end,1), :) = [];

            % which time lag changings allowed?
            [allowed_lags] = which_lags_are_allowed(lag, time_lag_tolerance);
                
            % find possible neighbours regarding to the time lag
            neighbour_index = [];
            for i = 1:length(allowed_lags),
                
                neighbour_index = [neighbour_index ; ...
                    find( ...
                        list_of_neighbours(end,1)-1 == list_of_peaks(:,1) & ...  
                        list_of_neighbours(end,2)+allowed_lags(i) == list_of_peaks(:,2) ) ];
            end
            
            
            % if neighbours found ...
            if ~isempty(neighbour_index),  

                % if neighbours, choose the one with the largest R2
                candidates = list_of_peaks(neighbour_index,:);
                candidates = candidates( ...
                    max(candidates( :, 3)) == candidates( :, 3),:);
                    
                % the impossible case of 2 maxima
                if length( candidates( :, 3) ) > 1,
                    candidates = candidates(1,:);
                end    


                % increase resp. decrease the "lag"
                lag = lag + candidates(1,2) - list_of_neighbours(end,2) ;

                % save new neighbour
                list_of_neighbours = [list_of_neighbours ; candidates ];
                
                % delete the neighbour in the peak list
                list_of_peaks( ...
                    list_of_peaks(:,1) == candidates(1,1) & ...
                    list_of_peaks(:,2) == candidates(1,2) & ...
                    list_of_peaks(:,3) == candidates(1,3) , : ) = [];

                % now, look for the next neighbour (iterative algorithm)
                [list_of_neighbours, list_of_peaks] = ...
                    find_neighbours_moving_backward( list_of_peaks, ...
                    time_lag_tolerance, list_of_neighbours, lag);
            end
        end
    end
    
    % prepare output
    if nargout == 2,
        list_of_peaks_out = list_of_peaks;
    end

    
end 

% ***********************************************************************
% *** it is a sub-routine for find_neighbours_moving_foreward and
% *** find_neighbours_moving_backward
function [allowed_lags] = which_lags_are_allowed(lag, tolerance)

    % proof input
    if nargin < 2,
        error('Input is not complete.');
    end
    
    % *** pre set
    if tolerance == 0,
        allowed_lags = 0;
    else
        allowed_lags = [-1 0 1]';
        
        % *** lag should not be overrun the tolerance range
        if lag + 1 > abs(tolerance),  
            allowed_lags = [-1 0]';
        elseif lag - 1 < -1*abs(tolerance),
            allowed_lags = [ 0 1]';
        end
        
    end
    
    
end

% ***********************************************************************