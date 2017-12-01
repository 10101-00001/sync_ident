%% written and developed by Uwe Altmann
%% please cite: Altmann, U. (2013). Synchronisation nonverbalen Verhaltens. Wiesbaden: VS Springer. ISBN 978-3-531-19815-6


%% ***********************************************************************
% *** convert the output of "select_sync_intervals" from frame in time
function [sync_intervals_new] = ...
    sync_intervals_frame2time(sync_intervals, X_axis_time, Y_axis_time)

    % proof the input
    if nargin < 3,
        error('Missing neseccary input.');
    end
    
    % now, convert
    if isempty(sync_intervals),
        
        sync_intervals_new = [];  % it is no work to do
        
    else
        
        % copy
        sync_intervals_new = sync_intervals;
        
        % modify the time vectors 
        sync_intervals_new(:, 1) = Y_axis_time( sync_intervals(:, 1) );
        sync_intervals_new(:, 2) = X_axis_time( sync_intervals(:, 2) );
        sync_intervals_new(:, 3) = X_axis_time( sync_intervals(:, 3) );  
        
    end
    
%% ***********************************************************************