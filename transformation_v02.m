%% written and developed by Uwe Altmann
%% please cite: Altmann, U. (2013). Synchronisation nonverbalen Verhaltens. Wiesbaden: VS Springer. ISBN 978-3-531-19815-6


%% ********************************************************
%% *** compute transformation

function [data] = transformation(data, which_transformation)

    % data: matix with 7 columns, line 1 is frame 1, line x is frame x,
    % first colum is a time vector, next columns are a motion energy 
    % time series (Pat all, The all, Background L, Background R, Pat head, The head

    % which_transformation:
    % 1 means no transformation
    % 2 means box-cox transformation with \lambda = 0 equal to log trans
    % 3 means anscombe transformation

    % Check input
    if nargin<1
        error('No input given! Please select a time series.')
    end
    
    if nargin<2,
        disp('Please select a transformation (1=no trans, 2=log trans, 3=anscombe trans)! Now, no transformation is applied, procedure is being continued.')
        which_transformation = 1;
    end

    if ~(which_transformation==1 | which_transformation==2 | which_transformation==3),
        error('The parameter which_transformation must have one of these values: 1=no trans, 2=log trans, 3=anscombe trans.')
    end
    
    
    % **************************************
    switch which_transformation,
        case 1, % no transfromation
            data = data;     
            
        case 2, %log
            data(:,2:end) = log( data(:,2:end) + 1 );
            
        case 3, % Anscombe transformation
            data(:,2:end) = 2*sqrt(data(:,2:end) + (3/8));  
            
    end
    
       
end
