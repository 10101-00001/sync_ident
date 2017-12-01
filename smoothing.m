%% written and developed by Uwe Altmann
%% please cite: Altmann, U. (2013). Synchronisation nonverbalen Verhaltens. Wiesbaden: VS Springer. ISBN 978-3-531-19815-6



function [data] = smoothing(data, which_smoothing, framenum) 

% 1 means no smoothing
% 2 means slight smoothing with cubic smoothing spline with parameter 0.9
% 3 means high smoothing with cubic smoothing spline with parameter 0.0005

% Check input
    if nargin<2,
    disp('Please select a smoothing procedure (1=no smoothing, 2=slight smoothing, 3=high smoothing)!')
        if nargin<1
            error('No input given! Please select a time series.')
        end
    end
    
          switch which_smoothing,
              case 1, % no smooth
                data=data;
              case 2,  % light smooth
                    me_prep=data(:,2:3);
                    for roi=1:2 
                        y_spline(:,roi) = csaps(framenum', me_prep(:,roi), 0.9);
                     end
                        y_new_1   = fnval( y_spline(:,1), framenum');
                        y_new_2   = fnval( y_spline(:,2), framenum');                    
                    data=[framenum' y_new_1 y_new_2];
               case 3,  % light smooth
                    me_prep=data(:,2:3);
                    for roi=1:2 
                        y_spline(:,roi) = csaps(framenum', me_prep(:,roi), 1-0.9995);
                    end
                        y_new_1   = fnval( y_spline(:,1), framenum');
                        y_new_2   = fnval( y_spline(:,2), framenum');                    
                data=[framenum' y_new_1 y_new_2];
               case 4,  % light smooth
                    me_prep=data(:,2:3);
                        y_med1   = moving_median( me_prep(:,1),3);
                        y_med2   = moving_median( me_prep(:,2),3);                    
                data=[framenum' y_med1 y_med2];
          end   
    

end