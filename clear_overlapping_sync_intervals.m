%% written and developed by Uwe Altmann
%% please cite: Altmann, U. (2013). Synchronisation nonverbalen Verhaltens. Wiesbaden: VS Springer. ISBN 978-3-531-19815-6

%% ***********************************************************************
% ***  clear overlapping sync intervals
% ***  input : list of sync intervals (it is the output of connect_peaks)
% ***  output: list of sync intervals
function [losi_new] = clear_overlapping_sync_intervals( losi )

    % pre set
    losi_new = [];
    
    % now, look for overlapping intervals
    if ~isempty( losi ),
       
        % at first, add a column with the length  (it is the 5. coloumn)
        losi = [losi [losi(:,3)-losi(:,2)+1 ]];
        
        % add a column with distance to the zero line (means time lag = 0)
        % (it is the 6. coloumn)
        losi = [losi abs(losi(:,1))];

        % adjustment of the beginning with the time lag
        %losi(:,2) = losi(:,2) + losi(:,6);
        
        % sort by R2, distance to zero line, the length, and beginning
        % the best candidate is in first line, the badest in the last line
        losi = sortrows( losi, [-4 -6 -5 2] );
        
        % *** find overlapping intervals and delete the interval with the
        % *** smaller R2
        while ~isempty( losi ),
        
            % *** save the first line, because it is the best candidate
            losi_new = [losi_new ; losi(1,:) ];
            losi(1,:) = [];
            
            % *** compare with the bader candidates 
            candidates_to_delete = ...
                (losi(:,2) >= losi_new(end,2) & ...   % overlapping begin
                    losi(:,2) <= losi_new(end,3)) | ... 
                (losi(:,3) >= losi_new(end,2) & ...   % overlapping end
                    losi(:,3) <= losi_new(end,3)) ; %| ... 
                
            % *** delete
            losi(candidates_to_delete,:) = [];
                
        end
        
        % *** find overlying intervals and delete the interval with the
        % *** smaller R2
        losi = losi_new;
        losi_new = [];
        
        while ~isempty( losi ),
        
            % *** save the first line, because it is the best candidate
            losi_new = [losi_new ; losi(1,:) ];
            losi(1,:) = [];
            
            % *** compare with the bader candidates 
            candidates_to_delete = ...
                (losi(:,2) >= losi_new(end,2) & ...   % old is within new
                    losi(:,3) <= losi_new(end,3)) | ... 
                (losi(:,2) <= losi_new(end,2) & ...   % new is within old
                    losi(:,3) >= losi_new(end,3)) ; 
                
            % *** delete
            losi(candidates_to_delete,:) = [];
                
        end
        
        
        % delete the addional columns
        losi_new(:,5:6) = [];
        
    end

%% ***********************************************************************