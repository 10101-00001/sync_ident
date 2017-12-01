%% written and developed by Désirée Thielemann and Uwe Altmann


function [] = analyse_time_series(TSname, method, trans, smooth)
%%First choose method, transformation and smoothing degree
% method: 'WCLC' or 'WCLR'
% trans: 'rawdata' or 'log_trans' or 'anscombe_trans'
% smooth: 'rawdata' or 'slight_smooth' or 'high_smooth'

%%%% Master File for computing list of sync intervals out of preprocessed
%%%% time series
%%%% time series should have 25 fps
%CAVE: Time series should be preprocessed --> moving median and
%standardization should be applied

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Step - Setting fixed parameters (may be adjusted)

% set maximum time lag to 5 seconds (= 125 frames)
% set bandwidth to 5 seconds (= 125 frames)
% calculate WCLR/WCLC in steps of 1
% add_noise = true to avoid convergence problems with respect to the WCLR
max_lag = 125; 
step = 1; 
add_noise=true;
bandwidth=125;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Step - Loading time series

% time series should have two columns (movement time series of 2 persons)
[pathstr, name, ext]=fileparts(TSname);
data=dlmread([pathstr, '\', name, ext]);
% length is set 1000 frames
data=data(1:1000,:);
framenum=1:(length(data));
data=[framenum' ...
     data(:,1) ...
     data(:,2) ]; 

 cd(pathstr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Step - Define a path for results

pathstr_results= ['Results_', num2str(bandwidth), '_' ,method, '_', ...
                   trans, '_', smooth, '_', datestr(now, 'yyyy-mm-dd')];

% Check if output directory of MEA already exists
if exist(pathstr_results) == 7 & exist([pathstr_results, '_new']) ~=7
    cd(pathstr_results);
    % Check if output file in folder already exists
    if exist([name,'_R2.mat']) == 2
        cd ..;
        warning('Analysis has been done before! - Rename folder using dialog box to continue analysis')
        disp('In case of renaming: _new will be added to the folder name')
        % Construct a questdlg with three options
        choice = questdlg('Do you want to rename the folder?', ...
            'Output directory already exists');
        % Handle response
        switch choice
            case 'Yes'              
                % Construct a new filename.
                pathstr_results = sprintf([pathstr_results,'_new']);
                mkdir(pathstr_results);
            case 'No'
                error(['Output directory ', pathstr_results, ...
                    ' already exists.' char(10) 'Copy existing folder in another folder or delete it.']);
            case 'Cancel'
                error(['Output directory ', pathstr_results, ...
                    ' already exists.' char(10) 'Copy existing folder in another folder or delete it.']);
        end 
    else
        cd ..;
    end
else
   mkdir(pathstr_results) 
end 


if exist(pathstr_results) == 7 & exist([pathstr_results, '_new']) == 7
    cd([pathstr_results, '_new']);
    if exist([name,'_R2.mat']) ~= 2
        pathstr_results=[pathstr_results, '_new'];
        cd ..;
    end
    if exist([name,'_R2.mat']) == 2
        cd ..;
        error('Analysis has been twice done before - Copy existing folder in another folder or delete it.')
    end   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Step - Apply transformation (1=raw data, 2=log trans, 3=anscombe trans)

if strcmp(trans,'rawdata')
    trans_num=1;
elseif strcmp(trans,'log_trans')
    trans_num=2;
elseif strcmp(trans,'anscombe_trans')
    trans_num=3;
end

data=transformation_v02(data, trans_num);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Step - Appy smoothing 1=rawdata, 2= slight_smooth, 3= high_smooth

if strcmp(smooth,'rawdata')
    smooth_num=1;
elseif strcmp(smooth,'slight_smooth')
    smooth_num=2;
elseif strcmp(smooth,'high_smooth')
    smooth_num=3;
end

data=smoothing(data, smooth_num, framenum);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. Step - Compute WCLC or WCLR (Please select the following parameters: 
%%%%%%%%%%%% data, bandwidth, step, max_lag, add_noise)

    % *** method selection
	switch method,
        case 'WCLC',  % WCLC
                disp('compute R2 for WCLC');
                [R2] = ...
                    compute_WCLC(data, bandwidth, step, max_lag, add_noise);
            % save R2
                save([pathstr, '\', pathstr_results, '\', name, '_WCLC_R2.mat'], 'R2');            
        case 'WCLR',  % WCLR 
                disp('compute R2 for WCLR');
                [R2] = ...
                    compute_WCLR(data, bandwidth, step, max_lag, add_noise);
            % save R2
                    save([pathstr, '\', pathstr_results, '\', name,'_WCLR_R2.mat'], 'R2');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7. Step Define a mesh grid and save

	X_axis_time = data( (1:step:(length(data(:, 1)) - bandwidth +1))'...
                        + fix(bandwidth/2), 1 );
    Y_axis_time = data( step:step:max_lag, 1 );
    Y_axis_time = [ -1*Y_axis_time(end:-1:1) ; 0; Y_axis_time ]';
    [Xaa, Xab]  = ndgrid(X_axis_time, Y_axis_time);
    
    save([pathstr, '\',pathstr_results, '\', name, '_X_axis_time.mat'], 'X_axis_time');     
    save([pathstr, '\',pathstr_results, '\', name, '_Y_axis_time.mat'], 'Y_axis_time');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%8. Step Compute losis and save

disp(['Analyzing file ', name , ' losis.']);
[losi, ~]=find_all_sync_intervals_v06(...
            R2, X_axis_time, Y_axis_time,1, 10);
save([pathstr, '\', pathstr_results, '\', name, '_losi.mat'], 'losi');

end