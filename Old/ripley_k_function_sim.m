function [rK, K, KNMax, KNMin, Kall] = ...
    ripley_k_function_sim(xy_pos, par)


    
    % KFUNCTION_MAIN computes Ripley's K function for soelected ROI in the data
    % [xK, K, KNMax, KNMin, Kall] = ...
    %     kfunction_main (dataXY_all, xlim1, xlim2, ylim1, ylim2, Klim, nSteps,
    %     envelopes, Nsimul, filename)
    %
    % xK - xcoordinates of Kfunction K, KNMax, KNMin are the simulated
    % envelopes of the K function (selected as max and min of Kall).
    % dataXY_all - input data, #points-by-2 matrix, each row correspond to xy
    % coordinate of the datapoint
    % xlim1, xlim2, ylim1, ylim2 - selected ROI
    % Klim - limits of the K function ([Klim_low, Klim_high]) should be less
    % then half of the smaler size of the ROI
    % nSteps - number of steps for computation of K function
    % envelopes - if set to 1 computes simulated envelopes for the same number
    % of datapoints
    % Nsimul -  number of simulations to compute envelops
    % filename - name of the file where data are written (in current director)
    % If ommited no file is created.
    
    if isfield(par,'box')
        xmin = par.box(1,1);
        xmax = par.box(1,2);
        ymin = par.box(1,3);
        ymax = par.box(1,4);
        box = [xmin xmax ymin ymax];
        dataXY = ROIdata(xy_pos, xmin, xmax, ymin, ymax, 0);
    else
        xmin = min(xy_pos(:,1));
        xmax = max(xy_pos(:,1));
        ymin = min(xy_pos(:,2));
        ymax = max(xy_pos(:,2));
        box = [xmin xmax ymin ymax];
        dataXY = xy_pos;
    end
        [par.nPoints,~,~] = size(dataXY);
% 
%     if ~isfield(par,'Klim')
%         disp('yes')
%         r = sqrt(xy_pos(:,1).^2 + xy_pos(:,2).^2);
%         par.Klim = [min(r) max(r)];
%     end
    par.rStep  = (par.Klim(2)-par.Klim(1))/par.nSteps;
    rK = (0:par.rStep: round(par.rStep*par.nSteps))';

    fprintf('Computing Ripley''s K-function for selected ROI... \n');
    K = ripley_k_function(dataXY, rK, box, 0);


    if par.envelopes
        Kall = zeros(length(rK), par.nSim);
        wb = waitbar(0,'Computing simulation envelopes...');
        for ii=1:par.nSim
            waitbar(ii/par.nSim,wb);
            simNoise = generateNoise(par.nPoints,box);
            KN = ripley_k_function(simNoise,rK, box, 0);
            Kall(:,ii) = KN;
            if ii==1 %first round
                KNMax = KN;
                KNMin = KN;
            else
                KNMax = max(KNMax, KN);
                KNMin = min(KNMin, KN);
            end
        end
    close (wb);
    else
        KNMax = [];
        KNMin = [];
    end




end