% run_HawkePoisR(3,'x2020-03-14', 'x2020-03-16', 'confirm')
function [] = run_HawkPrBklyFeatRsth_sim_v2(delta, d_pred_start, d_pred_end, type_case, muReg, alphaScale, betaShape, mobiShift, iteration)
maxNumCompThreads(1);
%% Read in job parameter
%delta = 3;  d_pred_start = 'x2020-03-14';   d_pred_end   = 'x2020-03-16';   type_case = 'confirm';
delta = str2num(delta);

Method = 'HawkPrBklyFeatRsth_sim_v2';
Method_cells = strsplit( Method, '_' );
Method_cells = cell2mat(Method_cells(1));

mdl_load = ['./mdl/' Method_cells '/Mdl_' type_case  '_pd_start_' d_pred_start '_muReg_' muReg '_alphaScale_' alphaScale '_betaShape_' betaShape '_mobiShift_' mobiShift '.mat']

results_save = ['./results/' Method '/results_' type_case '_delta_' num2str(delta) '_pd_start_' d_pred_start '_pd_end_' d_pred_end ...
		'_muReg_' muReg '_alphaScale_' alphaScale '_betaShape_' betaShape '_mobiShift_' mobiShift '_itr_' iteration '.mat'];

disp(exist(results_save, 'file'))
disp(results_save)
%
if exist(results_save, 'file') == 2
	exit
else
end

%% Save out
muReg = str2num(muReg);
alphaScale_in = str2num(alphaScale);
betaShape_in  = str2num(betaShape);
mobiShift_in = str2num(mobiShift);

%% parameter
iteration_cells = strsplit( iteration, '-' );
sim_start = str2num( cell2mat(iteration_cells(1)) );
sim_end   = str2num( cell2mat(iteration_cells(2)) );

n_sim = sim_end - sim_start +1;
sim_max = 5;
rng(123);

%% Read model
load(mdl_load, 'mus','alpha','beta','K0','mdl','VarNames','alpha_delta','beta_delta','mus_delta','K0_delta','theta_delta')

%% Save out
% Read-in mobility

%%%%%%% Call out python to impute
%mobi_out = strrep(mdl_load,'mdl','imputation_temp_dir');
%mobi_out = strrep(mobi_out,'mat','csv');

if strcmpi(type_case, 'confirm')
        mobi_input = ['../InputSplit/Output/mobility_confirm_gt_10.csv'];
else
        mobi_input = ['../InputSplit/Output/mobility_death_gt_1.csv'];
end
% Get last day as pred_end day
last_day = readtable( '../InputSplit/Output/NYT_confirmed_gt_10_offset.csv','ReadVariableNames',true, 'Range', '1:2').Properties.VariableNames(end);
last_day = last_day{:};
last_day = strrep(last_day,'_','-');

mobi_out = ['./imputation_temp_dir/mobi_' type_case  '_pd_start_' d_pred_start '_last_day_' last_day '.csv']

%% Save out
% Read-in mobility
mob = readtable( mobi_out,'ReadVariableNames',true);
%delete(mobi_out)

if strcmpi(type_case, 'confirm')
    NYT = readtable(['../InputSplit/Output/NYT_confirmed_gt_10_offset.csv'],'ReadVariableNames',true);
    UScensus = readtable(['../InputSplit/Output/US_Census_confirm_gt_10.csv']);
    %
    Berkerly = readtable(['../InputSplit/Output/Berkerly_Feat_confirm_gt_10.csv']);
else
    NYT = readtable(['../InputSplit/Output/NYT_deaths_gt_1_offset.csv'],'ReadVariableNames',true);
    UScensus = readtable(['../InputSplit/Output/US_Census_death_gt_1.csv']);
    %
    Berkerly = readtable(['../InputSplit/Output/Berkerly_Feat_death_gt_1.csv']);
end

%%%%%%%%%Pad to shift %%%%%%%%%%%%%%%
mob_head = mob(:,1:4);
mob_value = table2array(mob(:,5:end));

for pad = 1:mobiShift_in
    mob_value = [ mean(mob_value(:,1:7),2) mob_value ];
end

mob_value = mob_value(:, 1:size(mob(:,5:end),2));
mob(:,5:end) = array2table(mob_value);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
Mobi_Type_list = table2cell(mob(1:6,4));
Mobi_Date_list = mob.Properties.VariableNames(5:end);
Mobi_Key_list = table2cell(mob(1:6:end,1:3));
%
NYT_Date_list = NYT.Properties.VariableNames(4:end);
NYT_Key_list = table2cell(NYT(:,1:3));
%
UScensus_Type_list = UScensus.Properties.VariableNames(2:end);
UScensus_Key_list  = table2cell(UScensus(:,1));
%
Berkerly_Type_list = Berkerly.Properties.VariableNames(2:end);
Berkerly_Key_list  = table2cell(Berkerly(:,1));
%
mob_val = table2array( mob(:,5:end));
%
NYT_val = table2array( NYT(:,4:end));
NYT_val(isnan(NYT_val)) = 0;
%
UScensus_val = table2array( UScensus(:,2:end));
%
Berkerly_val = table2array( Berkerly(:,2:end));
%
%% Data Pre-Processing
covid = NYT_val; 

% Calculate the difference from accumulative sum of cases
% The first day is padded with zero (i.e., assusme the first day has no cases)
covid = [zeros(size(covid,1), 1 ) covid(:,2:end) - covid(:,1:end-1)];
covid(covid<=0) = 0;

% Train & Test Split
% change the name
d_pred_start = strrep(d_pred_start,'-','_');
d_pred_end = strrep(d_pred_end,'-','_');
%
n_tr_end = find(strcmpi(NYT.Properties.VariableNames, {d_pred_start} )) - 4;
covid_tr = covid(:, 1:n_tr_end);
covid_te = covid(:, n_tr_end+1:n_tr_end+delta);
%
n_tr_end = find(strcmpi(mob.Properties.VariableNames, {d_pred_start} )) - 5;
mob_tr = mob_val(:, 1:n_tr_end);
mob_te = mob_val(:, n_tr_end+1:n_tr_end+delta);
%

% Get number of counties and number of days
[n_cty, n_day_tr]=size(covid_tr);
[n_cty, n_day_te]=size(covid_te);
%
% Normalization 
mob_tr_reshape = reshape(mob_tr, 6, size(mob_tr,1)/6 * size(mob_tr,2) ).';
mob_te_reshape = reshape(mob_te, 6, size(mob_te,1)/6 * size(mob_te,2) ).';
%
census_in = UScensus_val(:,end);
census_tr = repmat(census_in, n_day_tr,1);
census_te = repmat(census_in, n_day_te,1);
%
berkerly_in = Berkerly_val(:,1:end);
berkerly_tr = repmat(berkerly_in, n_day_tr,1);
berkerly_te = repmat(berkerly_in, n_day_te,1);
%
Covar_tr = [mob_tr_reshape berkerly_tr];
Covar_te = [mob_te_reshape berkerly_te];
%
Covar_tr_mean = mean(Covar_tr,1);
Covar_tr_std = std(Covar_tr,1);
%
Covar_tr = (Covar_tr-Covar_tr_mean) ./ Covar_tr_std;
Covar_te = (Covar_te-Covar_tr_mean) ./ Covar_tr_std;
%

% The total number of days 
%T=n_day;%a = 10; b = 4; % hard code shape parameters of Weibull estimates

% Get Variable names
VarNamesOld = [ Mobi_Type_list; Berkerly_Type_list.'; {['Qprob']}];
VarNames=[];
% Rename
for i = 1:size(VarNamesOld,1)
    newStr = replace( VarNamesOld{i} , ' & ' , '_' );
    newStr = replace( newStr , ' ' , '_' );
    VarNames=[VarNames; {newStr}];
end

% Remove Workplace
Idx_remove = find( strcmp(VarNames,'_Workplace') );


Covar_tr = [Covar_tr(:, [(1:Idx_remove-1) (Idx_remove+1:size(Covar_tr,2))] )];
Covar_te = [Covar_te(:, [(1:Idx_remove-1) (Idx_remove+1:size(Covar_te,2))])];
VarNames = [VarNames([(1:Idx_remove-1) (Idx_remove+1:size(VarNames,1))])];

%% Get K0
Covar_all = [Covar_tr; Covar_te];
n_day = n_day_tr+n_day_te;

%% Predict
[ypred,yci] = predict(mdl,Covar_all);
fK0 = reshape(ypred, n_cty, n_day);

% Make fK0 stable
% fK0(find(fK0>1.3))=1.3;

T_sim = n_day;

% Simulation results
sim = zeros(n_cty, T_sim, n_sim);
sim_ct_red = zeros(n_cty, n_sim);
sim_ct_exp = zeros(n_cty, n_sim);

% Loop for simulation
for c = 1:ceil(n_cty*0.5)  
    
    tr_in = covid_tr(c,:);
    for itr = sim_start:sim_end
        rng(itr);
        % disp(['itr: ' num2str(itr)])
        pop = UScensus_val(c,1);
        
        explode=1; sim_ct = 0;
        while(explode==1 && sim_ct<sim_max)
            [times_sim, explode] = Hawkes_Sim_Corona_MEMv2( mus(c), alpha, beta, T_sim, fK0(c,:), T_sim-delta, explode, itr, pop,  tr_in);
            sim_ct = sim_ct +1;
        end
        %
        [Nt]=discrete_hawkes(times_sim,T_sim);
        %
        sim(c,:,itr-sim_start+1) = Nt;
        sim_ct_red(c,itr-sim_start+1 ) = sim_ct;
        sim_ct_exp(c,itr-sim_start+1) = (sim_ct==sim_max);
        
        disp(['itr: ' num2str(itr) ' cty: ' num2str(c) ' fK0max: ' num2str(max(fK0(c,:))) ' pop: ' num2str(pop) ' sim_ct: ' num2str(sim_ct) ])
    end
    

end

sim_out = sim(:,end-delta+1:end,:);
mean_sim = mean(sum( sim_out, 2),3);
gt = sum(covid_te,2);
save(results_save, 'sim_out', 'mean_sim', 'gt', 'sim_ct_red', 'sim_ct_exp');

end
