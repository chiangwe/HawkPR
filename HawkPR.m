function [] = HawkPR( InputPath_report, InputPath_mobility, InputPath_cencus, Delta, Alpha, Beta, OutputPath_mdl, OutputPath_pred)


%% Read in parameter

if strcmp(Alpha,'') && strcmp(Beta,'')
	disp('No shape and scale parameter for Weibull distribution provided. Use MLE to infer alpha and beta ... ')
else
	alphaScale_in = str2num(Alpha);
	betaShape_in  = str2num(Beta);
end 

if strcmp(Delta,'') 
        disp('No shift parameter for mobility provided.  It will set to zero ... ')
	mobiShift_in = 0;
else
	mobiShift_in = str2num(Delta);
end

% Read-in COVID data

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
% d_pred_end = strrep(d_pred_end,'-','_');
%
n_tr_end = find(strcmpi(NYT.Properties.VariableNames, {d_pred_start} )) - 4;
covid_tr = covid(:, 1:n_tr_end);
% covid_te = covid(:, n_tr_end+1:n_tr_end+delta);
%
n_tr_end = find(strcmpi(mob.Properties.VariableNames, {d_pred_start} )) - 5;
mob_tr = mob_val(:, 1:n_tr_end);
% mob_te = mob_val(:, n_tr_end+1:n_tr_end+delta);
%

% Get number of counties and number of days
[n_cty, n_day_tr]=size(covid_tr);
% [n_cty, n_day_te]=size(covid_te);
%
% Normalization
mob_tr_reshape = reshape(mob_tr, 6, size(mob_tr,1)/6 * size(mob_tr,2) ).';
% mob_te_reshape = reshape(mob_te, 6, size(mob_te,1)/6 * size(mob_te,2) ).';
%
census_in = UScensus_val(:,end);
census_tr = repmat(census_in, n_day_tr,1);
% census_te = repmat(census_in, n_day_te,1);
%
berkerly_in = Berkerly_val(:,1:end);
berkerly_tr = repmat(berkerly_in, n_day_tr,1);
%berkerly_te = repmat(berkerly_in, n_day_te,1);
%
Covar_tr = [mob_tr_reshape berkerly_tr];
%Covar_te = [mob_te_reshape berkerly_te];
%
Covar_tr_mean = mean(Covar_tr,1);
Covar_tr_std = std(Covar_tr,1);
%
Covar_tr = (Covar_tr-Covar_tr_mean) ./ Covar_tr_std;
%Covar_te = (Covar_te-Covar_tr_mean) ./ Covar_tr_std;
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
%Covar_te = [Covar_te(:, [(1:Idx_remove-1) (Idx_remove+1:size(Covar_te,2))])];
VarNames = [VarNames([(1:Idx_remove-1) (Idx_remove+1:size(VarNames,1))])];




%% Define Parameters
T = n_day_tr;
% Boundary correction, the number of days before the total number of days (n_day)
dry_correct = 5;

% EM step iterations
emiter = 200; %Nt=covid(:,2:end);
break_diff = 10^-4;
% Boundary correction: T-dry_correct
% Mobility has only 6 weeks so we take the less by min(T-dry_correct, size(mobi_in,2) )
day_for_tr = min(T-dry_correct, size(mob_tr,2) );

%% Initialize Inferred Parameters

if (alphaScale_in==0) && (betaShape_in==0)
    % Weibull Distribution (scale, shape) parameters as (alphas, betas)
    alpha  = 2;
    beta = 2;
else
    alpha  = alphaScale_in;
    beta = betaShape_in;
end
% K0 reproduction number, a fuction of time and mobility.
% Estimat for each county at each day.
K0 = 1*ones(n_cty, n_day_tr);

% p is the n_day by n_day matrix.
% p( i, j), i > j stands the probability of ONE SIGLE event at day i triggered by ALL events at day j
% i.e., Prob for (j_1, j_2, ...) triggered by each i
p=[];
for i = 1:n_cty
    p{i} = zeros(n_day_tr,n_day_tr);
end

% q is the n_day by n_day matrix.
% q( i, j), i > j stands the probability of ONE SIGLE event at day i triggered by ONE SIGLE event at day j
% i.e., Prob for each j triggered by each i
q=[];
for i = 1:n_cty
    q{i} = zeros(n_day_tr,n_day_tr);
end

% Mu is the back ground rate
mus=0.5*ones(n_cty,1);

% lam is the event intensity
lam = zeros(n_cty,T);

%% EM interation
alpha_delta = []; alpha_prev = [];
beta_delta = [];  beta_prev = [];
mus_delta = [];   mus_prev = [];
K0_delta = [];    K0_prev = [];
theta_delta = []; theta_prev = [];
for itr = 1:emiter
    
    %% E-step
    % county levelitr
    for c = 1:n_cty
        if( sum(covid_tr(c,:)) ~= 0)
            [p{c}, q{c}, lam(c,:)] = updatep( covid_tr(c,:) , p{c}, q{c}, K0(c,:), alpha, beta, mus(c) );
        end
    end
    
    %% M-step
    % Calculate Q, which stands for the average number (observed) of children generated by a SINGLE event j
    % Note taht the last "dry_correct" days of Q will be accurate
    % Since we haven't observed their children yet
    Q = [];
    for c = 1:n_cty
        Qprob = q{c} - diag(diag(q{c}));
        
        % Note that q is prob for one event to one event
        % The average number (observed) of children generated by j would be q(i,j)*t(i)
        n_i = covid_tr(c,:).';
        Q = [Q; sum( Qprob.* n_i, 1)];
        
    end
    
    
    %% Estimate K0 and Coefficients in Possion regression
    
    % parameters for possion regression
    opts = statset('glmfit');
    opts.MaxIter = 300;
    
    % boundaty correct
    glm_tr = Covar_tr(1: n_cty*day_for_tr ,:);
    glm_y = Q(:, 1:day_for_tr);
    glm_y = reshape(glm_y, prod(size(glm_y)), 1);
    
    % weight for observation, which is the number of evets at day j
    freqs = covid_tr(:, 1:day_for_tr);
    freqs = reshape(freqs, prod(size(freqs)), 1);
    
    mdl = fitglm( glm_tr, glm_y,'linear', 'Distribution', 'poisson', 'options', opts, 'VarNames', VarNames, 'Weights', freqs);
    
    %% Estimate K0
    %mdl = fitglm( Covar_tr(1: n_cty*day_for_tr ,6:7), glm_y,'linear', 'Distribution', 'poisson', 'options', opts, 'Weights', freqs);
    [ypred,yci] = predict(mdl,Covar_tr);
    K0 = reshape(ypred, n_cty, n_day_tr);

    %Bound K0
    K0 = smoothdata(K0, 2);
    %
    %% Estimate mu, the background rate
    
    for c = 1:n_cty
        mus(c) = sum(( diag(p{c}).' .* covid_tr(c,:) )) / (n_day_tr+muReg) ;
    end
    
    %% Take all the average
    %%mus = repmat( mean(mus), size(mus, 1), size(mus, 2) );
    
    %% Estimate alpha and beta in Weibull Distribution
    
    % Get all pairs
    combos = sortrows( nchoosek( (1:n_day_tr),2 ),2);
    
    % Get those within boundary
    combos = combos(find( prod(combos<=day_for_tr,2)),: );
    combos = combos(:,end:-1:1);
    
    % allocate memory first for spead-up
    obs_sample = zeros( size(combos,1)*n_cty ,1 );
    freq_sample = zeros( size(combos,1)*n_cty ,1 );
    
    for c = 1:n_cty
        
        nc = covid_tr(c,:);
        
        % obs vation is i-j
        obs =[combos(:,1) - combos(:,2)];
        
        % freq
        Prob = p{c};
        freq = Prob(  sub2ind(size(Prob), combos(:,1), combos(:,2) ) ) .* nc(combos(:,1)).';
        if(sum(isnan(freq))>0)
            disp(sum(isnan(freq)))
        end
        % Record obs, freq
        obs_sample( (c-1)*size(combos,1)+1:c*size(combos,1) ) = obs;
        freq_sample( (c-1)*size(combos,1)+1:c*size(combos,1) ) = freq;
        
    end
    
    if (alphaScale_in==0) && (betaShape_in==0)
        % Fit weibull
        [coef,~] = wblfit(obs_sample,[],[],freq_sample);
        alpha=coef(1); % scale
        beta=coef(2);  % shape
        if beta > 100
            beta = 100;
        end
        if alpha > 100
            alpha = 100;
        end
    else
        alpha  = alphaScale_in;
        beta = betaShape_in;
    end
    
    
    
    %% check for the convergence
    if(itr==1)
        % save the first value
        alpha_prev = alpha; beta_prev = beta; mus_prev = mus; K0_prev = K0; theta_prev = mdl.Coefficients.Estimate;
    else
        % Calculate the RMSR
        alpha_delta = [alpha_delta sqrt( (alpha-alpha_prev).^2 )];
        beta_delta  = [beta_delta  sqrt( (beta-beta_prev).^2 )];
        mus_delta   = [mus_delta sqrt( sum((mus_prev - mus).^2)/numel(mus) )];
        K0_delta    = [K0_delta sqrt( sum(sum((K0_prev - K0).^2))/numel(K0) )];
        theta_delta   = [theta_delta sqrt( sum((theta_prev - mdl.Coefficients.Estimate).^2)/numel(mdl.Coefficients.Estimate) )];
        % save the current
        alpha_prev = alpha; beta_prev = beta; mus_prev = mus; K0_prev = K0; theta_prev = mdl.Coefficients.Estimate;
    end
    %disp(max(K0(:)))
    [wblstat_M, wblstat_V]=wblstat(alpha, beta);
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%

    %%%%%%%%%%%%%%%%%%
    % Early Stop
    %%%%%%%%%%%%%%%%%%

    if (itr > 5)
        rule = all( alpha_delta(end-4:end) < break_diff) & all( beta_delta(end-4:end) < break_diff);
        rule = rule &  all( mus_delta(end-4:end) < break_diff) & all( K0_delta(end-4:end) < break_diff);
        rule = rule &  all( theta_delta(end-4:end) < break_diff);

        if( rule )
            break;
        end
    end
end
save(out_save,'mus','alpha','beta','K0','mdl','VarNames','alpha_delta','beta_delta','mus_delta','K0_delta','theta_delta')
end
