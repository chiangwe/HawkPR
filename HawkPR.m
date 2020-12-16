function [] = HawkPR( InputPath_report, InputPath_mobility, InputPath_demography, Delta, Alpha, Beta, EMitr, DaysPred, SimTimes, OutputPath_mdl, OutputPath_pred)
warning('off')

%% Read in parameter
EMitr = EMitr;
if strcmp(Alpha,'') && strcmp(Beta,'')
	disp('No shape and scale parameter for Weibull distribution provided. Use MLE to infer alpha and beta ... ')
    alphaScale_in = 0;
    betaShape_in  = 0;
else
	alphaScale_in = Alpha;
	betaShape_in  = Beta;
end 

if strcmp(Delta,'') 
        disp('No shift parameter for mobility provided.  It will set to zero ... ')
	mobiShift_in = 0;
else
	mobiShift_in = Delta;
end

% Read-in COVID data
NYT = readtable(InputPath_report,'ReadVariableNames',true);

% Read-in mobility
Mobi = readtable( InputPath_mobility,'ReadVariableNames',true);

% Read-in demographic
Demo = readtable(InputPath_demography,'ReadVariableNames',true);
Demo_val = table2array( Demo(:,4:end));

% Data pre-processing
covid = table2array( NYT(:,4:end));
covid(isnan(covid)) = 0;
covid = [zeros(size(covid,1), 1 ) covid(:,2:end) - covid(:,1:end-1)];
covid(covid<=0) = 0;


% Pad to shift 
mob_head = Mobi(:,1:4);
mob_val = table2array(Mobi(:,5:end));

for pad = 1:mobiShift_in
    mob_val = [ mean(mob_val(:,1:7),2) mob_val ];
end

% Get Key and Date

NYT_Date_list = NYT.Properties.VariableNames(4:end);
NYT_Key_list = table2cell(NYT(:,1:3));

Mobi_Type_list = table2cell(Mobi(1:6,4));
Mobi_Date_list = Mobi.Properties.VariableNames(5:end);
Mobi_Key_list = table2cell(Mobi(1:6:end,1:3));

Demo_Type_list = Demo.Properties.VariableNames(4:end);
Demo_Key_list  = table2cell(Demo(:,1));

% Get number of counties and number of days
[n_cty, n_day]=size(covid);
n_mobitype = size(mob_val,1)/n_cty;

disp(['There ' num2str(n_cty) ' counties, ' num2str(n_mobitype) ' types of Mobility indices, and ' num2str(n_day) ' days in the convid reports.' ])

% Train & Test Split
n_tr = size(covid,2);
mob_tr = mob_val(:, 1:n_tr);
mob_te = mob_val(:, n_tr+1:n_tr+DaysPred);


% Normalization
mob_tr_reshape = reshape(mob_tr, n_mobitype, size(mob_tr,1)/n_mobitype * size(mob_tr,2) ).';
mob_te_reshape = reshape(mob_te, n_mobitype, size(mob_te,1)/n_mobitype * size(mob_te,2) ).';
%
Demo_val_in = Demo_val;
Demo_val_tr = repmat(Demo_val_in, n_tr,1);
Demo_val_te = repmat(Demo_val_in, DaysPred,1);

covid_tr = covid;
%
Covar_tr = [mob_tr_reshape Demo_val_tr];
Covar_te = [mob_te_reshape Demo_val_te];
%
Covar_tr_mean = mean(Covar_tr,1);
Covar_tr_std = std(Covar_tr,1);
%
Covar_tr = (Covar_tr-Covar_tr_mean) ./ Covar_tr_std;
Covar_te = (Covar_te-Covar_tr_mean) ./ Covar_tr_std;
%

% Get Variable names
VarNamesOld = [ Mobi_Type_list; Demo_Type_list.'; {['Qprob']}];
VarNames=[];
% Rename
for i = 1:size(VarNamesOld,1)
    newStr = replace( VarNamesOld{i} , ' & ' , '_' );
    newStr = replace( newStr , ' ' , '_' );
    newStr = regexprep(newStr, '^_', '');
    VarNames=[VarNames; {newStr}];
end


%% Define Parameters
n_day_tr = n_day;
T = n_day_tr;
% Boundary correction, the number of days before the total number of days (n_day)
dry_correct = 14;

% EM step iterations
emiter = EMitr; %Nt=covid(:,2:end);
break_diff = 10^-3;
% Boundary correction: T-dry_correct
% Mobility has only 6 weeks so we take the less by min(T-dry_correct, size(mobi_in,2) )
day_for_tr = min(T-dry_correct, size(mob_tr,2) );

%% Initialize Inferred Parameters
if (alphaScale_in==0) && (betaShape_in==0)
    % Weibull Distribution (scale, shape) parameters as (alphas, betas)
    alpha  = 2;  beta = 2;
else
    alpha  = alphaScale_in; beta = betaShape_in;
end

% initial wbl values
wbl_val = repmat(sparse(tril(wblpdf( (1:n_day_tr).' - (1:n_day_tr), alpha, beta))), n_cty, 1);

% K0 reproduction number, a fuction of time and mobility.
% K0 is a n_county * n_day by n_day matrix.
% Estimat for each county at each day.
% each chunk n_day by n_day is a n_day times REPETITION for day j
K0 = ones(n_cty, n_day_tr);
K0_ext_j = repelem( K0, n_day_tr, 1);

% q is a n_county * n_day by n_day matrix.
% q( i, j), i > j stands the probability of ONE SIGLE event at day i triggered by ONE SIGLE event at day j
% i.e., Prob for each j triggered by each i
% each chunk n_day by n_day is for one county
q = sparse(n_cty * n_day_tr,n_day_tr);

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
    tic

    %% E-step

    q = K0_ext_j .* wbl_val .* ( covid_tr_ext_j(covid_tr, n_day_tr) >0);

    eye_mu = repmat( sparse(eye(n_day_tr)) , n_cty, 1) .* repelem(mus, n_day_tr, 1);
    eye_mu = full(eye_mu);

    lam = sum( q .* covid_tr_ext_j(covid_tr, n_day_tr) +  eye_mu, 2);
    lam_eq_zero = lam==0;

    q = q./lam;

    q(lam_eq_zero, :) = 0;

    lam = reshape( lam, n_day_tr, n_cty).';


    %% M-step
    % Calculate Q, which stands for the average number (observed) of children generated by a SINGLE event j
    % Note taht the last "dry_correct" days of Q will be accurate
    % Since we haven't observed their children yet
    Q = reshape( q .* covid_tr_ext_i(covid_tr, n_day_tr, n_cty), n_day_tr, n_day_tr*n_cty);
    Q = reshape( sum(Q,1), n_cty, n_day_tr);
    
    %% Estimate K0 and Coefficients in Possion regression
    %tic
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
    disp(mdl)
    %% Estimate K0
    %mdl = fitglm( Covar_tr(1: n_cty*day_for_tr ,6:7), glm_y,'linear', 'Distribution', 'poisson', 'options', opts, 'Weights', freqs);
    [ypred,yci] = predict(mdl,Covar_tr);
    K0 = reshape(ypred, n_cty, n_day_tr);

    %Bound K0
    K0 = smoothdata(K0, 2);
    K0_ext_j = repelem( K0, n_day_tr, 1);
    %
    %% Estimate mu, the background rate
    %toc
    %disp(toc)

    %tic
    lam_eq_zero = find(lam(:, 1:day_for_tr)==0);
    mus = (mus ./ lam(:, 1:day_for_tr));
    mus(lam_eq_zero) = 0;

    mus = sum(  mus .* covid_tr(:, 1:day_for_tr), 2) / day_for_tr;
    
	if (alphaScale_in==0) && (betaShape_in==0)
        %tic
        obs = sparse( tril( (1:day_for_tr).' - (1:day_for_tr), -1) );

        freq = covid_tr_ext_j(covid_tr, n_day_tr) .* covid_tr_ext_i(covid_tr, n_day_tr, n_cty) .* q;
        freq = full(freq);
        freq = sum( permute( reshape(freq, n_day_tr, n_cty, n_day_tr), [1, 3, 2]), 3);
        freq = freq(1:day_for_tr, 1:day_for_tr);

        Ind_ret = find( (obs > 0) & (freq > 0) ) ;
        obs = obs(Ind_ret);
        freq = freq(Ind_ret);
        %toc
        %disp(toc)
        % Fit weibull
        %tic
        [coef,~] = wblfit(obs,[],[],freq);
        %disp(coef)
        %toc
        %disp(toc)
        alpha=full(coef(1)); % scale
        beta=full(coef(2));  % shape
        wbl_val = repmat(sparse(tril(wblpdf( (1:n_day_tr).' - (1:n_day_tr), alpha, beta))), n_cty, 1);

        %disp(sum(obs.*freq)/sum(freq))
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
            disp(['Convergence Criterion Meet. Break out EM iteration ...'])
            break;
        end
    end
    t = toc;
    disp(['Iterattion ' num2str(itr) ', Elapse time: ' num2str(t)])
end
if(itr == emiter)
    disp(['Reach maximun EM iteration.'])
end
%save(OutputPath_mdl,'mus','alpha','beta','K0','mdl','VarNames','alpha_delta','beta_delta','mus_delta','K0_delta','theta_delta')

% Start Simulation
load(OutputPath_mdl, 'mus','alpha','beta','K0','mdl','VarNames','alpha_delta','beta_delta','mus_delta','K0_delta','theta_delta')

%% Get K0
Covar_all = [Covar_tr; Covar_te];
n_day = n_day_tr+DaysPred;
T_sim = n_day;
Tlow = T_sim-DaysPred;

%% Predict
[ypred,yci] = predict(mdl,Covar_all);
fK0 = reshape(ypred, n_cty, n_day);

% Make fK0 stable
fK0(find(fK0>4))=4;

% Simulation results
sim = zeros(n_cty, T_sim, SimTimes);

% Simulate offsprings
n_per_batch = 10^2;

K0_sim = fK0(:, Tlow+1:end);

for itr = 1:SimTimes
    rng(itr);

    % Calculate base rate
    base = zeros(n_cty, DaysPred);
    n_exh = zeros(n_cty, DaysPred);

    t_stamps = (Tlow+1: T_sim).' - (1:Tlow);
    intense = permute( repmat( wblpdf( t_stamps, alpha, beta), 1, 1, n_cty), [3, 1, 2]).* ...
            permute( repmat( fK0(:, 1:Tlow), 1, 1, DaysPred), [1, 3, 2]) .* ...
            permute( repmat( covid_tr(:, 1:Tlow), 1, 1, DaysPred), [1, 3, 2]);
    base = sum(intense, 3) + mus;
    n_exh = poissrnd( base );
    for itr_cty = 1:ceil(n_cty*0.5)
        for itr_d = 1:DaysPred
            max_d = DaysPred - itr_d;

            % Sample first
            if ( n_exh(itr_cty, itr_d) > n_per_batch)
                n_batch = floor( n_exh(itr_cty, itr_d)  / n_per_batch);
                cand = poissrnd( K0_sim(itr_cty, itr_d), 1, n_per_batch );

                n_mod = mod( n_exh(itr_cty, itr_d), n_per_batch);

                n_offs = sum( cand ) * n_batch + sum( poissrnd( K0_sim(itr_cty, itr_d), 1, n_mod ) );
            else
                n_offs = sum(poissrnd( K0_sim(itr_cty, itr_d), 1, n_exh(itr_cty, itr_d) ));
            end

            if ( n_offs > n_per_batch)
                n_batch = floor( n_offs / n_per_batch);
                n_mod = mod(n_offs, n_per_batch);

                sim_cand_wbl = ceil(wblrnd( alpha, beta, 1, n_per_batch));
                sim_cand_wbl = sim_cand_wbl( find( sim_cand_wbl <= max_d ) );
                sim_cand_wbl =  histcounts( sim_cand_wbl, (1: max_d+1));

                t_delta = ceil(wblrnd( alpha, beta, 1, n_mod));
                t_delta = t_delta( find( t_delta <= max_d ) );
                nt =  histcounts( t_delta, (1: max_d+1)) + sim_cand_wbl * n_batch;

            else
                t_delta = ceil(wblrnd( alpha, beta, 1, n_offs));
                t_delta = t_delta( find( t_delta <= max_d ) );
                nt =  histcounts( t_delta, (1: max_d+1));
            end

            n_exh(itr_cty, itr_d+1:end) = n_exh(itr_cty, itr_d+1:end) + nt;
 
        end
    end
    sim(:,:,itr) = [covid_tr n_exh];
end

sim_out = sim(:,end-DaysPred+1:end,:);

% Formate the output 
sim_mean = mean(sim_out,3);


% Get header 
Date_pred = datetime(datestr( replace(NYT_Date_list(end),'x','') ,'mm/dd/yyyy'))+days(1:DaysPred);
Date_pred = replace(string( datestr(Date_pred,'xyyyy/mm/dd')), '/' ,'_');

table_out  = array2table(sim_mean);
table_out.Properties.VariableNames = Date_pred;
table_out = [NYT(:,1:3) table_out];

writetable(table_out,OutputPath_pred)
end

function [covid_tr_ext_j] = covid_tr_ext_j(covid_tr, n_day_tr)
    covid_tr_ext_j = repelem(covid_tr, n_day_tr, 1);
end

function [covid_tr_ext_i] = covid_tr_ext_i(covid_tr, n_day_tr, n_cty)
    covid_tr_ext_i = repmat( cell2mat(mat2cell(covid_tr.', n_day_tr, ones(1, n_cty) ).'),  1, n_day_tr);
end