 %%
% 8-meter Simply Supported Fully Restrained Steel Beam with UDL
% Four-story steel office building.

% The basic variables include:
% 1.Corrosion rate; 2.Average metal–environment-specific time exponent; 
% 3.Relative humidity; 4.Temperature; 5.Fatigue damage;
% Material properties: 6.Yield strength and 7.Young%s modulus; 
% Sectional geometry: 8.Width, 9.Height, 10.Flange thickness, 11.Web 
% thickness, 12.Beam length, and 13.Beam bay width;
% Applied loads: 14.Dead load, 15.Live load, and 16.Accumulated damage from
% load history.

% Limite State: 
% 1.Shear resistance 2. Moment resistance 3.Deflection of the beam.

% For Time-Dependent Reliability Analysis:
% 1.The thickness of the I-beam's flange and web gradually diminishes due 
%   to corrosion.
% 
% 2.The effects of ambient humidity, temperature, and fatigue on the yield 
%   strength (f_y) and Young's modulus (E) were accounted for.
% 
% 3. For loading: 
%  Dead load imposed on the structure remains constant unless the usage
%     conditions alter.
%
%  Live load, combination of sustained loads and extraordinary loads. For
%    office buildings, live loads are typically considered statistically 
%    stable, with their mean value and coefficient of variation remaining 
%    constant over a 50-year period, unless the building's purpose changes. 
%     
%  But,In this analysis, the modelling was simplified by assuming the live
%  load acting on the building increases linearly by 0.25% annually. 
%  Consequently, at the 50-year, the live load will be increased by 12.5%.
%  
%  This approach simplifies the simulation of damage to the structure 
%  caused by extraordinary loads during its service life, thereby reducing 
%  the structure's reliability.

% For System Reliability:
% The correlation between the three failure modes was calculated, and by 
% treating the three failure modes as a series system, the system reliability 
% of the beam was determined through Monte Carlo simulation, along with its
% 95% confidence interval. Additionally, the System Failure Probability 
% Lower and Upper Bounds were verified using the Ditlevsen method."
 
% For Sensitivity Analysis:
% Quantified the influence of 3 key parameters (yield strength fy, corrosion 
% rate,and live load LL) on structural reliability over a 50-year horizon.
% Using a 5% perturbation method, computed normalized sensitivity coefficients 
% representing the change in reliability index (Δβ) per 1% parameter variation. 
% Critical findings revealed:
%
% Live load (mu_LL) consistently exhibited the strongest negative influence
% (up to -4.64 Δβ/% at t=50 years), indicating load increases accelerate 
% reliability degradation
%
% Yield strength (mu_fy) showed significant positive impact (+5.62 to +8.22
% Δβ/%), demonstrating material strength's protective role.
%
% Corrosion rate (mu_kcor1) effects were negligible (<0.08 Δβ/%) with 
% unstable sign, confirming indoor corrosion's minimal impact

% For Lifetime Prediction:
% Employed quadratic regression of time-dependent reliability indices 
% (β = a - bt - ct²) to forecast service life against 4 different target 
% reliability indices: 'High Safety (β=3.8)', 'Medium Safety (β=3.1)', 
% 'Low Safety (β=2.3), 'Critical (β=1.5)'.
% Delivering individual lifetime predictions for each failure mode with 95% 
% confidence interval calculations
%
% For Ultimate limit stat (ULS):EN 1990 recommends 4.2(CC1) to 4.7(CC2) to 
% 5.2(CC3) for 50-year reference period, depending on the consequence level.
% For Serviceability limit state (SLS),ISO 2394 recommends 3.8 for 50 years
% and 1.5 for irreversible case, 50-year reference period.
% target beta = 3.1 and 2.3 are choosen for reference only.

%%
% Monte Carlo Simulation Method with Integrated Lifetime Prediction
nSamples = 2e7;
current_age = 30;
design_life = 50;
t_points = 0:5:design_life; % Time points every 5 years
n_time = length(t_points);

%%
% Initialize arrays to store reliability indices
beta_shear = zeros(n_time, 1); Pf_shear = zeros(n_time, 1);
beta_moment = zeros(n_time, 1); Pf_moment = zeros(n_time, 1);
beta_deflection = zeros(n_time, 1); Pf_deflection = zeros(n_time, 1);
beta_system = zeros(n_time, 1); Pf_system = zeros(n_time, 1);
% Initialize arrays for storing all samples at each time point
Med_all = zeros(nSamples, n_time);
Mrd_all = zeros(nSamples, n_time);
Ved_all = zeros(nSamples, n_time);
Vrd_all = zeros(nSamples, n_time);
delta_all = zeros(nSamples, n_time);
delta_limit_all = zeros(nSamples, n_time);
% Initialize arrays for Ditlevsen bounds and confidence intervals
Pf_lower = zeros(n_time, 1);
Pf_upper = zeros(n_time, 1);
ci_low = zeros(n_time, 1);
ci_high = zeros(n_time, 1);
% Initialize results structure
results = struct();
% Initialize arrays for material properties tracking
mean_tw = zeros(n_time, 1);
mean_tf = zeros(n_time, 1);
mean_A = zeros(n_time, 1);
mean_fy = zeros(n_time, 1);
mean_E = zeros(n_time, 1);

%%
% Aging and Degradation
 mu_kcor1 =  0.0013; % Interior corrosion rate (mm/year) (BS EN ISO 9223(2012)) FOR C1
% mu_kcor1 = 0.0375; % For C3 corrosion level ####
cov_kcor1 = 0.3;
mu_B = 0.523; % Average metal-environment-specific time exponent (BS EN ISO 9223(2012) Table.3)
sd_B = 0.026;  cov_B = sd_B / mu_B;  % ISO CORRAG 
c = 0;     % Assume No cover (mm)

% Environment
RH = 0.55 + 0.2*rand(nSamples,1);   % Relative humidity 55%-75%
time_vec = (1:nSamples)';
T = 20 + 7*cos(2*pi*time_vec/365) + 2*randn(nSamples,1);   % Temperature
D_fat = 0.0005 * lognrnd(0,1,nSamples,1);   % Random fatigue damage

% Material
mu_fy = 275;      % Yield strength (MPa), Lognormal
cov_fy = 0.08;
mu_E = 210000;    % Young's modulus (MPa), Lognormal
cov_E = 0.03;

% Cross-section geometry (457 x 191 x 82 UKB, S275)
mu_b = 191.3;     % Width (mm), Normal
cov_b = 0.01;
mu_h = 460;       % Height (mm), Normal
cov_h = 0.01;
mu_tf = 16;       % Flange thickness (mm), Normal
cov_tf = 0.03;
mu_tw = 9.9;      % Web thickness (mm), Normal
cov_tw = 0.04;

% Loading
mu_L = 8000;      % Beam length (mm), Normal
cov_L = 0.005;
mu_d = 6000;      % Beam bay width (mm), Normal
cov_d = 0.005;
mu_DL = 3.7e-3;   % Dead load (N/mm^2), Lognormal
cov_DL = 0.15;
mu_LL = 3.3e-3;   % Live load (N/mm^2), Gumbel
cov_LL = 0.25;
increase_rate = 0.0025;   % 0.25% increase per year
cov_LL_final = 0.3;     % Final COV at 50 years
increase_rate_cov = (cov_LL_final - cov_LL) / 50;

% Sampling Functions
logn_sample = @(mu, cov, n) exp(log(mu) - 0.5*log(1 + cov^2) + sqrt(log(1 + cov^2))*randn(n,1));
norm_sample = @(mu, cov, n) mu + (mu*cov)*randn(n,1);

% Generate time-independent samples
samples_fy = logn_sample(mu_fy, cov_fy, nSamples);
samples_E = logn_sample(mu_E, cov_E, nSamples);
samples_b = norm_sample(mu_b, cov_b, nSamples);
samples_h = norm_sample(mu_h, cov_h, nSamples);
samples_tf = norm_sample(mu_tf, cov_tf, nSamples);
samples_tw = norm_sample(mu_tw, cov_tw, nSamples);
samples_L = norm_sample(mu_L, cov_L, nSamples);
samples_d = norm_sample(mu_d, cov_d, nSamples);
samples_DL = logn_sample(mu_DL, cov_DL, nSamples);
samples_kcor1 = logn_sample(mu_kcor1, cov_kcor1, nSamples);
samples_B = norm_sample(mu_B, cov_B, nSamples);
% Corrosion initiation times (fixed: ensure non-negative with proper handling when c=0)
if c == 0
    t0 = zeros(nSamples, 1); % No delay when no cover
else
    t0 = max(0, (c ./ samples_kcor1)); % Linear relationship when c>0
end

%% Time loop for reliability assessment
fprintf('Running Monte Carlo simulation...\n');
for i = 1:n_time
    t = t_points(i);
    fprintf('Processing time point: %d years\n', t);

    % Time-dependent live load (Gumbel Distribution)
    mu_LL_t = mu_LL * (1 + increase_rate * t);
    cov_LL_t = cov_LL + increase_rate_cov * t;
    scale_LL = (sqrt(6) * mu_LL_t * cov_LL_t) / pi;
    loc_LL = mu_LL_t - 0.5772 * scale_LL;
    samples_LL = evrnd(loc_LL, scale_LL, nSamples, 1);
    
    % Corrosion effects
    effective_time = max(0, t - t0); % Effective component corrosion time
    if effective_time > 0
        if effective_time <= 20
            % Power-law model for t ≤ 20 years
             D = samples_kcor1 .* (effective_time.^samples_B);
        else
            % Conservative linear model for t > 20 years
            D_20 = samples_kcor1 .* (20.^samples_B);  % Corrosion at 20 years
            k = samples_B .* (20.^(samples_B-1));      % Linear corrosion coefficient
            D = D_20 + samples_kcor1 .* k .* (effective_time - 20);
        end
        thickness_loss = 2 * D;  % Corrosion from both sides
    else
        thickness_loss = zeros(nSamples, 1);
    end
    tw_current = max(0.1, samples_tw - thickness_loss); % Minimum thickness 0.1mm
    tf_current = max(0.1, samples_tf - thickness_loss); % Minimum thickness 0.1mm
        
    % Section properties
    hw = samples_h - 2*tf_current;
    A = samples_b .* samples_h - (samples_b - samples_tw) .* (samples_h - 2*samples_tf);
    A_web = tw_current .* hw; 
    A_flange = 2 * samples_b .* tf_current; 
    Anet = A_web + A_flange; 

    % Material degradation/ fy degradation 
    alpha_corr = 0.005;   % Steel degradation factor (ACI 222R)
    fy_corr = samples_fy .* (1 - alpha_corr * (A - Anet) ./ A);

    % Strength loss due to fatigue (Miner's criterion)
    beta_fat = 0.6;   % Fatigue degradation index (Eurocode 3-1-9) 
    fy_fat = fy_corr .* (1 - D_fat).^beta_fat;

    % Long-term aging effects (exponential decay)
    k_aging = 1e-5; % ageing factor (Papakonstantinou 2021),(Corrosion Science, 2021)
    fy_t = fy_fat .* exp(-k_aging * t);

    % E degradation
    gamma = 0.001;  % stiffness degradation factor (Model Code 2010)
    E_corr = samples_E .* (1 - gamma * (A - Anet) ./ A);  % Corrosion impact
    Q = 35e3;   % Activation energy
    R = 8.314;  % Gas constant
    k_env = exp(-Q./(R*(T+273))) .* RH.^1.5;  % Environmental factor
    E_env = E_corr .* exp(-0.0001 * k_env .* t);  % Time-environment degradation, γ: stiffness degradation coefficient (0.0001)
    E_t = E_env .* (1 - 0.05 * D_fat);  % Fatigue impact
    
    % Second moment of area
    Iy_web = tw_current .* hw.^3 / 12;
    Iy_flange = 2 * (samples_b .* tf_current.^3 / 12 + samples_b .* tf_current .* (samples_h/2 - tf_current/2).^2);
    Iy = Iy_web + Iy_flange;
    
    % Plastic modulus
    Wpl = samples_b .* tf_current .* (samples_h - tf_current) + 0.25 * tw_current .* (samples_h - 2*tf_current).^2;
    
    % Loads and resistances
    Ved = (samples_DL + samples_LL) .* samples_d .* samples_L ./ 2;
    Vrd = fy_t .* samples_h .* tw_current ./ sqrt(3);
    Med = (samples_DL + samples_LL) .* samples_d .* samples_L.^2 ./ 8;
    Mrd = fy_t .* Wpl;
    deflection = 5 .* (samples_LL) .* samples_d .* (samples_L.^4) ./ (384 * E_t .* Iy);
    delta_limit = samples_L ./ 360;

    % Store values
    Ved_all(:, i) = Ved;
    Vrd_all(:, i) = Vrd;
    Med_all(:, i) = Med;
    Mrd_all(:, i) = Mrd;
    delta_all(:, i) = deflection;
    delta_limit_all(:, i) = delta_limit;
    
    % Store material properties
    mean_tw(i) = mean(tw_current);
    mean_tf(i) = mean(tf_current);
    mean_A(i) = mean(Anet);
    mean_fy(i) = mean(fy_t);
    mean_E(i) = mean(E_t); 

    % Limit state functions
    g_shear = Vrd - Ved;
    g_bending = Mrd - Med;
    g_deflection = delta_limit - deflection;
    
    % Individual reliability indices
    Pf_shear_temp = mean(g_shear <= 0);
    if Pf_shear_temp < 1e-7 || isnan(Pf_shear_temp)
    beta_shear(i) = 10;
    else
    beta_shear(i) = max(-10, min(10, -norminv(Pf_shear_temp)));
    end
    Pf_shear(i) = Pf_shear_temp;

    Pf_moment_temp = mean(g_bending <= 0);
    if Pf_moment_temp < 1e-7 || isnan(Pf_moment_temp)
    beta_moment(i) = 10;
    else
    beta_moment(i) = max(-10, min(10, -norminv(Pf_moment_temp)));
    end
    Pf_moment(i) = Pf_moment_temp;

    Pf_deflection_temp = mean(g_deflection <= 0);
    if Pf_deflection_temp < 1e-7 || isnan(Pf_deflection_temp)
    beta_deflection(i) = 10;
    else
    beta_deflection(i) = max(-10, min(10, -norminv(Pf_deflection_temp)));
    end
    Pf_deflection(i) = Pf_deflection_temp;

    % Correlation between failure modes
    G = [g_shear, g_bending, g_deflection];
    rho_failure_modes = corrcoef(G);
    
    % Store correlation matrix results
    results(i).time = t_points(i);
    results(i).rho_matrix = rho_failure_modes;
    results(i).mean_g = mean(G);
    results(i).std_g = std(G);
    
    % Generate correlated random variables
    betas_current = [beta_shear(i), beta_moment(i), beta_deflection(i)];
    R = mvnrnd(zeros(1, 3), rho_failure_modes, nSamples);
    
    % Transform to standard normal space
    g_std = betas_current - R;
    
    % Series system failure probability
    system_failure = any(g_std < 0, 2);
    Pf_system(i) = mean(system_failure);
    
    % Calculate system reliability index
    if Pf_system(i) < 1e-7 || isnan(Pf_system(i))
        beta_system(i) = 10;
    elseif Pf_system(i) >= 1
        beta_system(i) = -10;
    else
        beta_system(i) = max(-10, min(10, -norminv(Pf_system(i))));
    end
    
    % Calculate Ditlevsen bounds
    [Pf_lower(i), Pf_upper(i)] = ditlevsen_bounds(betas_current, rho_failure_modes);
    
    % 95% confidence interval
    se = sqrt(Pf_system(i)*(1-Pf_system(i))/nSamples);
    ci_low(i) = max(0, Pf_system(i) - 1.96*se);
    ci_high(i) = min(1, Pf_system(i) + 1.96*se);
end

%% Lifetime Prediction Analysis
fprintf('\nPerforming lifetime prediction analysis...\n');

% Define target reliability indices for different scenarios
target_betas = [3.8, 3.1, 2.3]; % High, Medium, Low reliability targets
target_names = {'High Safety (β=3.8)', 'Medium Safety (β=3.1)', 'Low Safety (β=2.3)'};

% Initialize lifetime prediction results
lifetime_results = struct();
lifetime_results.target_betas = target_betas;
lifetime_results.target_names = target_names;
lifetime_results.shear_lifetimes = zeros(length(target_betas), 1);
lifetime_results.moment_lifetimes = zeros(length(target_betas), 1);
lifetime_results.deflection_lifetimes = zeros(length(target_betas), 1);
lifetime_results.system_lifetimes = zeros(length(target_betas), 1);
lifetime_results.shear_ci = zeros(length(target_betas), 2);
lifetime_results.moment_ci = zeros(length(target_betas), 2);
lifetime_results.deflection_ci = zeros(length(target_betas), 2);
lifetime_results.system_ci = zeros(length(target_betas), 2);

% Predict lifetimes for each failure mode and target reliability
for j = 1:length(target_betas)
    target_beta = target_betas(j);
    
    % Shear failure mode
    [life_shear, ci_shear] = predict_lifetime(t_points, beta_shear, target_beta);
    lifetime_results.shear_lifetimes(j) = life_shear;
    lifetime_results.shear_ci(j, :) = ci_shear;
    
    % Moment failure mode
    [life_moment, ci_moment] = predict_lifetime(t_points, beta_moment, target_beta);
    lifetime_results.moment_lifetimes(j) = life_moment;
    lifetime_results.moment_ci(j, :) = ci_moment;
    
    % Deflection failure mode
    [life_deflection, ci_deflection] = predict_lifetime(t_points, beta_deflection, target_beta);
    lifetime_results.deflection_lifetimes(j) = life_deflection;
    lifetime_results.deflection_ci(j, :) = ci_deflection;
    
    % System failure mode
    [life_system, ci_system] = predict_lifetime(t_points, beta_system, target_beta);
    lifetime_results.system_lifetimes(j) = life_system;
    lifetime_results.system_ci(j, :) = ci_system;
    
    fprintf('Target β = %.1f: System lifetime = %.1f years (95%% CI: %.1f-%.1f years)\n', ...
        target_beta, life_system, ci_system(1), ci_system(2));
end

%% Sensitivity Analysis - Enhanced with validation
fprintf('\nPerforming sensitivity analysis...\n');
parameters = {'mu_fy', 'mu_kcor1', 'mu_LL'};
delta = 0.05; 
beta_system_base = beta_system; % Store baseline beta_system
sensitivity = zeros(length(parameters), n_time); % Sensitivity matrix

% Store original parameter values
mu_fy_orig = mu_fy;
mu_kcor1_orig = mu_kcor1;
mu_LL_orig = mu_LL;

for p = 1:length(parameters)
    param_name = parameters{p};
    fprintf('Analyzing sensitivity to %s...\n', param_name);
      
    % Get base parameter value
    if strcmp(param_name, 'mu_fy')
        param_base = mu_fy_orig;
        mu_fy = param_base * (1 + delta);
        samples_fy_sens = logn_sample(mu_fy, cov_fy, nSamples);
    elseif strcmp(param_name, 'mu_kcor1')
        param_base = mu_kcor1_orig;
        mu_kcor1 = param_base * (1 + delta);
        samples_kcor1_sens = logn_sample(mu_kcor1, cov_kcor1, nSamples);
        if c == 0
            t0_sens = zeros(nSamples, 1);
        else
            t0_sens = max(0, (c ./ samples_kcor1_sens));
        end
   elseif strcmp(param_name, 'mu_LL')
        param_base = mu_LL_orig;
        mu_LL = param_base * (1 + delta);
    end
    
    % Re-run Monte Carlo for each time point with perturbed parameter
    for i = 1:n_time
        t = t_points(i);
        
        % Time-dependent live load
        mu_LL_t = mu_LL * (1 + increase_rate * t);
        cov_LL_t = cov_LL + increase_rate_cov * t;
        scale_LL = (sqrt(6) * mu_LL_t * cov_LL_t) / pi;
        loc_LL = mu_LL_t - 0.5772 * scale_LL;
        samples_LL_sens = evrnd(loc_LL, scale_LL, nSamples, 1);
        
        % Use appropriate samples based on parameter being varied
        if strcmp(param_name, 'mu_fy')
            samples_fy_current = samples_fy_sens;
            samples_kcor1_current = samples_kcor1;
            t0_current = t0;
        elseif strcmp(param_name, 'mu_kcor1')
            samples_fy_current = samples_fy;
            samples_kcor1_current = samples_kcor1_sens;
            t0_current = t0_sens;
        else % mu_LL case
            samples_fy_current = samples_fy;
            samples_kcor1_current = samples_kcor1;
            t0_current = t0;
        end
        
        % Corrosion effect
        effective_time = max(0, t - t0); % Effective component corrosion time
    if effective_time > 0
        if effective_time <= 20
            % Power-law model for t ≤ 20 years
             D = samples_kcor1 .* (effective_time.^samples_B);
        else
            % Conservative linear model for t > 20 years
            D_20 = samples_kcor1 .* (20.^samples_B);  % Corrosion at 20 years
            k = samples_B .* (20.^(samples_B-1));            % Linear corrosion coefficient
            D = D_20 + samples_kcor1 .* k .* (effective_time - 20);
        end
        thickness_loss = 2 * D;  % Corrosion from both sides
    else
        thickness_loss = zeros(nSamples, 1);
    end
        tw_current = max(0.1, samples_tw - thickness_loss);
        tf_current = max(0.1, samples_tf - thickness_loss);
        
        % Section properties
        hw = samples_h - 2*tf_current;
        A = samples_b .* samples_h - (samples_b - samples_tw) .* (samples_h - 2*samples_tf);
        A_web = tw_current .* hw; 
        A_flange = 2 * samples_b .* tf_current; 
        Anet = A_web + A_flange; 

        % Material degradation
        alpha_corr = 0.005;
        fy_corr = samples_fy_current .* (1 - alpha_corr * (A - Anet) ./ A);
        beta_fat = 0.6;
        fy_fat = fy_corr .* (1 - D_fat).^beta_fat;
        k_aging = 1e-5;
        fy_t = fy_fat .* exp(-k_aging * t);
        
        % E degradation
        gamma = 0.001;
        E_corr = samples_E .* (1 - gamma * (A - Anet) ./ A);
        Q = 35e3;
        R = 8.314;
        k_env = exp(-Q./(R*(T+273))) .* RH.^1.5;
        E_env = E_corr .* exp(-0.0001 * k_env .* t);
        E_t = E_env .* (1 - 0.05 * D_fat);
        
        % Second moment of area
        Iy_web = tw_current .* hw.^3 / 12;
        Iy_flange = 2 * (samples_b .* tf_current.^3 / 12 + samples_b .* tf_current .* (samples_h/2 - tf_current/2).^2);
        Iy = Iy_web + Iy_flange;
        
        % Plastic modulus
        Wpl = samples_b .* tf_current .* (samples_h - tf_current) + 0.25 * tw_current .* (samples_h - 2*tf_current).^2;
        
        % Loads and resistances
        Ved = (samples_DL + samples_LL_sens) .* samples_d .* samples_L ./ 2;
        Vrd = fy_t .* samples_h .* tw_current ./ sqrt(3);
        Med = (samples_DL + samples_LL_sens) .* samples_d .* samples_L.^2 ./ 8;
        Mrd = fy_t .* Wpl;
        delta_calc = 5 .* (samples_LL_sens) .* samples_d .* (samples_L.^4) ./ (384 * E_t .* Iy);
        delta_limit = samples_L ./ 360;
        
        % Limit state functions
        g_shear = Vrd - Ved;
        g_bending = Mrd - Med;
        g_deflection = delta_limit - delta_calc;
        
        % Check for numerical issues
        valid_idx = ~(isinf(g_shear) | isnan(g_shear) | isinf(g_bending) | isnan(g_bending) | isinf(g_deflection) | isnan(g_deflection));
        
        % Individual reliability indices for sensitivity case
        Pf_shear_sens = mean(g_shear(valid_idx) <= 0);
        if Pf_shear_sens < 1e-7 || isnan(Pf_shear_sens)
           beta_shear_sens = 10;
        else
           beta_shear_sens = max(-10, min(10, -norminv(Pf_shear_sens)));
        end
        
        Pf_moment_sens = mean(g_bending(valid_idx) <= 0);
        if Pf_moment_sens < 1e-7 || isnan(Pf_moment_sens)
            beta_moment_sens = 10;
        else
            beta_moment_sens = max(-10, min(10, -norminv(Pf_moment_sens)));
        end
        
        Pf_deflection_sens = mean(g_deflection(valid_idx) <= 0);
        if Pf_deflection_sens < 1e-7 || isnan(Pf_deflection_sens)
            beta_deflection_sens = 10;
        else
            beta_deflection_sens =  max(-10, min(10, -norminv(Pf_deflection_sens)));
        end
        
        % Correlation between failure modes
        G_sens = [g_shear(valid_idx), g_bending(valid_idx), g_deflection(valid_idx)];
        rho_failure_modes_sens = corrcoef(G_sens);
        
        % Generate correlated random variables
        betas_current_sens = [beta_shear_sens, beta_moment_sens, beta_deflection_sens];
        R_sens = mvnrnd(zeros(1, 3), rho_failure_modes_sens, sum(valid_idx));
        
        % Transform to standard normal space
        g_std_sens = betas_current_sens - R_sens;
        
        % Series system failure probability
        system_failure_sens = any(g_std_sens < 0, 2);
        Pf_system_sens = mean(system_failure_sens);
        
        % Calculate system reliability index for sensitivity case
        if Pf_system_sens < 1e-7 || isnan(Pf_system_sens)
            beta_system_sens = 10;
        elseif Pf_system_sens >= 1
            beta_system_sens = -10;
        else
            beta_system_sens = max(-10, min(10, -norminv(Pf_system_sens)));
        end

         if abs(beta_system_base(i)) < 1e-5
            sensitivity(p, i) = 0; 
        else
            % S = [Δβ/β_base] / [Δθ/θ_base]
            relative_beta_change = (beta_system_sens - beta_system_base(i)) / beta_system_base(i);
            relative_param_change = delta; % Since θ_new = θ_base*(1+delta)
            sensitivity(p, i) = relative_beta_change / relative_param_change;
        end
    end
    % Reset parameters to original values
    mu_fy = mu_fy_orig;
    mu_kcor1 = mu_kcor1_orig;
    mu_LL = mu_LL_orig;
end


    %% Result check:
mu_A = mean(A);var_A = var(A);sd_A = std(A);
mu_Iy = mean(Iy);var_Iy = var(Iy);sd_Iy = std(Iy);
mu_Wpl = mean(Wpl);var_Wpl = var(Wpl);sd_Wpl = std(Wpl);
mu_Ved = mean(Ved);var_Ved = var(Ved);sd_Ved = std(Ved);
mu_Vrd = mean(Vrd);var_Vrd = var(Vrd);sd_Vrd = std(Vrd);
mu_Med = mean(Med);var_Med = var(Med);sd_Med = std(Med);
mu_Mrd = mean(Mrd);var_Mrd = var(Mrd);sd_Mrd = std(Mrd);
mu_delta = mean(delta);var_delta = var(delta);sd_delta = std(delta);
mu_delta_limit = mean(delta_limit);var_delta_limit = var(delta_limit);sd_delta_limit = std(delta_limit);

%% Result Print

% Ensure all arrays are of consistent length
n_time = length(t_points);

% Check and fix array lengths
if length(beta_shear) < n_time
    fprintf('Fix: Length of beta_shear (%d) < number of time points (%d), padding with NaN\n', length(beta_shear), n_time);
    beta_shear = [beta_shear; nan(n_time - length(beta_shear), 1)];
end

if length(Pf_shear) < n_time
    fprintf('Fix: Length of Pf_shear (%d) < number of time points (%d), padding with NaN\n', length(Pf_shear), n_time);
    Pf_shear = [Pf_shear; nan(n_time - length(Pf_shear), 1)];
end

% Repeat the same check for other arrays
arrays_to_check = {'beta_moment', 'Pf_moment', 'beta_deflection', 'Pf_deflection', ...
                   'Pf_system', 'beta_system', 'ci_low', 'ci_high', 'Pf_lower', 'Pf_upper'};
for i = 1:length(arrays_to_check)
    var_name = arrays_to_check{i};
    if length(eval(var_name)) < n_time
        fprintf('Fix: Length of %s (%d) < number of time points (%d), padding with NaN\n',...
                var_name, length(eval(var_name)), n_time);
        eval([var_name ' = [' var_name '; nan(n_time - length(' var_name '), 1)];']);
    end
end

% Output results at each time point
for i = 1:n_time
    fprintf('\n\n===== Time Point: %d years =====\n', t_points(i));
    
    % ===== 1. Output basic reliability indicators =====
    fprintf('\n=== Probability of failure and Reliability Index ===\n');
    
    % Shear limit state
    if isnan(beta_shear(i)) || isnan(Pf_shear(i))
        fprintf('Shear: Data missing\n');
    else
        if Pf_shear(i) < 1e-7
            fprintf('Shear: β = %.4f (Pf < 1e-7)\n', beta_shear(i));
        else
            fprintf('Shear: β = %.4f (Pf = %.4e)\n', beta_shear(i), Pf_shear(i));
        end
    end
    
   % Moment limit state
    if isnan(beta_moment(i)) || isnan(Pf_moment(i))
        fprintf('Moment: Data missing\n');
    else
        if Pf_moment(i) < 1e-7
            fprintf('Moment: β = %.4f (Pf < 1e-7)\n', beta_moment(i));
        else
            fprintf('Moment: β = %.4f (Pf = %.4e)\n', beta_moment(i), Pf_moment(i));
        end
    end
    
    % Deflection limit state
    if isnan(beta_deflection(i)) || isnan(Pf_deflection(i))
        fprintf('Displacement: Data missing\n');
    else
        if Pf_deflection(i) < 1e-7
            fprintf('Displacement: β = %.4f (Pf < 1e-7)\n', beta_deflection(i));
        else
            fprintf('Displacement: β = %.4f (Pf = %.4e)\n', beta_deflection(i), Pf_deflection(i));
        end
    end
    
    % ===== 2. Output correlation matrix =====
    if i <= length(results) && isfield(results(i), 'rho_matrix')
        failure_modes = {'Shear', 'Moment', 'Displacement'};
        n_modes = length(failure_modes);
        
        fprintf('\n=== Failure Modes Correlations Matrix ===\n');
        fprintf('%15s', ' ');  
        for j = 1:n_modes
            fprintf('%15s', failure_modes{j});
        end
        fprintf('\n%s\n', repmat('-', 1, 15*(n_modes+1)));
        
        rho_matrix = results(i).rho_matrix;
        for j = 1:n_modes
            fprintf('%15s', failure_modes{j});  
            for k = 1:n_modes
                if k < j
                    fprintf('%15s', ' ');
                else
                    if j <= size(rho_matrix, 1) && k <= size(rho_matrix, 2)
                        fprintf('%15.4f', rho_matrix(j, k));
                    else
                        fprintf('%15s', 'N/A');
                    end
                end
            end
            fprintf('\n');
        end
    else
        fprintf('\n=== Warning: Correlation matrix data missing===\n');
    end
    
    % ===== 3. Output system reliability results =====
    fprintf('\n=== System Reliability Analysis Results ===\n');
    
    % Monte Carlo results
    if isnan(Pf_system(i))
        fprintf('System Failure Probability (Monte Carlo): Data Lost\n');
    else
        fprintf('System Failure Probability (Monte Carlo): %.4e\n', Pf_system(i));
    end
    
    % Confidence interval
    if isnan(ci_low(i)) || isnan(ci_high(i))
        fprintf('95%% Confidence Interval: Data Lost\n');
    else
        fprintf('95%% Confidence Interval: [%.4e, %.4e]\n', ci_low(i), ci_high(i));
    end
    
    % System reliability index
    if isnan(beta_system(i))
        fprintf('System Reliability Index: Data Lost\n');
    else
        fprintf('System Reliability Index: β_system = %.4f\n', beta_system(i));
    end
    
    % ===== Ditlevsen Boundary =====
    fprintf('\n=== Ditlevsen Boundary Analysis ===\n');
    
    if isnan(Pf_lower(i))
        fprintf('Lower Bound: Data Lost\n');
    else
        fprintf('Lower Bound: %.4e\n', Pf_lower(i));
    end
    
    if isnan(Pf_upper(i))
        fprintf('Upper Bound: Data Lost\n');
    else
        fprintf('Upper Bound: %.4e\n', Pf_upper(i));
    end
    
    % ===== 5. Degradation trend analysis =====
    if i > 1
        time_interval = t_points(i) - t_points(i-1);
        
        % Change rate of system failure probability
        if ~isnan(Pf_system(i)) && ~isnan(Pf_system(i-1))
            Pf_change = (Pf_system(i) - Pf_system(i-1)) / time_interval;
            fprintf('\nAnnual ΔPf_system: %+.4e/year\n', Pf_change);
        end
        
        % Change rate of system reliability index
        if ~isnan(beta_system(i)) && ~isnan(beta_system(i-1))
            beta_change = (beta_system(i) - beta_system(i-1)) / time_interval;
            fprintf('Annual Δβ_system: %+.4f/year\n', beta_change);
        end
    end
end

% Save diagnostic information
diagnostic_info = struct();
diagnostic_info.t_points = t_points;
diagnostic_info.array_lengths = struct(...
    'beta_shear', length(beta_shear), ...
    'Pf_shear', length(Pf_shear), ...
    'beta_moment', length(beta_moment), ...
    'Pf_moment', length(Pf_moment), ...
    'beta_deflection', length(beta_deflection), ...
    'Pf_deflection', length(Pf_deflection), ...
    'Pf_system', length(Pf_system), ...
    'beta_system', length(beta_system), ...
    'ci_low', length(ci_low), ...
    'ci_high', length(ci_high), ...
    'Pf_lower', length(Pf_lower), ...
    'Pf_upper', length(Pf_upper), ...
    'results', length(results) ...
);


%% Detailed Results Summary
fprintf('\n========== COMPREHENSIVE ANALYSIS RESULTS ==========\n');
fprintf('Monte Carlo Samples: %d\n', nSamples);
fprintf('Design Life: %d years\n', design_life);
fprintf('\n--- RELIABILITY INDICES AT KEY TIME POINTS ---\n');
fprintf('Time (years)\tShear β\t\tMoment β\tDeflection β\tSystem β\n');
for i = 1:n_time
    fprintf('%d\t\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\n', ...
        t_points(i), beta_shear(i), beta_moment(i), beta_deflection(i), beta_system(i));
end

fprintf('\n--- LIFETIME PREDICTIONS ---\n');
for j = 1:length(target_betas)
    fprintf('\n%s:\n', target_names{j});
    fprintf('  Shear failure mode: %.1f years (95%% CI: %.1f-%.1f)\n', ...
        lifetime_results.shear_lifetimes(j), lifetime_results.shear_ci(j,1), lifetime_results.shear_ci(j,2));
    fprintf('  Moment failure mode: %.1f years (95%% CI: %.1f-%.1f)\n', ...
        lifetime_results.moment_lifetimes(j), lifetime_results.moment_ci(j,1), lifetime_results.moment_ci(j,2));
    fprintf('  Deflection failure mode: %.1f years (95%% CI: %.1f-%.1f)\n', ...
        lifetime_results.deflection_lifetimes(j), lifetime_results.deflection_ci(j,1), lifetime_results.deflection_ci(j,2));
    fprintf('  SYSTEM FAILURE: %.1f years (95%% CI: %.1f-%.1f)\n', ...
        lifetime_results.system_lifetimes(j), lifetime_results.system_ci(j,1), lifetime_results.system_ci(j,2));
end

% Display sensitivity results with physical interpretation
fprintf('\n=== SENSITIVITY ANALYSIS RESULTS ===\n');
fprintf('Sensitivity coefficients (Δβ per %% parameter change):\n');
fprintf('Expected signs: mu_fy(+), mu_kcor1(-), mu_LL(-)\n\n');
fprintf('Parameter\t');
for i = 1:n_time
    fprintf('t=%d\t\t', t_points(i));
end
fprintf('\n');

for p = 1:length(parameters)
    fprintf('%s\t\t', parameters{p});
    for i = 1:n_time
        if strcmp(parameters{p}, 'mu_kcor1') && sensitivity(p, i) > 0
            fprintf('%.6f*\t', sensitivity(p, i)); % Mark unexpected positive values
        elseif strcmp(parameters{p}, 'mu_fy') && sensitivity(p, i) < 0
            fprintf('%.6f*\t', sensitivity(p, i)); % Mark unexpected negative values
        elseif strcmp(parameters{p}, 'mu_LL') && sensitivity(p, i) > 0
            fprintf('%.6f*\t', sensitivity(p, i)); % Mark unexpected positive values
        else
            fprintf('%.6f\t', sensitivity(p, i));
        end
    end
    fprintf('\n');
end
fprintf('* indicates potentially unexpected sign\n');

fprintf('\n--- MATERIAL DEGRADATION SUMMARY ---\n');
fprintf('Property degradation over %d years:\n', design_life);
fprintf('  Web thickness: %.2f mm → %.2f mm (%.1f%% loss)\n', ...
    mean_tw(1), mean_tw(end), (mean_tw(1)-mean_tw(end))/mean_tw(1)*100);
fprintf('  Flange thickness: %.2f mm → %.2f mm (%.1f%% loss)\n', ...
    mean_tf(1), mean_tf(end), (mean_tf(1)-mean_tf(end))/mean_tf(1)*100);
fprintf(' Cross-section area: %.2f mm^2 → %.2f mm^2 (%.1f%% loss)\n', ...
    mean_A(1), mean_A(end), (mean_A(1)-mean_A(end))/mean_A(1)*100);
fprintf('  Yield strength: %.0f MPa → %.0f MPa (%.1f%% loss)\n', ...
    mean_fy(1), mean_fy(end), (mean_fy(1)-mean_fy(end))/mean_fy(1)*100);
fprintf('  Youngs modulus: %.0f MPa → %.0f MPa (%.1f%% loss)\n', ...
    mean_E(1), mean_E(end), (mean_E(1)-mean_E(end))/mean_E(1)*100);

%% Example usage of remaining life assessment
current_structure_age = 30; % Example: 30-year-old structure
remaining_life_assessment(t_points, beta_system, current_structure_age);

%% Plot Figures

fprintf('\nGenerating comprehensive plots...\n');

% Create comprehensive reliability and lifetime visualization
figure('Position', [100, 100, 1400, 1000]);

% Subplot 1: Reliability indices over time
subplot(2, 3, 1);
plot(t_points, beta_shear, 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
hold on;
plot(t_points, beta_moment, 'r-s', 'LineWidth', 2, 'MarkerSize', 6);
plot(t_points, beta_deflection, 'g-^', 'LineWidth', 2, 'MarkerSize', 6);
plot(t_points, beta_system, 'k-d', 'LineWidth', 3, 'MarkerSize', 8);
yline(3.8, '--', 'High Safety', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);
yline(3.1, '--', 'Medium Safety', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);
yline(2.3, '--', 'Low Safety', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);
xlabel('Time (years)');
ylabel('Reliability Index β');
title('Reliability Degradation Over Time');
legend('Shear', 'Moment', 'Deflection', 'System', 'Location', 'best');
grid on;

% Subplot 2: Failure probabilities
subplot(2, 3, 2);
semilogy(t_points, Pf_shear, 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
hold on;
semilogy(t_points, Pf_moment, 'r-s', 'LineWidth', 2, 'MarkerSize', 6);
semilogy(t_points, Pf_deflection, 'g-^', 'LineWidth', 2, 'MarkerSize', 6);
semilogy(t_points, Pf_system, 'k-d', 'LineWidth', 3, 'MarkerSize', 8);
xlabel('Time (years)');
ylabel('Failure Probability');
title('Failure Probability Evolution');
legend('Shear', 'Moment', 'Deflection', 'System', 'Location', 'best');
grid on;

% Subplot 3: Material degradation
subplot(2, 3, 3);
yyaxis left;
plot(t_points, mean_tw, 'b-o', 'LineWidth', 2);
hold on;
plot(t_points, mean_tf, 'r-s', 'LineWidth', 2);
ylabel('Thickness (mm)');
yyaxis right;
plot(t_points, mean_fy/1000, 'g-^', 'LineWidth', 2);
ylabel('Yield Strength (GPa)');
xlabel('Time (years)');
title('Material Property Degradation');
legend('Web Thickness', 'Flange Thickness', 'Yield Strength', 'Location', 'best');
grid on;

% Subplot 4: Sensitivity analysis
subplot(2, 3, 4);
plot(t_points, sensitivity(1,:), 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
hold on;
plot(t_points, sensitivity(2,:), 'r-s', 'LineWidth', 2, 'MarkerSize', 6);
plot(t_points, sensitivity(3,:), 'g-^', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Time (years)');
ylabel('Sensitivity Coefficient');
title('Parameter Sensitivity Analysis');
legend('Yield Strength', 'Corrosion Rate', 'Live Load', 'Location', 'best');
grid on;

% Subplot 5: Lifetime prediction summary
subplot(2, 3, 5);
bar_data = [lifetime_results.shear_lifetimes, lifetime_results.moment_lifetimes, ...
           lifetime_results.deflection_lifetimes, lifetime_results.system_lifetimes];
bar(bar_data);
set(gca, 'XTickLabel', {'β=3.8', 'β=3.1', 'β=2.3'});
xlabel('Target Reliability Level');
ylabel('Predicted Lifetime (years)');
title('Predicted Lifetimes by Failure Mode');
legend('Shear', 'Moment', 'Deflection', 'System', 'Location', 'best');
grid on;

% Subplot 6: System lifetime with confidence intervals
subplot(2, 3, 6);
errorbar(1:length(target_betas), lifetime_results.system_lifetimes, ...
 lifetime_results.system_lifetimes - lifetime_results.system_ci(:,1), ...
 lifetime_results.system_ci(:,2) - lifetime_results.system_lifetimes, ...
'ko-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'red');
set(gca, 'XTick', 1:length(target_betas), 'XTickLabel', {'β=3.8', 'β=3.1', 'β=2.3'});
xlabel('Target Reliability Level');
ylabel('Predicted System Lifetime (years)');
title('System Lifetime with 95% Confidence Intervals');
grid on;
sgtitle('Comprehensive Structural Reliability and Lifetime Analysis', 'FontSize', 16, 'FontWeight', 'bold');

 % Failure modeds correlation matrix
figure;
heatmap({'Shear', 'Moment', 'displacement'}, {'Shear', 'Moment', 'displacement'}, rho_failure_modes,...
            'Colormap', jet, 'ColorLimits', [-1, 1], 'Title', 'Failure modes correlation matrix');
    

figure;
[S, AX, BigAx] = plotmatrix(G);
title(BigAx, 'Safety Redundancy Matrix');
labels = {'Shear', 'Moment', 'displacement'};
n = size(G, 2); 
for i = 1:n
     ylabel(AX(i,1), labels{i}, 'FontWeight', 'bold');
     xlabel(AX(n,i), labels{i}, 'FontWeight', 'bold');
    for j = 1:n
        if j > 1
            set(AX(i,j), 'YTickLabel', []);
        end
        if i < n
            set(AX(i,j), 'XTickLabel', []);
        end
    end
end
set(AX, 'FontSize', 10, 'Box', 'on');
set(gcf, 'Position', [100 100 800 600]); 


%% Med-Mrd-time

figure('Position', [100, 100, 1000, 800], 'Name', 'Moment Resistance-Stress Distribution Evolution');
hold on;

mu_Med = zeros(1, n_time);
sd_Med = zeros(1, n_time);
mu_Mrd = zeros(1, n_time);
sd_Mrd = zeros(1, n_time);

for i = 1:n_time
    mu_Med(i) = mean(Med_all(:, i));
    sd_Med(i) = std(Med_all(:, i));
    mu_Mrd(i) = mean(Mrd_all(:, i));
    sd_Mrd(i) = std(Mrd_all(:, i));
end
% Set colour map for time representation
colors = parula(n_time);  % Use parula colour map

% Loop through each time point to plot PDF curves
for i = 1:n_time
    % Calculate lognormal distribution parameters
    sigma_ln_M_Ed = sqrt(log(1 + (sd_Med(i)/mu_Med(i))^2));
    mu_ln_M_Ed = log(mu_Med(i)) - 0.5*sigma_ln_M_Ed^2;
    sigma_ln_M_Rd = sqrt(log(1 + (sd_Mrd(i)/mu_Mrd(i))^2));
    mu_ln_M_Rd = log(mu_Mrd(i)) - 0.5*sigma_ln_M_Rd^2;
    
    % Create abscissa range
    M_range = linspace(0, max([Med_all(:); Mrd_all(:)]), 1000);
    % Calculate PDFs
    pdf_Med = lognpdf(M_range, mu_ln_M_Ed, sigma_ln_M_Ed);
    pdf_Mrd = lognpdf(M_range, mu_ln_M_Rd, sigma_ln_M_Rd);
    % Plot curves
    plot(M_range, pdf_Med, 'Color', colors(i, :), 'LineWidth', 1.5, 'LineStyle', '-');
    plot(M_range, pdf_Mrd, 'Color', colors(i, :), 'LineWidth', 1.5, 'LineStyle', '--');
end

% legend
legend_str = cell(1, 2*n_time);
for i = 1:n_time
    legend_str{2*i-1} = sprintf('M_{Ed} (%d years)', t_points(i));
    legend_str{2*i} = sprintf('M_{Rd} (%d years)', t_points(i));
end
legend(legend_str, 'Location', 'northeastoutside', 'FontSize', 8);

% labels and title
xlabel('Bending Moment (N·mm)', 'FontSize', 12);
ylabel('Probability Density', 'FontSize', 12);
title('Evolution of Moment Reliability Distribution', 'FontSize', 14);

% failure probability annotation
annotation('textbox', [0.51, 0.71, 0.2, 0.18], 'String', ...
    sprintf('Distribution Trends:\n• Resistance shifts left (degradation)\n• Stress shifts right (loading increase)\n• Overlap area = failure probability'), ...
    'BackgroundColor', 'white', 'FontSize', 10);

line_x = linspace(min(M_range), max(M_range), 100);
line_y = zeros(size(line_x));
plot(line_x, line_y, 'k--', 'LineWidth', 1.5);
text(mean(M_range), -0.0001, 'Safe boundary', 'FontSize', 10, 'HorizontalAlignment', 'center');

% Add colour bar for time representation
colormap(parula);
c = colorbar;
c.Label.String = 'Time (years)';
c.Ticks = linspace(0, 1, n_time);
c.TickLabels = arrayfun(@num2str, t_points, 'UniformOutput', false);

grid on;
set(gca, 'FontSize', 11, 'Box', 'on');
hold off;
%% Ved-Vrd-Time

figure('Position', [100, 100, 1000, 800], 'Name', 'Shear Resistance-Stress Distribution Evolution');
hold on;

mu_Ved = zeros(1, n_time);
sd_Ved = zeros(1, n_time);
mu_Vrd = zeros(1, n_time);
sd_Vrd = zeros(1, n_time);

for i = 1:n_time
    mu_Ved(i) = mean(Ved_all(:, i));
    sd_Ved(i) = std(Ved_all(:, i));
    mu_Vrd(i) = mean(Vrd_all(:, i));
    sd_Vrd(i) = std(Vrd_all(:, i));
end

colors = parula(n_time); 

for i = 1:n_time
    sigma_ln_V_Ed = sqrt(log(1 + (sd_Ved(i)/mu_Ved(i))^2));
    mu_ln_V_Ed = log(mu_Ved(i)) - 0.5*sigma_ln_V_Ed^2;
    sigma_ln_V_Rd = sqrt(log(1 + (sd_Vrd(i)/mu_Vrd(i))^2));
    mu_ln_V_Rd = log(mu_Vrd(i)) - 0.5*sigma_ln_V_Rd^2;
    
   
    V_range = linspace(0, max([Ved_all(:); Vrd_all(:)]), 1000);
    pdf_Ved = lognpdf(V_range, mu_ln_V_Ed, sigma_ln_V_Ed);
    pdf_Vrd = lognpdf(V_range, mu_ln_V_Rd, sigma_ln_V_Rd);
    plot(V_range, pdf_Ved, 'Color', colors(i, :), 'LineWidth', 1.5, 'LineStyle', '-');
    plot(V_range, pdf_Vrd, 'Color', colors(i, :), 'LineWidth', 1.5, 'LineStyle', '--');
end

legend_str = cell(1, 2*n_time);
for i = 1:n_time
    legend_str{2*i-1} = sprintf('V_{Ed} (%d years)', t_points(i));
    legend_str{2*i} = sprintf('V_{Rd} (%d years)', t_points(i));
end
legend(legend_str, 'Location', 'northeastoutside', 'FontSize', 8);

xlabel('Shear Force (N)', 'FontSize', 12);
ylabel('Probability Density', 'FontSize', 12);
title('Evolution of Shear Reliability Distribution', 'FontSize', 14);

% Add failure probability annotation
annotation('textbox', [0.51, 0.71, 0.2, 0.18], 'String', ...
    sprintf('Distribution Trends:\n• Shear resistance shifts left (web corrosion)\n• Shear stress shifts right (loading increase)\n• Overlap area = failure probability'), ...
    'BackgroundColor', 'white', 'FontSize', 10);

colormap(parula);
c = colorbar;
c.Label.String = 'Time (years)';
c.Ticks = linspace(0, 1, n_time);
c.TickLabels = arrayfun(@num2str, t_points, 'UniformOutput', false);

grid on;
set(gca, 'FontSize', 11, 'Box', 'on');
hold off;

%% Deflection - Deflection limit - time
figure('Position', [100, 100, 1000, 800], 'Name', 'Deflection Serviceability Distribution Evolution');
hold on;

mu_delta = zeros(1, n_time);
sd_delta = zeros(1, n_time);
mu_delta_limit = zeros(1, n_time);
sd_delta_limit = zeros(1, n_time);

for i = 1:n_time
    mu_delta(i) = mean(delta_all(:, i));
    sd_delta(i) = std(delta_all(:, i));
    mu_delta_limit(i) = mean(delta_limit_all(:, i));
    sd_delta_limit(i) = std(delta_limit_all(:, i));
end

colors = parula(n_time); 

for i = 1:n_time
    sigma_ln_delta = sqrt(log(1 + (sd_delta(i)/mu_delta(i))^2));
    mu_ln_delta = log(mu_delta(i)) - 0.5*sigma_ln_delta^2;
    sigma_ln_delta_limit = sqrt(log(1 + (sd_delta_limit(i)/mu_delta_limit(i))^2));
    mu_ln_delta_limit = log(mu_delta_limit(i)) - 0.5*sigma_ln_delta_limit^2;
    
    delta_range = linspace(0, max([delta_all(:); delta_limit_all(:)]), 1000);

    pdf_delta = lognpdf(delta_range, mu_ln_delta, sigma_ln_delta);
    pdf_delta_limit = lognpdf(delta_range, mu_ln_delta_limit, sigma_ln_delta_limit);

    plot(delta_range, pdf_delta, 'Color', colors(i, :), 'LineWidth', 1.5, 'LineStyle', '-');
    plot(delta_range, pdf_delta_limit, 'Color', colors(i, :), 'LineWidth', 1.5, 'LineStyle', '--');
end

legend_str = cell(1, 2*n_time);
for i = 1:n_time
    legend_str{2*i-1} = sprintf('Actual Deflection (%d years)', t_points(i));
    legend_str{2*i} = sprintf('Deflection Limit (%d years)', t_points(i));
end
legend(legend_str, 'Location', 'northeastoutside', 'FontSize', 8);

xlabel('Deflection (mm)', 'FontSize', 12);
ylabel('Probability Density', 'FontSize', 12);
title('Evolution of Deflection Serviceability Distribution', 'FontSize', 14);

annotation('textbox', [0.11, 0.71, 0.2, 0.18], 'String', ...
    sprintf('Distribution Trends:\n• Actual deflection shifts right (stiffness degradation)\n• Deflection limit remains stable\n• Overlap area = serviceability failure probability'), ...
    'BackgroundColor', 'white', 'FontSize', 10);

colormap(parula);
c = colorbar;
c.Label.String =  'Time (years)';
c.Ticks = linspace(0, 1, n_time);
c.TickLabels = arrayfun(@num2str, t_points, 'UniformOutput', false);

grid on;
set(gca, 'FontSize', 11, 'Box', 'on');
hold off;


%% 1. Bending Limit State (Med-Mrd)
num_plot_samples = min(500, nSamples); % Plot 500 samples
selected_samples = randperm(nSamples, num_plot_samples);
figure('Position', [100, 100, 800, 600], 'Name', 'Bending Limit State');
hold on;

% Plot mean evolution curve
plot3(t_points, mu_Med, mu_Mrd, 'r-', 'LineWidth', 1);

% Plot sample points
for i = 1:n_time
    scatter3(...
        t_points(i) * ones(length(selected_samples), 1), ... % X: Time
        Med_all(selected_samples, i), ... % Y: Med
        Mrd_all(selected_samples, i), ... % Z: Mrd
        10, ... % Point size
        'b', ... % Colour
        'filled', ... % Solid points
        'MarkerFaceAlpha', 0.3); % Transparency
end

% Add safety boundary surface
[t_grid, M_grid] = meshgrid(linspace(0, design_life, 10), linspace(min(Mrd_all(:)), max(Mrd_all(:)), 10));
Med_bound = M_grid; % M_r = M_d boundary
surf(t_grid, M_grid, Med_bound, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);

% Labelling and formatting
xlabel('Time (years)');
ylabel('Load Effect M_{Ed} (N·mm)');
zlabel('Resistance M_{Rd} (N·mm)');
title('Evolution of Bending Limit State');
legend('Mean evolution', 'Sample points', 'Safety boundary (M_r = M_d)', 'Location', 'best');
view(45, 30);
grid on;
set(gca, 'FontSize', 12, 'Box', 'on');
hold off;

%% 2. Shear Limit State (Ved-Vrd)
figure('Position', [100, 100, 800, 600], 'Name', 'Shear Limit State');
hold on;

plot3(t_points, mu_Ved, mu_Vrd, 'g-', 'LineWidth', 1);

for i = 1:n_time
    scatter3(...
        t_points(i) * ones(length(selected_samples), 1), ... % X: Time
        Ved_all(selected_samples, i), ... % Y: Ved
        Vrd_all(selected_samples, i), ... % Z: Vrd
        10, ...
        'm', ... % 
        'filled', ...
        'MarkerFaceAlpha', 0.3);
end

[t_grid, V_grid] = meshgrid(linspace(0, design_life, 10), linspace(min(Vrd_all(:)), max(Vrd_all(:)), 10));
Ved_bound = V_grid; % V_r = V_d 
surf(t_grid, V_grid, Ved_bound, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);

xlabel('Time (years)');
ylabel('Shear Force V_{Ed} (N)');
zlabel('Shear Resistance V_{Rd} (N)');
title('Evolution of Shear Limit State');
legend('Mean evolution', 'Sample points', 'Safety boundary (V_r = V_d)', 'Location', 'best');
view(45, 30);
grid on;
set(gca, 'FontSize', 12, 'Box', 'on');
hold off;

%% 3. Deflection Serviceability Limit State (Delta-Delta Limit)
figure('Position', [100, 100, 800, 600], 'Name', 'Deflection Serviceability Limit State');
hold on;

plot3(t_points, mu_delta, mu_delta_limit, 'b-', 'LineWidth', 1);

for i = 1:n_time
    scatter3(...
        t_points(i) * ones(length(selected_samples), 1), ... % X: Time
        delta_all(selected_samples, i), ... % Y: Actual deflection
        delta_limit_all(selected_samples, i), ... % Z: Deflection limit
        10, ...
        'c', ... % Cyan colour
        'filled', ...
        'MarkerFaceAlpha', 0.3);
end

[t_grid, D_grid] = meshgrid(linspace(0, design_life, 10), linspace(min(delta_all(:)), max(delta_all(:)), 10));
delta_bound = D_grid; % Δ = Δ_limit 
surf(t_grid, D_grid, delta_bound, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);

% Labelling and formatting
xlabel('Time (years)');
ylabel('Actual Deflection Δ (mm)');
zlabel('Deflection Limit Δ_{limit} (mm)');
title('Evolution of Deflection Serviceability');
legend('Mean evolution', 'Sample points', 'Serviceability boundary (Δ = Δ_{limit})', 'Location', 'best');
view(45, 30);
grid on;
set(gca, 'FontSize', 12, 'Box', 'on');
hold off;

%% Corrosion Calculation 

function D = calc_corrosion_ISO9224(t, metal_type, samples_kcor1, use_conservative)
% Calculate metal loss under long-term corrosion according to ISO 9224:2012
% Inputs:
%   t: Exposure time in years (scalar or vector)
%   metal_type: Metal type ('steel', 'zinc', 'copper', 'aluminium')
%   r_corr: First-year corrosion rate (μm/year or g/m²/year)
%   use_conservative: true for upper-limit estimation (B2 values), false for average (B1 values)
% Output:
%   D: Total corrosion loss (same units as r_corr)

% Define time exponent b values from Table 2
b_values = struct(...
    'steel',     struct('B1', 0.523, 'B2', 0.575), ...
    'zinc',      struct('B1', 0.813, 'B2', 0.873), ...
    'copper',    struct('B1', 0.667, 'B2', 0.726), ...
    'aluminium', struct('B1', 0.728, 'B2', 0.807));

% Validate metal type input
valid_metals = {'steel', 'zinc', 'copper', 'aluminium'};
if ~ismember(metal_type, valid_metals)
    error('Invalid metal type. Valid options: steel, zinc, copper, aluminium');
end

% Select appropriate b-value (B1 for average, B2 for conservative)
b_data = b_values.(metal_type);
if use_conservative
    b = b_data.B2;  % Upper-limit estimation
else
    b = b_data.B1;  % Average estimation
end

% Initialize output array
D = zeros(size(t));

% Calculate corrosion loss for each time point
for i = 1:numel(t)
    if t(i) <= 20
        % Equation (1): Power-law model for t ≤ 20 years
        % D = r_corr * t^b
        D(i) = samples_kcor1 * (t(i)^b);
        
    else
        % Equation (3): Linear model for t > 20 years (conservative)
        % D = r_corr * [20^b + b * 20^(b-1) * (t - 20)]
        D_20 = kcor1 * (20^b);        % Corrosion loss at 20 years
        k = b * (20^(b-1));            % Linear corrosion rate coefficient
        D(i) = D_20 + samples_kcor1 * k * (t(i) - 20);  % Total loss
    end
end
end

%% Lifetime Prediction Function (Integrated) - Fixed for indexing issues
function [life, ci] = predict_lifetime(t_data, beta_data, target_beta)
    % Ensure inputs are column vectors and same length
    t_data = t_data(:);
    beta_data = beta_data(:);
    
    if length(t_data) ~= length(beta_data)
        fprintf('Error: t_data and beta_data must have same length\n');
        life = NaN;
        ci = [NaN, NaN];
        return;
    end
    
    % Remove any NaN or Inf values
    valid_idx = isfinite(beta_data) & isfinite(t_data);
    
    if sum(valid_idx) < 2
        fprintf('Warning: Insufficient valid data points for interpolation\n');
        life = NaN;
        ci = [NaN, NaN];
        return;
    end
    
    t_data = t_data(valid_idx);
    beta_data = beta_data(valid_idx);
    
    % Handle duplicate beta values by averaging corresponding time values
    [beta_unique, ~, idx] = unique(beta_data, 'stable');
    if length(beta_unique) < length(beta_data)
        fprintf('Warning: Found duplicate beta values, averaging time values\n');
        t_unique = zeros(size(beta_unique));
        for i = 1:length(beta_unique)
            t_unique(i) = mean(t_data(idx == i));
        end
        beta_data = beta_unique;
        t_data = t_unique;
    end
    
    % Final check for sufficient data
    if length(beta_data) < 2
        fprintf('Warning: Insufficient unique data points after processing\n');
        life = NaN;
        ci = [NaN, NaN];
        return;
    end
    
    % Check if target_beta is within the data range
    beta_min = min(beta_data);
    beta_max = max(beta_data);
    
    if target_beta > beta_max
        fprintf('Warning: Target β=%.2f is above maximum β=%.2f in data\n', target_beta, beta_max);
        life = 0;
        ci = [0, 0];
        return;
    end
    
    if target_beta < beta_min
        fprintf('Warning: Target β=%.2f is below minimum β=%.2f in data\n', target_beta, beta_min);
        % Extrapolate to find lifetime
        if length(t_data) >= 2
            % Use linear extrapolation
            p = polyfit(t_data, beta_data, 1);
            if abs(p(1)) > 1e-10 % Check if slope is not zero
                life = (target_beta - p(2)) / p(1);
                life = max(max(t_data), life); % Ensure it's beyond current data
            else
                life = 100; % Default long lifetime
            end
        else
            life = 100;
        end
        ci = [life, inf];
        return;
    end
    
    try
        % Sort data by time for proper interpolation (time should be monotonic)
        [t_sorted, sort_idx] = sort(t_data);
        beta_sorted = beta_data(sort_idx);
        
        % Check if beta decreases monotonically with time (expected behavior)
        if ~issorted(beta_sorted, 'descend')
            fprintf('Warning: Beta values are not monotonically decreasing with time\n');
        end
        
        % Use interpolation with sorted data
        life = interp1(beta_sorted, t_sorted, target_beta, 'linear', 'extrap');
        
        % If interpolation gives unreasonable result, try polynomial fitting
        if isnan(life) || life < 0
            if length(t_data) >= 3
                try
                    % Fit polynomial to time-ordered data
                    p = polyfit(t_sorted, beta_sorted, min(2, length(t_sorted)-1));
                    
                    if length(p) == 3 % Quadratic
                        coeffs = [p(1), p(2), p(3) - target_beta];
                        roots_t = roots(coeffs);
                        real_roots = real(roots_t(abs(imag(roots_t)) < 1e-10));
                        valid_roots = real_roots(real_roots >= min(t_sorted) & real_roots <= max(t_sorted)*2);
                        
                        if ~isempty(valid_roots)
                            life = min(valid_roots);
                        else
                            % Use linear extrapolation as last resort
                            p_linear = polyfit(t_sorted, beta_sorted, 1);
                            life = (target_beta - p_linear(2)) / p_linear(1);
                        end
                    else % Linear
                        life = (target_beta - p(2)) / p(1);
                    end
                catch
                    % Final fallback: simple linear extrapolation
                    if length(t_sorted) >= 2
                        slope = (beta_sorted(end) - beta_sorted(1)) / (t_sorted(end) - t_sorted(1));
                        intercept = beta_sorted(1) - slope * t_sorted(1);
                        life = (target_beta - intercept) / slope;
                    else
                        life = t_sorted(1);
                    end
                end
            else
                % Linear interpolation/extrapolation for 2 points
                slope = (beta_sorted(2) - beta_sorted(1)) / (t_sorted(2) - t_sorted(1));
                intercept = beta_sorted(1) - slope * t_sorted(1);
                life = (target_beta - intercept) / slope;
            end
        end
        
        % Calculate confidence interval
        if length(t_data) >= 3
            try
                residuals = beta_sorted - polyval(polyfit(t_sorted, beta_sorted, min(2, length(t_sorted)-1)), t_sorted);
                se = std(residuals);
            catch
                se = std(beta_data) * 0.1;
            end
        else
            se = std(beta_data) * 0.1; % Conservative estimate
        end
        
        % Convert to time domain uncertainty
        if length(t_data) >= 2
            dbeta_dt = abs((beta_sorted(end) - beta_sorted(1)) / (t_sorted(end) - t_sorted(1)));
            if dbeta_dt > 1e-6
                dt_se = se / dbeta_dt;
            else
                dt_se = 5; % Default uncertainty
            end
        else
            dt_se = 5;
        end
        
        ci_low = max(0, life - 1.96*dt_se);
        ci_high = life + 1.96*dt_se;
        ci = [ci_low, ci_high];
        
    catch ME
        fprintf('Error in lifetime prediction: %s\n', ME.message);
        fprintf('Using robust fallback method\n');
        
        % Most robust fallback: find closest points and interpolate
        [~, closest_idx] = min(abs(beta_data - target_beta));
        life = t_data(closest_idx);
        ci = [max(0, life-5), life+5];
    end
    
    % Ensure reasonable bounds
    if isnan(life) || life < 0
        life = 0;
    end
    if any(isnan(ci)) || ci(1) < 0
        ci = [max(0, life-5), life+5];
    end
    
    % Only create plot if we have reasonable results
    if isfinite(life) && life >= 0 && length(t_data) >= 2
        try
            figure('Position', [200, 200, 800, 600]);
            
            % Plot data points
            scatter(t_data, beta_data, 100, 'blue', 'filled', 'DisplayName', 'Computed β');
            hold on;
            
            % Create fitted curve for visualization
            t_fine = linspace(min(t_data), max(max(t_data), life*1.1), 100);
            if length(t_data) >= 3
                try
                    [t_plot, sort_idx] = sort(t_data);
                    beta_plot = beta_data(sort_idx);
                    beta_fine = interp1(t_plot, beta_plot, t_fine, 'pchip', 'extrap');
                    plot(t_fine, beta_fine, 'b-', 'LineWidth', 2, 'DisplayName', 'Fitted Curve');
                catch
                    % Use linear interpolation
                    beta_fine = interp1(t_data, beta_data, t_fine, 'linear', 'extrap');
                    plot(t_fine, beta_fine, 'b-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
                end
            else
                beta_fine = interp1(t_data, beta_data, t_fine, 'linear', 'extrap');
                plot(t_fine, beta_fine, 'b-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
            end
            
            % Add target line and predicted life
            yline(target_beta, 'r--', sprintf('Target β = %.2f', target_beta), ...
                  'LineWidth', 2, 'LabelHorizontalAlignment', 'left');
            
            if isfinite(life) && life > 0
                xline(life, 'g--', sprintf('Predicted Life = %.1f years', life), ...
                      'LineWidth', 2, 'LabelVerticalAlignment', 'bottom');
                
                % Highlight intersection point
                scatter(life, target_beta, 200, 'red', 'filled', 'DisplayName', 'Predicted Lifetime');
                
                % Add confidence interval
                if all(isfinite(ci)) && ci(2) > ci(1)
                    fill([ci(1), ci(2), ci(2), ci(1)], [target_beta-0.05, target_beta-0.05, target_beta+0.05, target_beta+0.05], ...
                         'green', 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', '95% CI');
                end
            end
            
            xlabel('Time (years)', 'FontSize', 12);
            ylabel('Reliability Index β', 'FontSize', 12);
            title(sprintf('Lifetime Prediction: %.1f years (95%% CI: %.1f-%.1f years)', ...
                  life, ci(1), ci(2)), 'FontSize', 14, 'FontWeight', 'bold');
            legend('Location', 'best');
            grid on;
            
            % Set axis limits
            xlim([0, max(max(t_data)*1.1, life*1.2)]);
            ylim([min(min(beta_data)-0.2, target_beta-0.3), max(beta_data)+0.2]);
            
        catch plot_error
            fprintf('Warning: Could not create lifetime plot: %s\n', plot_error.message);
        end
    end
end


%% Ditlevsen Boundary Analysis Function
function [Pf_lower, Pf_upper] = ditlevsen_bounds(betas, rho_failure_modes)
    n = length(betas);
    Pf_individual = normcdf(-betas);
    
    P_joint = zeros(n);
    for i = 1:n
        for j = (i+1):n
            P_joint(i,j) = mvncdf([-betas(i), -betas(j)], [0,0], [1, rho_failure_modes(i,j); rho_failure_modes(i,j), 1]);
        end
    end
    
    Pf_upper = Pf_individual(1);
    Pf_lower = Pf_individual(1);
    
    for k = 2:n
        sum_term = 0;
        for j = 1:(k-1)
            sum_term = sum_term + max(0, P_joint(j,k) - Pf_upper);
        end
        Pf_upper = Pf_upper + Pf_individual(k) - sum_term;
        
        sum_term = 0;
        for j = 1:(k-1)
            sum_term = sum_term + P_joint(j,k);
        end
        Pf_lower = Pf_lower + max(0, Pf_individual(k) - sum_term);
    end
    
    Pf_lower = max(0, min(1, Pf_lower));
    Pf_upper = max(0, min(1, Pf_upper));
end


%% Additional Analysis: Remaining Life Assessment
function remaining_life_assessment(t_points, beta_system, current_age)
    fprintf('\n--- REMAINING LIFE ASSESSMENT ---\n');
    fprintf('Current structure age: %.1f years\n', current_age);
    
    if current_age > max(t_points)
        fprintf('Warning: Current age exceeds analysis range\n');
        return;
    end
    
    % Interpolate current reliability index
    current_beta = interp1(t_points, beta_system, current_age);
    fprintf('Current system reliability index: %.3f\n', current_beta);
    
    % Predict remaining life for different targets
    target_betas = [3.8, 3.1, 2.3, 1.5];
    target_names = {'High Safety (β=3.8)','Medium Safety (β=3.1)', 'Low Safety (β=2.3)', 'Critical (β=1.5)'};
    
    for i = 1:length(target_betas)
        if current_beta > target_betas(i)
            % Predict when target will be reached with linear extrapolation
            future_time = interp1(beta_system, t_points, target_betas(i), 'linear', 'extrap');
            remaining_years = future_time - current_age;
            fprintf('%s: %.1f years remaining\n', target_names{i}, remaining_years);
        else
            fprintf('%s: Target already exceeded\n', target_names{i});
        end
    end
end

%% 
