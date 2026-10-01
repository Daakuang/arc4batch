function [performance,fig,SimData]=arc_controller(sys,par,dk)


% ARC_PID Advanced Regulatory Control for Batch Polymerization
% Implements the "Inverted ARC" / "Valve Position Control" (VPC) structure.
%
% Structure:
% 1. Pressure Control Loop (P -> FA): Manipulates Feed A to control Pressure.
% 2. Quality/Temperature Control Loop (Cascade):
%    - Master (Outer): T_reactor -> T_jacket_sp
%    - Slave (Inner): T_jacket -> Cooling Duty (F_cw_req / alpha)
% 3. Economic/VPC Loop (VPC -> FB): Manipulates Feed B to keep Cooling Duty near maximum (thermal ceiling).
%
% History:
% Refactored to align with paper/V1.tex terminology.
% Update: Fixed Time_hours error. Restored Hybrid Control Logic.

%% 1. Initialization and Parameters

if isfield(par, 'arc_variant')
    arc_variant = char(par.arc_variant);
else
    arc_variant = 'full_arc';
end

valid_variants = {'full_arc', 'fixed_psp', 'recipe_fb', 'recipe_only'};
assert(any(strcmp(arc_variant, valid_variants)), ...
    'Unsupported arc_variant "%s".', arc_variant);

use_fb_feedback = any(strcmp(arc_variant, {'full_arc', 'fixed_psp'}));
use_psp_econ = any(strcmp(arc_variant, {'full_arc', 'recipe_fb'}));

% State and inputs
xk = par.x0;
uk = par.u0;
% dk comes from arguments

% --- Setpoints and Constraints ---
Trsp    = 273 + 78;   % [K] Reactor temperature set point
Psp     = 1500000;    % [Pa] Reactor gas pressure set point
Mmax_A  = 3250;       % [kg] Maximum feed mass of A
Mlim_B  = 0.0020;     % [kg] Feeding APS mass limit for changing the B feed rate
Pact_Tc = 1000000;    % [Pa] Pressure threshold for activation

% --- Controller Parameters ---

% 1. Pressure Controller (PI) - Manipulates FA
Kp_P = 4;
Ti_P = 20;
Td_P = 0.0;
ui_P(1) = 0;   % Integrator state

% 2. Temperature Master Controller (Outer Loop) (PI) - T -> Tjsp
Kp_T_master = 35;
Ti_T_master = 300;
Td_T_master = 0;
Tt_T_master = 60;   % Back-calculation tracking time constant
ui_T_master(1) = 0; ui_T_master(2) = 0;

% 3. Temperature Slave Controller (Inner Loop) (PD) - Tj -> F_heating/cooling
Kp_T_slave = 2;
Kd_T_slave = 0;
k_sw_slave = 1;  % Smooth switching slope in slave nonlinear mapping
% Scaling/Gain for feedforward (if any)
K_FF_P = 5; % Feedforward from FA to Tjsp (compensate pressure control effect)

% 4. Economic / VPC Controller (PI) - CoolingDemand -> FB
% "Inverted ARC": Adjust FB to keep 'alpha' (cooling demand) high (e.g., -1).
if isfield(par, 'K_VPC'), K_VPC = par.K_VPC; else K_VPC = 8; end
if isfield(par, 'Ti_VPC'), Ti_VPC = par.Ti_VPC; else Ti_VPC = 4e1; end
if isfield(par, 'VPC_Setpoint'), VPC_Setpoint = par.VPC_Setpoint; else VPC_Setpoint = -0.95; end
if isfield(par, 'ui_VPC_max'), ui_VPC_max = par.ui_VPC_max; else ui_VPC_max = 1.0; end
if isfield(par, 'ui_VPC_min'), ui_VPC_min = par.ui_VPC_min; else ui_VPC_min = -1.5; end
ui_VPC(1) = 0; % Integral state
val_FB_clip_min = -40;
val_FB_clip_max = 10;
val_FB_to_FB = 1e-7;

% Psp Safety Logic Parameters
Psp_Nominal = 1500000; % Base Pressure Setpoint
if isfield(par, 'K_Psp_Relax'), K_Psp_Relax = par.K_Psp_Relax; else K_Psp_Relax = 1e4; end
if isfield(par, 'VPC_Psp_Overload_Gain'), VPC_Psp_Overload_Gain = par.VPC_Psp_Overload_Gain; else VPC_Psp_Overload_Gain = 0.0; end

% Shared-controller mapping: one VPC output drives both FB and Psp relaxation.

% --- Actuator Model Parameters ---
% First-order lag time constants (seconds)
tau_FA = 5.0;       % Time constant for FA actuator
tau_FB = 5.0;       % Time constant for FB actuator
tau_Alpha = 3.0;    % Time constant for Alpha actuator

% Rate limits (max change per second)
du_max_FA = par.ubu(1) * 0.05;      % Max 20% of range per second
du_max_FB = par.ubu(2) * 0.05;      % Max 20% of range per second
du_max_Alpha = par.ubu(3) * 0.05;   % Max 50% of range per second

% --- Simulation Setup ---
if isfield(par, 'max_steps')
    sim_steps = par.max_steps;
else
    sim_steps = 30000; % Increased to cover 8h+ if needed
end
t_sim = 1:sim_steps;
Ramp_Steps = 3300 / par.tf; % Duration of S-curve startup

% Dynamic-limit calibration for outer loop from actual inner-loop mapping.
inner_cmd_min = max(par.lbu(3), -1);  % actual minimum command after clamp
inner_cmd_max = par.ubu(3);
alpha_raw_for_inner_min = invert_slave_alpha_cmd(inner_cmd_min, par.lbu(3), par.ubu(3), k_sw_slave);
alpha_raw_for_inner_max = invert_slave_alpha_cmd(inner_cmd_max, par.lbu(3), par.ubu(3), k_sw_slave);

% Preallocation
[C_R, C_A, C_B, C_C, C_E, Tr, NgA, P, Tjacket, RA, Mn, Mw, PDI, MAin, MBin, Vf, kt, X_t, k70] = deal(zeros(1, sim_steps));
Fgl = zeros(1, sim_steps);
% Control signals
[FA, FB, u_cool_demand, Tjsp_trace, Psp_trace, val_FB_trace] = deal(zeros(1, sim_steps));
[e_P, e_T_master, e_T_slave] = deal(zeros(1, sim_steps));
[ui_P_trace, ui_T_master_trace, ui_VPC_trace] = deal(zeros(1, sim_steps));
[switch_P_active_trace] = deal(false(1, sim_steps));
[FB_recipe_cmd_trace, FB_feedback_cmd_trace] = deal(zeros(1, sim_steps));
[Psp_econ_active_trace, FB_feedback_active_trace, FB_channel_active_trace] = deal(false(1, sim_steps));
[cooling_slack_trace, cooling_distance_trace] = deal(zeros(1, sim_steps));
[Tjsp_upper_dyn_trace, Tjsp_lower_dyn_trace] = deal(nan(1, sim_steps));
solve_times = zeros(1, sim_steps);
u_cmd_actual_trace = zeros(3,sim_steps);

% Switch flags
switch_T_active = 0;
switch_P_active = 0;
switch_emergency_shutdown = 0;

% Noise generation
rng(par.seed, 'twister');
mNoise = par.paired_noise;

%% 2. Main Simulation Loop
for i = 1:sim_steps

    % --- System Simulation Step ---
    yk = full(sys.m(xk, uk, dk));

    % --- Start CPU Timer (Control Step) ---
    t_start_cpu = tic;

    % Unpack outputs (State observation)
    C_R(i) = yk(1);    C_A(i) = yk(2);    C_B(i) = yk(3);    C_C(i) = yk(4);
    C_E(i) = yk(5);    Tr(i) = yk(6);     NgA(i) = yk(7);    P(i) = yk(8);
    Tjacket(i) = yk(9); RA(i) = yk(10);   Mn(i) = yk(11);    Mw(i) = yk(12);
    PDI(i) = yk(13);   MAin(i) = yk(14);  MBin(i) = yk(15);
    Vf(i) = yk(16);    kt(i) = yk(17);    X_t(i) = yk(18);
    if length(yk) >= 21
        k70(i) = yk(21);  % k70_effective (Gel Effect modified)
    end
    if length(yk) >= 22
        Fgl(i) = yk(22);
    else
        % Backward-compatible fallback if model output has not been extended yet.
        Fgl(i) = 0.0562 * (P(i) * (1/700) - 1000 * C_A(i));
    end

    % Add Measurement Noise
    noise_scale = 1.0;
    if isfield(par, 'noise_scale'), noise_scale = par.noise_scale; end

    Tr(i)      = Tr(i) + noise_scale * 0.1/3 * mNoise(1,i);
    Tjacket(i) = Tjacket(i) + noise_scale * 0.1/3 * mNoise(2,i);
    P(i)       = P(i) + noise_scale * 5e3/3 * mNoise(3,i);

    % Filter Measurements (First-order lag)
    if i >= 2 && noise_scale > 0
        alpha_filter = 1/(1 + 3/par.tf);
        Tr(i)      = (1 - alpha_filter)*Tr(i-1) + alpha_filter*Tr(i);
        Tjacket(i) = (1 - alpha_filter)*Tjacket(i-1) + alpha_filter*Tjacket(i);
        P(i)       = (1 - alpha_filter)*P(i-1) + alpha_filter*P(i);
    end

    % ================= CONTROL LOGIC =================

    % --- 1. Loop: FA -> P (Pressure Control) ---
    if i < 2
        FA(i) = par.ubu(1);
        ui_P(i) = 0;
    else
        if MAin(i) >= Mmax_A-1e-6
            e_P(i) = 0;
            ui_P(i) = max(min(ui_P(i-1), par.ubu(1)), par.lbu(1));
            FA(i) = 0;
        else
            e_P(i) = (Psp - P(i)) / 1e5;
            if switch_P_active
                ui_P(i) = max(min(ui_P(i-1) + Kp_P * par.tf / Ti_P * e_P(i), par.ubu(1)), par.lbu(1));
            else
                ui_P(i) = ui_P(i-1);
            end
            term_P = Kp_P * e_P(i);
            term_D = Kp_P * Td_P * (e_P(i) - e_P(i-1)) / par.tf;
            FA(i) = term_P + term_D + ui_P(i);
        end
    end
    FA(i) = max(par.lbu(1), min(FA(i), par.ubu(1)));
    uk(1) = FA(i);
    ui_P_trace(i) = ui_P(i);

    % --- 2. Loop: Cascade Temperature Control (Tr -> Tjacket -> u_cool) ---
    if i >= 2
        % -- Master Loop (Tr -> Tjsp) --

        % Trajectory Generation (S-curve)
        Tr_start = 273 + 25; % Assumed ambient/start temp. ideally use Tr(1) but Par.x0 is better?
        % Tr(1) is available.
        Tr_start_val = Tr(1)+5; % Capture start temp

        Trsp_final = 273 + 78;
        % Ramp_Steps defined in Init section

        if i <= Ramp_Steps
            tau = (i-1) / Ramp_Steps;
            % Minimum-jerk trajectory (5th order polynomial)
            % beta = 10*tau^3 - 15*tau^4 + 6*tau^5
            beta = tau^3 * (10 - 15*tau + 6*tau^2);
            Trsp_curr = Tr_start_val + (Trsp_final - Tr_start_val) * beta;
        else
            Trsp_curr = Trsp_final;
        end

        % Trace key variables
        % Note: Trsp variable was constant before, now we use local Trsp_curr

        e_T_master(i) = Trsp_curr - Tr(i);

        % Hard safety bounds for Tjsp
        Tjsp_max_hard = 400; % [K]
        Tjsp_min_hard = 280; % [K]

        % Integrator limits
        ui_T_master_max = 60;   % Max integrator output (Tjsp offset)
        ui_T_master_min = -60;  % Min integrator output

        term_P_Master = Kp_T_master * e_T_master(i);
        term_D_Master = Kp_T_master * Td_T_master * (e_T_master(i) - e_T_master(i-1)) / par.tf;

        % Feedforward (Pre-calculate)
        FF_term = 0;
        if P(i) >= Psp || switch_P_active == 1
            FF_term = -K_FF_P * (FA(i) / 0.694)*0;
            switch_P_active = 1;
        end

        % Continuous Control (Trajectory Tracking)
        switch_T_active = 1; % Always active

        % Dynamic outer-loop bounds from actual inner-loop mapping (including
        % slave nonlinear scaling), converted back to Tjsp bounds.
        kp_in_safe = Kp_T_slave;
        if abs(kp_in_safe) < 1e-9
            Tjsp_upper_dynamic = Tjsp_max_hard;
            Tjsp_lower_dynamic = Tjsp_min_hard;
        else
            b1 = Tjacket(i) + alpha_raw_for_inner_max / kp_in_safe;
            b2 = Tjacket(i) + alpha_raw_for_inner_min / kp_in_safe;
            Tjsp_upper_dynamic = max(b1, b2);
            Tjsp_lower_dynamic = min(b1, b2);
        end

        % Intersect dynamic feasibility with hard safety bounds
        Tjsp_upper = min(Tjsp_max_hard, Tjsp_upper_dynamic);
        Tjsp_lower = max(Tjsp_min_hard, Tjsp_lower_dynamic);
        if Tjsp_lower > Tjsp_upper
            Tjsp_lower = Tjsp_min_hard;
            Tjsp_upper = Tjsp_max_hard;
        end
        Tjsp_upper_dyn_trace(i) = Tjsp_upper;
        Tjsp_lower_dyn_trace(i) = Tjsp_lower;

        % Calculate output (pre-saturation) using previous integral
        Tjsp_ideal = term_P_Master + term_D_Master + ui_T_master(i-1) + 340 + FF_term;
        Tjsp_val = max(min(Tjsp_ideal, Tjsp_upper), Tjsp_lower);

        % Back-calculation anti-windup:
        % I_dot = Ki*e + (u_sat - u_unsat)/Tt
        Ki_T_master = Kp_T_master / max(Ti_T_master, 1e-9);
        Tt_safe = max(Tt_T_master, 1e-6);
        int_increment = (Ki_T_master * e_T_master(i) + (Tjsp_val - Tjsp_ideal) / Tt_safe) * par.tf;

        ui_T_master(i) = ui_T_master(i-1) + int_increment;
        ui_T_master(i) = max(min(ui_T_master(i), ui_T_master_max), ui_T_master_min);
        ui_T_master_trace(i) = ui_T_master(i);
        Tjsp_trace(i) = Tjsp_val;

        % -- Slave Loop (Tjsp -> u_cool / alpha) --
        e_T_slave(i) = Tjsp_val - Tjacket(i);

        alpha_val = Kp_T_slave * e_T_slave(i) + Kd_T_slave * (e_T_slave(i) - e_T_slave(i-1)) / par.tf;
        % Smooth switching logic
        k_sw = k_sw_slave;

        % Switch 1: Scale by 1/7 when alpha < 0
        s1 = 0.5 * (1 - tanh(k_sw * alpha_val));
        gain1 = 1 + (1/7 - 1) * s1;
        alpha_val = alpha_val * gain1;

        % Switch 2: Scale by 1/2 when alpha (scaled) < -1
        s2 = 0.5 * (1 - tanh(k_sw * (alpha_val + 1)));
        gain2 = 1 + (0.5 - 1) * s2;
        alpha_val = alpha_val * gain2;
        alpha_val = min(max(alpha_val, par.lbu(3)-10), par.ubu(3));

        % Demand must precede the master dynamic limit as well as the
        % physical slave limit; otherwise overload magnitude is erased.
        av = Kp_T_slave * (Tjsp_ideal - Tjacket(i));
        av = av * (1 + (1/7-1)*0.5*(1-tanh(k_sw*av)));
        av = av * (1 + (0.5-1)*0.5*(1-tanh(k_sw*(av+1))));
        u_cool_demand(i) = min(max(av,par.lbu(3)-10),par.ubu(3));
        cooling_slack_trace(i) = max(alpha_val + 1, 0);
        cooling_distance_trace(i) = abs(alpha_val + 1);
        uk(3) = max(alpha_val, -1);
    else
        % i = 1: Temperature loop not yet active, initialize to VPC setpoint
        u_cool_demand(i) = VPC_Setpoint;  % so VPC error starts near zero
        cooling_slack_trace(i) = max(u_cool_demand(i) + 1, 0);
        cooling_distance_trace(i) = abs(u_cool_demand(i) + 1);
        uk(3) = -1;  % Start with max cooling direction
    end

    % --- 3. Economic Loop / Safety Override ---
    val_FB = 0;  % shared VPC output driving both FB and Psp-relaxation
    FB_recipe_cmd = compute_fb_recipe_cmd(P(i), MAin(i), MBin(i), par.ubu(2));
    FB_feedback_cmd = 0;
    e_VPC = u_cool_demand(i) - VPC_Setpoint;
    if i < 2
        ui_VPC(i) = 0;
        if use_fb_feedback || use_psp_econ
            val_FB = K_VPC * e_VPC + ui_VPC(i);
        end
    elseif MAin(i) >= Mmax_A-1e-6
        ui_VPC(i) = ui_VPC(i-1); % Hold
        if use_fb_feedback || use_psp_econ
            val_FB = K_VPC * e_VPC + ui_VPC(i);
        end
    elseif MBin(i) >= Mlim_B
        % Initiator limit reached: stop FB, keep VPC for Psp channel only
        ui_VPC(i) = ui_VPC(i-1); % Hold
        if use_psp_econ
            val_FB = K_VPC * e_VPC + ui_VPC(i);
        end
    else
        % PI Control for VPC Loop
        % Error: cooling demand - setpoint
        % If Demand > Setpoint (e.g., -0.5 > -1), spare capacity exists -> Increase FB.
        if use_fb_feedback || use_psp_econ
            % Integral Update (Standard)
            ui_VPC(i) = ui_VPC(i-1) + (K_VPC / Ti_VPC * e_VPC) * par.tf;
            ui_VPC(i) = min(max(ui_VPC(i),ui_VPC_min),ui_VPC_max);
            % PI Output
            val_FB = K_VPC * e_VPC + ui_VPC(i);
            val_FB = min(max(val_FB, val_FB_clip_min), val_FB_clip_max);
        else
            ui_VPC(i) = ui_VPC(i-1);
        end
    end

    % Enforce the declared economic output range also in dose-hold modes.
    val_FB = min(max(val_FB,val_FB_clip_min),val_FB_clip_max);

    if MBin(i) >= Mlim_B
        % Hard stop: no more initiator once limit reached
        FB_feedback_cmd = 0;
        FB_recipe_cmd = 0;
        FB(i) = 0;
    elseif use_fb_feedback
        FB_feedback_cmd = min(max(val_FB * val_FB_to_FB, 0), par.ubu(2));
        FB(i) = FB_feedback_cmd;
    else
        FB(i) = FB_recipe_cmd;
    end
    uk(2) = FB(i);
    ui_VPC_trace(i) = ui_VPC(i);
    FB_recipe_cmd_trace(i) = FB_recipe_cmd;
    FB_feedback_cmd_trace(i) = FB_feedback_cmd;
    FB_feedback_active_trace(i) = use_fb_feedback && (FB_feedback_cmd > 1e-12);
    FB_channel_active_trace(i) = use_fb_feedback && (FB_feedback_cmd > 1e-12);

    % Psp update from the same shared VPC output.
    if switch_T_active && use_psp_econ
        Psp_limit = Psp_Nominal + K_Psp_Relax * min(val_FB, 0);
        Psp = min(max(Psp_limit,1.1e6), Psp_Nominal);
    else
        Psp = Psp_Nominal;
    end
    Psp_trace(i) = Psp;
    val_FB_trace(i) = val_FB;
    Psp_econ_active_trace(i) = switch_T_active && use_psp_econ && (val_FB < -1e-9);
    switch_P_active_trace(i) = logical(switch_P_active);

    % --- Actuator Model with First-Order Lag and Rate Limiting ---
    % u_cmd = commanded input from controller
    % u_impl = implemented (actual) input to plant
    % Implements: saturation + first-order lag + rate limit + final saturation

    % if e_VPC < -1 || switch_emergency_shutdown == 1
    %     uk(1:2) = 0; % Emergency shutdown if cooling demand is far above setpoint (safety override)
    %     switch_emergency_shutdown = 1;
    % end
    % if switch_emergency_shutdown == 1 && e_VPC > -0.8
    %     switch_emergency_shutdown = 0; % Reset shutdown if condition clears
    % end

    u_cmd = uk;  % Store commanded values
    u_cmd_actual_trace(:,i) = u_cmd;

    if i == 1
        u_impl_prev = par.initial_actuator;
    else
        u_impl_prev = U_n(:, i-1);
    end
    u_impl = zeros(3,1);
    % Apply the same actuator recurrence at every step, including startup.
        Ts = par.tf;  % Sampling time

        % FA actuator
        u_sat = min(max(u_cmd(1), par.lbu(1)), par.ubu(1));  % Step 1: saturation
        du_raw = (Ts / tau_FA) * (u_sat - u_impl_prev(1));   % Step 2: first-order lag
        du_step_max = du_max_FA * Ts;                         % Max change per step
        du = min(max(du_raw, -du_step_max), du_step_max);    % Step 3: rate limit
        u_impl(1) = min(max(u_impl_prev(1) + du, par.lbu(1)), par.ubu(1)); % Step 4: final sat

        % FB actuator
        u_sat = min(max(u_cmd(2), par.lbu(2)), par.ubu(2));
        du_raw = (Ts / tau_FB) * (u_sat - u_impl_prev(2));
        du_step_max = du_max_FB * Ts;
        du = min(max(du_raw, -du_step_max), du_step_max);
        u_impl(2) = min(max(u_impl_prev(2) + du, par.lbu(2)), par.ubu(2));

        % Alpha actuator
        u_sat = min(max(u_cmd(3), par.lbu(3)), par.ubu(3));
        du_raw = (Ts / tau_Alpha) * (u_sat - u_impl_prev(3));
        du_step_max = du_max_Alpha * Ts;
        du = min(max(du_raw, -du_step_max), du_step_max);
        u_impl(3) = min(max(u_impl_prev(3) + du, par.lbu(3)), par.ubu(3));
    % Dose meters enforce the same final clipping used by both NMPC variants.
    u_impl(1) = min(u_impl(1), max(0, Mmax_A-xk(16))/par.tf);
    u_impl(2) = min(u_impl(2), max(0, Mlim_B-xk(17))/par.tf);
    uk = u_impl;  % Use implemented (actual) input for plant simulation
    solve_times(i) = toc(t_start_cpu);

    % Dynamics Update
    U_n(:, i) = uk;
    Fk = sys.F('x0', xk, 'p', [uk; dk]);
    xk = full(Fk.xf);
    X_n(:, i) = xk;


    % Safety Stop: Explosion Prevention
    if xk(12) >= 373.15
        fprintf('Temperature trip: true reactor state %.4f K >= 373.15 K.\n', xk(12));
        % Fill remaining arrays with current value or NaN to indicate stop?
        % Actually, breaking is enough; the results calculation script should handle partial data or take the end value.
        break;
    end

    if xk(16) >= Mmax_A-1e-3 && xk(14) <= Pact_Tc, break; end

end
sim_end_idx = i;

%% 3. Visualization (Standardized)
Tr_traj = X_n(12, 1:sim_end_idx);
P_traj = X_n(14, 1:sim_end_idx);
Psp = 1.5e6;
% Metrics Calculation (Moved BEFORE Table usage)
[~, ind_T_exc] = max(Tr_traj >= Trsp); if isempty(ind_T_exc), ind_T_exc=1; end
Tmax = max(abs(Tr_traj(ind_T_exc:end) - Trsp));
Tmean = mean(abs(Tr_traj(ind_T_exc:end) - Trsp));
[~, iPs] = max(P_traj >= Psp); if isempty(iPs), iPs=1; end
[~, iPe] = max(MAin(1:sim_end_idx) >= Mmax_A); if isempty(iPe), iPe=sim_end_idx; end
Pmax = max(abs(P_traj(iPs:iPe) - Psp)/1e5);
Pmean = mean(abs(P_traj(iPs:iPe) - Psp)/1e5);

mP_final = X_n(5, end);
PDI_final = PDI(sim_end_idx);
Mn_final = Mn(sim_end_idx);
Time_hours = sim_end_idx / 3600; % Defined here

% IAE and ISE Calculations
IAE_P = sum(abs(e_P(1:sim_end_idx))) * par.tf;
ISE_P = sum(e_P(1:sim_end_idx).^2) * par.tf;
IAE_T = sum(abs(e_T_master(1:sim_end_idx))) * par.tf;
ISE_T = sum(e_T_master(1:sim_end_idx).^2) * par.tf;

performance = table(Time_hours, mP_final, PDI_final, Mn_final, Tmean, Tmax, Pmean, Pmax, IAE_P, ISE_P, IAE_T, ISE_T);

% Reconstruct Reference Trajectory for Plotting
Trsp_final = 273 + 78;
Tr_start_val = Tr(1);
Trsp_ref = zeros(1, sim_end_idx);
% Ramp_Steps from Init section

for t_Step = 1:sim_end_idx
    if t_Step <= Ramp_Steps
        tau = (t_Step-1) / Ramp_Steps;
        % Minimum-jerk trajectory (5th order polynomial)
        beta = tau^3 * (10 - 15*tau + 6*tau^2);
        Trsp_ref(t_Step) = Tr_start_val + (Trsp_final - Tr_start_val) * beta;
    else
        Trsp_ref(t_Step) = Trsp_final;
    end
end

% Only plot if caller explicitly requests the fig output (nargout == 2)
% When called as [perf, ~, SimData] the fig is unused, skip plotting.
if nargout == 2 || (isfield(par, 'plot_fig') && par.plot_fig)
    fig = figure(101); clf;
    t = 1:sim_end_idx;

    subplot(4,2,1); plot(t, FA(1:sim_end_idx)); title('FA'); grid on;
    subplot(4,2,2); plot(t, FB(1:sim_end_idx)); title('FB'); grid on;
    subplot(4,2,3); plot(t, u_cool_demand(1:sim_end_idx), 'b', t, -1*ones(size(t)), 'r--'); title('Cooling Demand'); grid on;
    subplot(4,2,4); plot(t, P(1:sim_end_idx), 'b', t, Psp_trace(1:sim_end_idx), 'k--'); title('Pressure'); grid on;
    subplot(4,2,5);
    plot(t, Tr(1:sim_end_idx), 'b', 'LineWidth', 1.5); hold on;
    plot(t, Trsp_ref, 'g--', 'LineWidth', 1.5);
    title('Temperature'); grid on; ylim([330 360]); legend('Tr', 'Tr_{sp}');
    subplot(4,2,6); plot(t, RA(1:sim_end_idx)); title('Reaction Rate'); grid on;
    subplot(4,2,7); plot(t, Mw(1:sim_end_idx)); title('weight-average molecular weight'); grid on;
    subplot(4,2,8); plot(t, k70(1:sim_end_idx)); title('k70'); grid on;
else
    fig = [];
end

if nargout > 2
    ctrl_params = struct( ...
        'arc_variant', arc_variant, ...
        'use_fb_feedback', logical(use_fb_feedback), ...
        'use_psp_econ', logical(use_psp_econ), ...
        'Kp_P', Kp_P, ...
        'Ti_P', Ti_P, ...
        'Td_P', Td_P, ...
        'Kp_T_master', Kp_T_master, ...
        'Ti_T_master', Ti_T_master, ...
        'Td_T_master', Td_T_master, ...
        'Tt_T_master', Tt_T_master, ...
        'Kp_T_slave', Kp_T_slave, ...
        'Kd_T_slave', Kd_T_slave, ...
        'k_sw_slave', k_sw_slave, ...
        'K_FF_P', K_FF_P, ...
        'K_VPC', K_VPC, ...
        'Ti_VPC', Ti_VPC, ...
        'VPC_Setpoint', VPC_Setpoint, ...
        'ui_VPC_min', ui_VPC_min, ...
        'ui_VPC_max', ui_VPC_max, ...
        'val_FB_clip_min', val_FB_clip_min, ...
        'val_FB_clip_max', val_FB_clip_max, ...
        'val_FB_to_FB', val_FB_to_FB, ...
        'Psp_Nominal', Psp_Nominal, ...
        'K_Psp_Relax', K_Psp_Relax, ...
        'VPC_Psp_Overload_Gain', VPC_Psp_Overload_Gain, ...
        'tau_FA', tau_FA, ...
        'tau_FB', tau_FB, ...
        'tau_Alpha', tau_Alpha, ...
        'du_max_FA', du_max_FA, ...
        'du_max_FB', du_max_FB, ...
        'du_max_Alpha', du_max_Alpha, ...
        'FA_max', par.ubu(1), ...
        'FB_max', par.ubu(2), ...
        'alpha_min', par.lbu(3), ...
        'alpha_max', par.ubu(3), ...
        'Mmax_A', Mmax_A, ...
        'Mlim_B', Mlim_B, ...
        'Pact_Tc', Pact_Tc, ...
        'sample_time', par.tf);

    t_seconds = (1:sim_end_idx) * par.tf;
    thermal_active_mask = build_thermal_active_mask( ...
        t_seconds, 3300, true(1, sim_end_idx), 1200);
    productive_mask = MAin(1:sim_end_idx) < Mmax_A;
    boundary_mask = U_n(3, 1:sim_end_idx) <= -0.999;
    near_boundary_05_mask = cooling_distance_trace(1:sim_end_idx) <= 0.05;
    near_boundary_10_mask = cooling_distance_trace(1:sim_end_idx) <= 0.10;

    SimData.t = 1:sim_end_idx;
    SimData.Time_hours = (1:sim_end_idx) * par.tf / 3600;
    SimData.variant_name = arc_variant;
    SimData.Tr = Tr(1:sim_end_idx);
    SimData.Trsp = Trsp_ref(1:sim_end_idx);
    SimData.P = P(1:sim_end_idx);
    SimData.Psp = Psp_trace(1:sim_end_idx);
    SimData.Tjacket = Tjacket(1:sim_end_idx);  % Jacket temperature
    SimData.FA = FA(1:sim_end_idx);
    SimData.FB = FB(1:sim_end_idx);
    SimData.Fgl = Fgl(1:sim_end_idx);
    SimData.val_FB = val_FB_trace(1:sim_end_idx);
    SimData.ui_P = ui_P_trace(1:sim_end_idx);
    SimData.ui_T_master = ui_T_master_trace(1:sim_end_idx);
    SimData.ui_VPC = ui_VPC_trace(1:sim_end_idx);
    SimData.switch_P_active = switch_P_active_trace(1:sim_end_idx);
    SimData.FB_recipe_cmd = FB_recipe_cmd_trace(1:sim_end_idx);
    SimData.FB_feedback_cmd = FB_feedback_cmd_trace(1:sim_end_idx);
    SimData.Psp_econ_active_mask = Psp_econ_active_trace(1:sim_end_idx);
    SimData.FB_feedback_active_mask = FB_feedback_active_trace(1:sim_end_idx);
    SimData.FB_channel_active_mask = FB_channel_active_trace(1:sim_end_idx);
    SimData.alpha = u_cool_demand(1:sim_end_idx);
    SimData.cooling_slack = cooling_slack_trace(1:sim_end_idx);
    SimData.cooling_distance = cooling_distance_trace(1:sim_end_idx);
    SimData.thermal_active_mask = thermal_active_mask;
    SimData.productive_mask = productive_mask;
    SimData.boundary_mask = boundary_mask;
    SimData.near_boundary_05_mask = near_boundary_05_mask;
    SimData.near_boundary_10_mask = near_boundary_10_mask;
    SimData.Tjsp = Tjsp_trace(1:sim_end_idx);
    SimData.Tjsp_upper_dyn = Tjsp_upper_dyn_trace(1:sim_end_idx);
    SimData.Tjsp_lower_dyn = Tjsp_lower_dyn_trace(1:sim_end_idx);
    % Control signals (for experiment framework)
    SimData.u_cmd = [FA(1:sim_end_idx); FB(1:sim_end_idx); u_cool_demand(1:sim_end_idx)];
    SimData.u_impl = U_n(:, 1:sim_end_idx);
    SimData.u_cmd_actual = u_cmd_actual_trace(:,1:sim_end_idx);  % Actual applied (after actuator)
    SimData.Mn = Mn(1:sim_end_idx);
    SimData.PDI = PDI(1:sim_end_idx);
    SimData.mP_final = mP_final;
    SimData.kt = kt(1:sim_end_idx);  % Termination rate (k7*lambda0^2)
    SimData.Vf = Vf(1:sim_end_idx);  % Free volume fraction
    SimData.k70 = k70(1:sim_end_idx);
    SimData.solve_times = solve_times(1:sim_end_idx);
    SimData.ctrl_params = ctrl_params;
    SimData.initial_state = par.x0;
    SimData.X = X_n(:, 1:sim_end_idx); % Full state for theoretical analysis
    performance.solve_times = solve_times(1:sim_end_idx);  % Pre-exponential factor (Gel Effect modified)
end

end

function FB_recipe_cmd = compute_fb_recipe_cmd(P_meas, MAin, MBin, FB_max)
P_recipe_min = 1.4e6;
MA_recipe_cutoff = 2950;
MB_mid = 0.0012;
MB_max = 0.0020;
FB_high = min(4e-7, FB_max);
FB_low = min(1e-7, FB_max);

if P_meas < P_recipe_min || MAin > MA_recipe_cutoff || MBin >= MB_max
    FB_recipe_cmd = 0;
elseif MBin < MB_mid
    FB_recipe_cmd = FB_high;
else
    FB_recipe_cmd = FB_low;
end
end

function alpha_raw = invert_slave_alpha_cmd(target_cmd, lbu_alpha, ubu_alpha, k_sw)
% Find raw slave-controller output that produces a desired commanded alpha.
target = min(max(target_cmd, max(lbu_alpha, -1)), ubu_alpha);
lo = -200;
hi = 200;
f_lo = slave_alpha_cmd_from_raw(lo, lbu_alpha, ubu_alpha, k_sw) - target;
f_hi = slave_alpha_cmd_from_raw(hi, lbu_alpha, ubu_alpha, k_sw) - target;
if f_lo > 0
    alpha_raw = lo;
    return;
end
if f_hi < 0
    alpha_raw = hi;
    return;
end
for it = 1:60
    mid = 0.5 * (lo + hi);
    f_mid = slave_alpha_cmd_from_raw(mid, lbu_alpha, ubu_alpha, k_sw) - target;
    if f_mid >= 0
        hi = mid;
    else
        lo = mid;
    end
end
alpha_raw = 0.5 * (lo + hi);
end

function alpha_cmd = slave_alpha_cmd_from_raw(alpha_raw, lbu_alpha, ubu_alpha, k_sw)
s1 = 0.5 * (1 - tanh(k_sw * alpha_raw));
gain1 = 1 + (1/7 - 1) * s1;
alpha_scaled = alpha_raw * gain1;
s2 = 0.5 * (1 - tanh(k_sw * (alpha_scaled + 1)));
gain2 = 1 + (0.5 - 1) * s2;
alpha_scaled = alpha_scaled * gain2;
alpha_scaled = min(max(alpha_scaled, lbu_alpha - 10), ubu_alpha);
alpha_cmd = max(alpha_scaled, -1);
end

function mask = build_thermal_active_mask(t_seconds, t_start, base_mask, ~)
% build_thermal_active_mask  Mark samples as thermally active.
%   mask(k) = true iff t_seconds(k) >= t_start AND base_mask(k).
    mask = (t_seconds >= t_start) & base_mask;
end
