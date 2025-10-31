% Authors: Teni Adekunle Awomuti, Purav Alva, Tate Dailey
% Course: ME 4053 – Wind Turbine Analysis
% Date: October 30th 2025

% Group Parameters:
%   Q1 (Deliverable 1): Test case -> U = 10 m/s, rpm = 14, pitch = 0 deg
%   Q2 (Deliverable 2): U = 6.6 m/s, TSR λ = 6.09 (compute rpm from λ)
%   Q3 (Deliverable 3): U = 7.1 m/s (sweep β and λ, rpm cap 15.5)
%   Q4 (Deliverable 4): Rated control around U = 19.4 m/s (rpm fixed 15.5)
%   Q5 (Deliverable 5): Tower analysis uses Fz from D4 (here kept at 1.65e5 N)
%
% Dependencies in same folder:
%   BladeProfile.csv
%   DU91-W2-250.csv, DU93-W-210.csv, DU96-W-180.csv, DU97-W-300.csv
%   towerSpecs.csv
%
% Outputs:
%   D1_section_outputs.csv
%   D4_pitch_vs_wind.csv
% -------------------------------------------------------------------------

clear; clc; close all;

%% --------------------- GLOBAL SETTINGS (shared) ------------------------
rho = 1.1;              % air density [kg/m^3]
mu  = 1.8e-5;           % dynamic viscosity [Pa·s]
B   = 3;                % blades
bladeFile  = 'BladeProfile.csv';
polarFiles = {'DU91-W2-250.csv','DU93-W-210.csv','DU96-W-180.csv','DU97-W-300.csv'};
% BEM numerical params
tol = 1e-5; max_it = 300; relax = 0.3; hub_cut_frac = 0.05; cap_a = 0.95;

% Small helper for polars key
mapKey = @(s) lower(strrep(s,'-',''));

%% -------------------------- DELIVERABLE 1 ------------------------------
% Baseline test case: U = 10 m/s, rpm = 14, pitch = 0 deg
U_inf = 10; rpm = 14; pitch_deg = 0;

% Geometry
T = readtable(bladeFile);
r  = T.DistanceFromCenterOfRotation/1000;     % [m]
c  = T.ChordLength/1000;                      % [m]
tw = deg2rad(T.BladeTwist);                   % [rad]
airfoils = strtrim(string(T.Airfoil));

% Sort + cut hub
[r, idx] = sort(r); c = c(idx); tw = tw(idx); airfoils = airfoils(idx);
R = max(r); A = pi*R^2; omega = rpm*2*pi/60; lambda = omega*R/U_inf;
mask = r > hub_cut_frac*R;
r = r(mask); c = c(mask); tw = tw(mask); airfoils = airfoils(mask);
N = numel(r);

% Load polars
polars = struct();
for k = 1:numel(polarFiles)
    f = polarFiles{k};
    if ~isfile(f), warning('Missing polar file %s',f); continue; end
    D = readtable(f);
    [~,base] = fileparts(f);
    key = mapKey(base);
    polars.(key).alpha = D.AoA;
    polars.(key).Cl    = D.CL;
    polars.(key).Cd    = D.CD;
end

% Spanwise BEM
dr = diff(r); dr = [dr; dr(end)];
pitch = deg2rad(pitch_deg);
a_vec = zeros(N,1); ap_vec = zeros(N,1);
dT = zeros(N,1); dQ = zeros(N,1);
Cl_v = zeros(N,1); Cd_v = zeros(N,1);

for i = 1:N
    ri = r(i); ci = c(i); twi = tw(i);
    air = airfoils(i); isCircle = strcmpi(air,'circle');
    polKey = mapKey(air); hasPolar = isfield(polars,polKey);
    a = 0.3; ap = 0;
    for it = 1:max_it
        Vax = U_inf*(1-a); Vtan = omega*ri*(1+ap);
        phi = atan2(Vax,Vtan);
        alpha_deg = rad2deg(phi - (twi + pitch));
        if ~isCircle && hasPolar
            P = polars.(polKey);
            Cl = interp1(P.alpha,P.Cl,alpha_deg,'linear','extrap');
            Cd = interp1(P.alpha,P.Cd,alpha_deg,'linear','extrap');
        elseif ~isCircle && ~hasPolar
            [Cl,Cd] = fallbackPolar(alpha_deg);
        else
            Re = rho*U_inf*ci/mu; Cl = 0; Cd = cylinderCD(Re);
        end
        Cn = Cl*cos(phi) + Cd*sin(phi);
        Ct_loc = Cl*sin(phi) - Cd*cos(phi);
        sigma = B*ci/(2*pi*ri);
        f = (B/2)*(R-ri)/(ri*max(abs(sin(phi)),1e-6));
        F = (2/pi)*acos(exp(-f)); F = max(F,1e-6);
        a_new  = 1/(1+(4*F*sin(phi)^2)/(sigma*max(Cn,1e-6)));
        ap_new = 1/((4*F*sin(phi)*cos(phi))/(sigma*max(Ct_loc,1e-6))-1);
        a  = (1-relax)*a  + relax*a_new;  a  = min(max(a,0),cap_a);
        ap = (1-relax)*ap + relax*ap_new; ap = min(max(ap,-0.5),1.0);
        if abs(a-a_new)<tol && abs(ap-ap_new)<tol, break; end
    end
    Vrel = hypot(U_inf*(1-a), omega*ri*(1+ap));
    dL = 0.5*rho*Vrel^2*ci*Cl; dD = 0.5*rho*Vrel^2*ci*Cd;
    dFx = dL*cos(phi) + dD*sin(phi);
    dFt = (~isCircle)*(dL*sin(phi) - dD*cos(phi));
    dT(i) = B*dFx; dQ(i) = B*dFt*ri;
    Cl_v(i)=Cl; Cd_v(i)=Cd; a_vec(i)=a; ap_vec(i)=ap;
end

% Totals
T_total = trapz(r,dT);
Q_total = trapz(r,dQ);
P1 = omega*Q_total;
Cp1 = P1/(0.5*rho*A*U_inf^3);
Ct1 = T_total/(0.5*rho*A*U_inf^2);

% Print & plot & save
disp('--- Deliverable 1 ---');
disp(['TSR = ', num2str(lambda)]);
disp(['Cp = ', num2str(Cp1), ', Ct = ', num2str(Ct1)]);

figure;
subplot(1,2,1);
plot(r,dQ,'-o','LineWidth',1.5); grid on; box on;
xlabel('r (m)'); ylabel('dQ (N·m)');
title('Deliverable 1 – Torque Distribution');
subplot(1,2,2);
plot(r,Cl_v,'-o','LineWidth',1.5); hold on;
plot(r,Cd_v,'-s','LineWidth',1.5); grid on; box on;
xlabel('r (m)'); ylabel('Coeff');
legend('C_l','C_d','Location','best');
title('Deliverable 1 – Section Cl & Cd');

sec = table(r,c,rad2deg(tw),airfoils,a_vec,ap_vec,Cl_v,Cd_v,dT,dQ, ...
    'VariableNames',{'r_m','chord_m','twist_deg','airfoil','a','ap','Cl','Cd','dT','dQ'});
writetable(sec, fullfile(pwd,'D1_section_outputs.csv'));

%% -------------------------- DELIVERABLE 2 ------------------------------
% Group parameters: U = 6.6 m/s, TSR λ = 6.09
U_inf = 6.6; lambda = 6.09;
% Geometry (for R)
T = readtable(bladeFile);
R_all = max(T.DistanceFromCenterOfRotation)/1000;
omega = lambda*U_inf/R_all; rpm = omega*60/(2*pi);

% Rebuild geometry arrays
r  = T.DistanceFromCenterOfRotation/1000; c  = T.ChordLength/1000;
tw = deg2rad(T.BladeTwist); airfoils = strtrim(string(T.Airfoil));
[r,idx] = sort(r); c = c(idx); tw = tw(idx); airfoils = airfoils(idx);
R = max(r); A = pi*R^2;
mask = r > hub_cut_frac*R; r=r(mask); c=c(mask); tw=tw(mask); airfoils=airfoils(mask);

% Load polars
polars = struct();
for k = 1:numel(polarFiles)
    f = polarFiles{k};
    if ~isfile(f), warning('Missing polar file %s',f); continue; end
    D = readtable(f);
    [~,base] = fileparts(f);
    key = mapKey(base);
    polars.(key).alpha = D.AoA; polars.(key).Cl = D.CL; polars.(key).Cd = D.CD;
end

pitch_range = -3:0.3:20;
Cp_list = zeros(length(pitch_range),1);
for p = 1:length(pitch_range)
    Cp_list(p) = computeCp(U_inf, pitch_range(p), omega, r,c,tw,airfoils,polars,mapKey,B,rho,mu,R,A,tol,max_it,relax,cap_a);
end
[maxCp2,idx2] = max(Cp_list); bestPitch2 = pitch_range(idx2);
disp('--- Deliverable 2 ---');
disp(['Max Cp = ', num2str(maxCp2), ' at pitch = ', num2str(bestPitch2), ' deg']);

figure;
plot(pitch_range,Cp_list,'-o','LineWidth',1.5); grid on; box on;
xlabel('Pitch \beta (deg)'); ylabel('C_p');
title(sprintf('Deliverable 2 – C_p vs Pitch (U=%.2f m/s, \\lambda=%.2f, rpm=%.2f)',U_inf,lambda,rpm));

%% -------------------------- DELIVERABLE 3 ------------------------------
% Group parameter: U = 7.1 m/s
U_inf = 7.1;
pitch_range = -2:1:20; lambda_range = 3:0.5:9; rpm_max = 15.5;

% Geometry
T = readtable(bladeFile);
r  = T.DistanceFromCenterOfRotation/1000; c  = T.ChordLength/1000;
tw = deg2rad(T.BladeTwist); airfoils = strtrim(string(T.Airfoil));
[r,idx] = sort(r); c=c(idx); tw=tw(idx); airfoils=airfoils(idx);
R = max(r); A = pi*R^2; mask=r>hub_cut_frac*R;
r=r(mask); c=c(mask); tw=tw(mask); airfoils=airfoils(mask);

% Polars
polars = struct();
for k = 1:numel(polarFiles)
    f = polarFiles{k};
    if ~isfile(f), warning('Missing polar file %s',f); continue; end
    D = readtable(f);
    [~,base] = fileparts(f);
    key = mapKey(base);
    polars.(key).alpha = D.AoA; polars.(key).Cl = D.CL; polars.(key).Cd = D.CD;
end

Cp_mat = nan(length(pitch_range), length(lambda_range));
for ip = 1:length(pitch_range)
    for il = 1:length(lambda_range)
        lambda = lambda_range(il);
        omega = lambda*U_inf/R;
        rpm = omega*60/(2*pi);
        if rpm > rpm_max, continue; end
        [Cp_mat(ip,il),~] = runBEM_once(U_inf,r,c,tw,airfoils,polars,mapKey,B,rho,mu,pitch_range(ip),omega,R,A,tol,max_it,relax,cap_a);
    end
end
[maxCp3, idx3] = max(Cp_mat(:));
[ip3, il3] = ind2sub(size(Cp_mat),idx3);
disp('--- Deliverable 3 ---');
disp(['Max Cp = ', num2str(maxCp3), ' at pitch = ', num2str(pitch_range(ip3)), ', TSR = ', num2str(lambda_range(il3))]);

figure;
[X,Y] = meshgrid(lambda_range,pitch_range);
surf(X,Y,Cp_mat,'EdgeColor','none'); grid on; box on; colorbar;
xlabel('TSR \lambda'); ylabel('Pitch \beta (deg)'); zlabel('C_p');
title(sprintf('Deliverable 3 – C_p(\\lambda, \\beta) Surface (U=%.2f m/s)',U_inf));

%% -------------------------- DELIVERABLE 4 ------------------------------
% Rated control: rpm fixed 15.5, power 2.5 MW, sweep U = 11..25 m/s
P_rated = 2.5e6; rpm_fix = 15.5; omega_fix = rpm_fix*2*pi/60;
U_list = 11:1:25;  % includes 19.4 m/s

% Geometry (once)
T = readtable(bladeFile);
r_all = T.DistanceFromCenterOfRotation/1000; c_all = T.ChordLength/1000;
tw_all = deg2rad(T.BladeTwist); af_all = strtrim(string(T.Airfoil));
[r_all,idx] = sort(r_all); c_all=c_all(idx); tw_all=tw_all(idx); af_all=af_all(idx);
R = max(r_all); A = pi*R^2;
mask = r_all > hub_cut_frac*R;
r = r_all(mask); c = c_all(mask); tw = tw_all(mask); afs = af_all(mask);

% Polars (once)
polars = struct();
for k = 1:numel(polarFiles)
    f = polarFiles{k};
    if ~isfile(f), warning('Missing polar file %s',f); continue; end
    D = readtable(f);
    [~,base] = fileparts(f);
    key = mapKey(base);
    polars.(key).alpha = D.AoA; polars.(key).Cl = D.CL; polars.(key).Cd = D.CD;
end

Pitch_req = nan(size(U_list));
Cp_at_req = nan(size(U_list));
Ct_at_req = nan(size(U_list));
P_at_req  = nan(size(U_list));
T_at_req  = nan(size(U_list));
lambda_at_req = nan(size(U_list));

for k = 1:numel(U_list)
    U = U_list(k);
    getPowerAtPitch = @(pitch_deg) bem_power_at_pitch(U,pitch_deg,omega_fix, ...
        r,c,tw,afs,polars,mapKey,B,rho,mu,R,A,tol,max_it,relax,cap_a);
    [Cp0,Ct0,P0,T0] = getPowerAtPitch(0);
    [CpH,CtH,PH,TH] = getPowerAtPitch(25);
    if P0 <= P_rated
        Pitch_req(k)=0; Cp_at_req(k)=Cp0; Ct_at_req(k)=Ct0; P_at_req(k)=P0; T_at_req(k)=T0;
        lambda_at_req(k)=(omega_fix*R)/U; continue;
    end
    if PH > P_rated
        Pitch_req(k)=NaN; Cp_at_req(k)=NaN; Ct_at_req(k)=NaN; P_at_req(k)=PH; T_at_req(k)=TH;
        lambda_at_req(k)=(omega_fix*R)/U; continue;
    end
    lo=0; hi=25; CpM=NaN; CtM=NaN; PM=NaN; Tmid=NaN;
    for iter=1:40
        mid=0.5*(lo+hi);
        [CpM,CtM,PM,Tmid] = getPowerAtPitch(mid);
        if PM > P_rated, lo = mid; else, hi = mid; end
        if abs(hi-lo) < 1e-3, break; end
    end
    Pitch_req(k)=hi; Cp_at_req(k)=CpM; Ct_at_req(k)=CtM; P_at_req(k)=PM; T_at_req(k)=Tmid;
    lambda_at_req(k)=(omega_fix*R)/U;
end

disp('--- Deliverable 4 ---');
disp('U (m/s)   Pitch (deg)   Cp     Ct     P (MW)   TSR    Thrust (N)');
for k=1:numel(U_list)
    fprintf('%6.1f   %11.3f  %6.4f %6.4f  %7.3f  %6.3f  %10.1f\n', ...
        U_list(k), Pitch_req(k), Cp_at_req(k), Ct_at_req(k), ...
        P_at_req(k)/1e6, lambda_at_req(k), T_at_req(k));
end

T_out = table(U_list(:),Pitch_req(:),Cp_at_req(:),Ct_at_req(:), ...
              P_at_req(:),lambda_at_req(:),T_at_req(:), ...
              'VariableNames',{'U_ms','pitch_deg','Cp','Ct','P_W','lambda','ThrustN'});
writetable(T_out, fullfile(pwd,'D4_pitch_vs_wind.csv'));

%% -------------------------- DELIVERABLE 5 ------------------------------
% Tower structural & fatigue (uses fixed Fz from D4)
WindTurbineTower_Analysis();   % prints max deflection explicitly

%% ============================= HELPERS ================================
function [Cl,Cd] = fallbackPolar(alpha_deg)
a = deg2rad(alpha_deg);
Cl = 2*pi*a; Cl = max(min(Cl,1.2),-1.2);
Cd = 0.01 + 0.02*(abs(alpha_deg)/10);
end

function C_D = cylinderCD(Re)
% Empirical drag for smooth cylinders across regimes.
if numel(Re) > 1, C_D = arrayfun(@cylinderCD, Re); return; end
if Re < 2e5
    C_D = 11*Re.^(-0.75) + 0.9*(1-exp(-1000./max(Re,1))) + 1.2*(1-exp(-(Re./4500).^0.7));
elseif Re <= 5e5
    C_D = 10.^(0.32*tanh(44.45-8*log10(Re))-0.239);
else
    C_D = 0.1*log10(Re)-0.253;
end
end

function [Cp,Ct] = runBEM_once(U_inf,r,c,tw,airfoils,polars,mapKey, ...
                               B,rho,mu,pitch_deg,omega,R,A,tol,max_it,relax,cap_a)
N = numel(r);
pitch = deg2rad(pitch_deg);
dT = zeros(N,1); dQ = zeros(N,1);
for i = 1:N
    ri = r(i); ci = c(i); twi = tw(i); air = airfoils(i);
    isCircle = strcmpi(air,'circle');
    polKey = mapKey(air); hasPolar = isfield(polars,polKey);
    a=0.3; ap=0;
    for it=1:max_it
        Vax = U_inf*(1-a); Vtan = omega*ri*(1+ap); phi = atan2(Vax,Vtan);
        alpha_deg = rad2deg(phi - (twi + pitch));
        if ~isCircle && hasPolar
            P = polars.(polKey);
            Cl = interp1(P.alpha,P.Cl,alpha_deg,'linear','extrap');
            Cd = interp1(P.alpha,P.Cd,alpha_deg,'linear','extrap');
        elseif ~isCircle && ~hasPolar
            [Cl,Cd] = fallbackPolar(alpha_deg);
        else
            Re = rho*U_inf*ci/mu; Cl = 0; Cd = cylinderCD(Re);
        end
        Cn = Cl*cos(phi) + Cd*sin(phi);
        Ct_loc = Cl*sin(phi) - Cd*cos(phi);
        sigma = B*ci/(2*pi*ri);
        f = (B/2)*(R-ri)/(ri*max(abs(sin(phi)),1e-6));
        F = (2/pi)*acos(exp(-f)); F = max(F,1e-6);
        a_new = 1/(1+(4*F*sin(phi)^2)/(sigma*max(Cn,1e-6)));
        ap_new= 1/((4*F*sin(phi)*cos(phi))/(sigma*max(Ct_loc,1e-6))-1);
        a  = (1-relax)*a  + relax*a_new;  a  = min(max(a,0),cap_a);
        ap = (1-relax)*ap + relax*ap_new; ap = min(max(ap,-0.5),1.0);
        if abs(a-a_new)<tol && abs(ap-ap_new)<tol, break; end
    end
    Vrel = hypot(U_inf*(1-a), omega*ri*(1+ap));
    dL = 0.5*rho*Vrel^2*ci*Cl; dD = 0.5*rho*Vrel^2*ci*Cd;
    dFx = dL*cos(phi) + dD*sin(phi);
    dFt = (~isCircle)*(dL*sin(phi) - dD*cos(phi));
    dT(i) = B*dFx; dQ(i) = B*dFt*ri;
end
T_total = trapz(r,dT); Q_total = trapz(r,dQ);
P = omega*Q_total;
Cp = P/(0.5*rho*A*U_inf^3);
Ct = T_total/(0.5*rho*A*U_inf^2);
end

function Cp = computeCp(U,pitch_deg,omega,r,c,tw,airfoils,polars,mapKey, ...
                         B,rho,mu,R,A,tol,max_it,relax,cap_a)
[Cp,~] = runBEM_once(U,r,c,tw,airfoils,polars,mapKey, ...
                     B,rho,mu,pitch_deg,omega,R,A,tol,max_it,relax,cap_a);
end

function [Cp,Ct,P,T_total] = bem_power_at_pitch(U,pitch_deg,omega, ...
    r,c,tw,airfoils,polars,mapKey,B,rho,mu,R,A,tol,max_it,relax,cap_a)
N=numel(r); pitch=deg2rad(pitch_deg);
dT=zeros(N,1); dQ=zeros(N,1);
for i=1:N
    ri=r(i); ci=c(i); twi=tw(i); air=airfoils(i);
    isCircle=strcmpi(air,'circle'); polKey=mapKey(air); hasPolar=isfield(polars,polKey);
    a=0.3; ap=0;
    for it=1:max_it
        Vax=U*(1-a); Vtan=omega*ri*(1+ap); phi=atan2(Vax,Vtan);
        alpha_deg=rad2deg(phi-(twi+pitch));
        if ~isCircle && hasPolar
            Ptab=polars.(polKey);
            Cl=interp1(Ptab.alpha,Ptab.Cl,alpha_deg,'linear','extrap');
            Cd=interp1(Ptab.alpha,Ptab.Cd,alpha_deg,'linear','extrap');
        elseif ~isCircle && ~hasPolar
            [Cl,Cd]=fallbackPolar(alpha_deg);
        else
            Re=rho*U*ci/mu; Cl=0; Cd=cylinderCD(Re);
        end
        Cn=Cl*cos(phi)+Cd*sin(phi); Ct_loc=Cl*sin(phi)-Cd*cos(phi);
        sigma=B*ci/(2*pi*ri);
        f=(B/2)*(R-ri)/(ri*max(abs(sin(phi)),1e-6));
        F=(2/pi)*acos(exp(-f)); F=max(F,1e-6);
        a_new=1/(1+(4*F*sin(phi)^2)/(sigma*max(Cn,1e-6)));
        ap_new=1/((4*F*sin(phi)*cos(phi))/(sigma*max(Ct_loc,1e-6))-1);
        a=(1-relax)*a+relax*a_new; a=min(max(a,0),cap_a);
        ap=(1-relax)*ap+relax*ap_new; ap=min(max(ap,-0.5),1.0);
        if abs(a-a_new)<tol && abs(ap-ap_new)<tol, break; end
    end
    Vrel=hypot(U*(1-a),omega*ri*(1+ap));
    dL=0.5*rho*Vrel^2*ci*Cl; dD=0.5*rho*Vrel^2*ci*Cd;
    dFx=dL*cos(phi)+dD*sin(phi);
    dFt=(~isCircle)*(dL*sin(phi)-dD*cos(phi));
    dT(i)=B*dFx; dQ(i)=B*dFt*ri;
end
T_total=trapz(r,dT); Q_total=trapz(r,dQ);
P=omega*Q_total; 
Cp=P/(0.5*rho*A*U^3); 
Ct=T_total/(0.5*rho*A*U^2);
end

function WindTurbineTower_Analysis()
% Steel tower with distributed wind load + hub thrust.
rho = 1.1; mu = 1.8e-5;
U_ref = 19.4; x_ref = 77.7; e = 1/7;
E = 200e9; Fz = 1.65e5; Hub_H = 80.4;

% Geometry from CSV (mm -> m)
data = readmatrix('towerSpecs.csv');
x_data = data(:,1)/1000;  OD = data(:,2)/1000;  t = data(:,3)/1000;

% Interpolate to a fine grid
x = linspace(0.001, Hub_H, 500)'; 
D = interp1(x_data, OD, x, 'linear', 'extrap');
thk = interp1(x_data, t,  x, 'linear', 'extrap');

% Wind profile and loads
U = (x/x_ref).^e .* U_ref;
Re = (rho.*U.*D)./mu;
C_D = cylinderCD(Re);
q = 0.5.*rho.*C_D.*U.^2.*D;      % N/m
I = (pi/64).*(D.^4 - (D - 2*thk).^4);

% Shear, moment, slope, deflection (cantilever base at x=0)
V_rxn = trapz(x,q) + Fz;
M_rxn = trapz(x, q.*x) + Fz*((Hub_H-77.7)/2 + 77.7);
V = cumtrapz(x, -q) + V_rxn;
M = cumtrapz(x, V) - M_rxn;
theta = cumtrapz(x, M./(E.*I));
delta = cumtrapz(x, theta);

% Bending stress
sigma_b = (M .* (D/2))./I; sigma_b_MPa = sigma_b/1e6;
sigma_max = max(abs(sigma_b_MPa));

% Base axial stress from thrust
D_o = D(1); t_b = thk(1); D_i = D_o - 2*t_b;
A_base = (pi/4)*(D_o^2 - D_i^2);
sigma_axial_base_MPa = (Fz / A_base)/1e6;

% Safety factors
Sy = 345; Su = 450; Ssy = 0.5*Sy;
sigma_combined_base = max(sigma_b_MPa) + sigma_axial_base_MPa;
sigma_vm_base = sqrt(sigma_combined_base.^2);
FS_VM = Sy / max(sigma_vm_base,1e-12);
FS_MaxNormal = Sy / max(abs(sigma_combined_base),1e-12);
FS_MaxShear  = Ssy / max(abs(sigma_combined_base)/2,1e-12);

% Simple Goodman fatigue estimate
min_Stress = sigma_max * cosd(315 - 157.5);
mean_Stress = (sigma_max + min_Stress)/2;
alternating_Stress = abs(sigma_max - mean_Stress)/2;
sigma_e = 0.5*Su;
FS_goodman = 1 / ((alternating_Stress/sigma_e) + (mean_Stress/Su));

% Print concise summary (now includes deflection line for your report)
disp('--- Deliverable 5 ---');
fprintf('Max tower deflection = %.6f m\n', max(abs(delta)));
disp(['Max bending stress = ', num2str(sigma_max), ' MPa']);
disp(['FS_yield ~ ', num2str(Sy / max(sigma_max,1e-12))]);
disp(['FS_ultimate ~ ', num2str(Su / max(sigma_max,1e-12))]);
disp(['Goodman FS ~ ', num2str(FS_goodman)]);

% Plots
figure;
subplot(5,1,1); plot(x,q,'LineWidth',1.5); grid on; box on;
ylabel('Load (N/m)'); title('Deliverable 5 – Distributed Drag');
subplot(5,1,2); plot(x,V,'LineWidth',1.5); grid on; box on;
ylabel('Shear (N)'); title('Shear Force V');
subplot(5,1,3); plot(x,M,'LineWidth',1.5); grid on; box on;
ylabel('Moment (N·m)'); title('Bending Moment M');
subplot(5,1,4); plot(x,theta,'LineWidth',1.5); grid on; box on;
ylabel('Slope (rad)'); title('Slope \theta');
subplot(5,1,5); plot(x,delta,'LineWidth',1.5); grid on; box on;
ylabel('Deflection (m)'); xlabel('Height (m)'); title('Deflection \delta');

figure;
plot(x, sigma_b_MPa, 'LineWidth', 1.5); grid on; box on;
xlabel('Height (m)'); ylabel('\sigma_{b,max} (MPa)');
title('Deliverable 5 – Max Bending Stress Along Tower');
end
