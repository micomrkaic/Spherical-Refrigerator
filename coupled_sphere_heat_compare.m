% coupled_sphere_heat_compare.m
%
% Compare two spherical cooling models:
%
%   Model 1: water + glass + stagnant-air conduction PDE
%   Model 2: water + glass PDE + Robin convection boundary condition
%
% Geometry:
%   0 < rho < r       water
%   r < rho < r+d     glass
%   r+d < rho < R     stagnant air, model 1 only
%
% Robin model:
%   -k_g dT/drho = h (T_surface - Tf) at rho = r+d
%
% Method:
%   finite-volume discretization
%   generalized eigenproblem K v = lambda C v
%   modal reconstruction theta(t) = sum c_n v_n exp(-lambda_n t)
%
% This version avoids MATLAB-only functions like yline/xline.

clear; clc; close all;

% ============================================================
% USER INPUTS
% ============================================================

r  = 0.050;       % water radius, meters
d  = 0.003;       % glass thickness, meters
R  = 0.150;       % fridge-wall radius for stagnant-air model, meters

T0 = 293.15;      % initial water temperature, K, about 20 C
Tf = 273.15;      % fridge bulk air/wall temperature, K, about 4 C

% Convective heat-transfer coefficient for Robin model
h = 8.0;          % W/m^2/K; try 2, 5, 8, 12, 20

% Material properties
% water
kw   = 0.60;      % W/m/K
rhow = 1000;      % kg/m^3
cw   = 4180;      % J/kg/K

% glass
kg   = 1.05;      % W/m/K
rhog = 2500;      % kg/m^3
cg   = 800;       % J/kg/K

% air
ka   = 0.026;     % W/m/K
rhoa = 1.2;       % kg/m^3
ca   = 1005;      % J/kg/K

% Numerical controls
nt = 500;         % number of output time points

% ============================================================
% DEUTERANOPIA-SAFE COLOR PALETTE
% Okabe-Ito-inspired palette
% ============================================================

col_black   = [0.000, 0.000, 0.000];
col_orange  = [0.902, 0.624, 0.000];
col_sky     = [0.337, 0.706, 0.914];
col_blue    = [0.000, 0.447, 0.698];
col_yellow  = [0.941, 0.894, 0.259];
col_purple  = [0.800, 0.475, 0.655];
col_verm    = [0.835, 0.369, 0.000];
col_teal    = [0.000, 0.620, 0.451];

% Use high-contrast subset heavily:
%   black, blue, orange, purple, teal
palette = [
    col_black;
    col_blue;
    col_orange;
    col_purple;
    col_teal;
    col_sky;
    col_verm;
    col_yellow
];

% ============================================================
% GEOMETRY CHECKS
% ============================================================

b = r + d;
airgap = R - b;

if r <= 0
    error("Need r > 0.");
end

if d <= 0
    error("Need d > 0.");
end

if R <= b
    error("Need R > r+d for the stagnant-air model.");
end

if h < 0
    error("Need h >= 0.");
end

% ============================================================
% MODEL 1: WATER + GLASS + STAGNANT AIR PDE
% ============================================================

fprintf("\n============================================================\n");
fprintf("MODEL 1: water + glass + stagnant-air conduction PDE\n");
fprintf("============================================================\n");

% Separate grids so interfaces are exact.
Nw1 = max(60, ceil(80 * r / R));
Ng1 = max(16, ceil(8 + 8 * d / min(r, airgap)));
Na1 = max(60, ceil(100 * airgap / R));

Nw1 = max(Nw1, 40);
Ng1 = max(Ng1, 12);
Na1 = max(Na1, 40);

edges_w1 = linspace(0, r, Nw1+1);
edges_g1 = linspace(r, b, Ng1+1);
edges_a1 = linspace(b, R, Na1+1);

edges1 = [edges_w1, edges_g1(2:end), edges_a1(2:end)]';

rho1 = 0.5 * (edges1(1:end-1) + edges1(2:end));
dr1  = edges1(2:end) - edges1(1:end-1);

N1 = length(rho1);

k1 = zeros(N1,1);
Cvol1 = zeros(N1,1);

for j = 1:N1
    if rho1(j) < r
        k1(j) = kw;
        Cvol1(j) = rhow * cw;
    elseif rho1(j) < b
        k1(j) = kg;
        Cvol1(j) = rhog * cg;
    else
        k1(j) = ka;
        Cvol1(j) = rhoa * ca;
    end
end

V1 = (4*pi/3) * (edges1(2:end).^3 - edges1(1:end-1).^3);
Cdiag1 = Cvol1 .* V1;

K1 = sparse(N1,N1);

% Internal conductive faces
for j = 1:N1-1
    rf = edges1(j+1);
    Af = 4*pi*rf^2;

    Rleft  = (dr1(j)/2)   / (k1(j)   * Af);
    Rright = (dr1(j+1)/2) / (k1(j+1) * Af);

    G = 1 / (Rleft + Rright);

    K1(j,j)       = K1(j,j)       + G;
    K1(j+1,j+1)   = K1(j+1,j+1)   + G;
    K1(j,j+1)     = K1(j,j+1)     - G;
    K1(j+1,j)     = K1(j+1,j)     - G;
end

% Outer Dirichlet wall at rho = R: theta(R,t)=0
rf = R;
Af = 4*pi*rf^2;
Gwall = 1 / ((dr1(N1)/2) / (k1(N1)*Af));
K1(N1,N1) = K1(N1,N1) + Gwall;

Cmat1 = spdiags(Cdiag1, 0, N1, N1);

nmodes1 = min(N1-2, max(50, ceil(0.45 * N1)));

fprintf("Grid cells: water %d, glass %d, air %d, total %d\n", ...
        Nw1, Ng1, Na1, N1);
fprintf("Using %d eigenmodes\n", nmodes1);

opts.tol = 1e-10;
opts.maxit = 3000;
opts.disp = 0;

[Vmode1, Lambda1] = eigs(K1, Cmat1, nmodes1, "sm", opts);

lambda1 = diag(Lambda1);
[lambda1, idx1] = sort(lambda1);
Vmode1 = Vmode1(:,idx1);

keep1 = find(lambda1 > 0);
lambda1 = lambda1(keep1);
Vmode1 = Vmode1(:,keep1);
nmodes1 = length(lambda1);

% C-orthonormalize
for n = 1:nmodes1
    normC = sqrt(Vmode1(:,n)' * (Cdiag1 .* Vmode1(:,n)));
    Vmode1(:,n) = Vmode1(:,n) / normC;
end

theta01 = zeros(N1,1);
theta01(rho1 < r) = T0 - Tf;

coef1 = zeros(nmodes1,1);
for n = 1:nmodes1
    coef1(n) = Vmode1(:,n)' * (Cdiag1 .* theta01);
end

tau1_stag = 1 / lambda1(1);

fprintf("Dominant eigenvalue = %.6e 1/s\n", lambda1(1));
fprintf("Dominant time scale = %.3f hours\n", tau1_stag/3600);

% ============================================================
% MODEL 2: WATER + GLASS + ROBIN CONVECTION
% ============================================================

fprintf("\n============================================================\n");
fprintf("MODEL 2: water + glass + Robin convection boundary\n");
fprintf("============================================================\n");

% Domain only goes to b = r+d.
% Need good glass resolution because thin glass controls boundary flux.
Nw2 = max(80, ceil(120 * r / b));
Ng2 = max(24, ceil(12 + 20 * d / r));

edges_w2 = linspace(0, r, Nw2+1);
edges_g2 = linspace(r, b, Ng2+1);

edges2 = [edges_w2, edges_g2(2:end)]';

rho2 = 0.5 * (edges2(1:end-1) + edges2(2:end));
dr2  = edges2(2:end) - edges2(1:end-1);

N2 = length(rho2);

k2 = zeros(N2,1);
Cvol2 = zeros(N2,1);

for j = 1:N2
    if rho2(j) < r
        k2(j) = kw;
        Cvol2(j) = rhow * cw;
    else
        k2(j) = kg;
        Cvol2(j) = rhog * cg;
    end
end

V2 = (4*pi/3) * (edges2(2:end).^3 - edges2(1:end-1).^3);
Cdiag2 = Cvol2 .* V2;

K2 = sparse(N2,N2);

% Internal conductive faces
for j = 1:N2-1
    rf = edges2(j+1);
    Af = 4*pi*rf^2;

    Rleft  = (dr2(j)/2)   / (k2(j)   * Af);
    Rright = (dr2(j+1)/2) / (k2(j+1) * Af);

    G = 1 / (Rleft + Rright);

    K2(j,j)       = K2(j,j)       + G;
    K2(j+1,j+1)   = K2(j+1,j+1)   + G;
    K2(j,j+1)     = K2(j,j+1)     - G;
    K2(j+1,j)     = K2(j+1,j)     - G;
end

% Robin convection boundary at outer glass surface rho = b:
%
%   -k_g dtheta/drho = h theta
%
% Total heat loss:
%
%   Q = h A_b theta_surface
%
% Finite-volume implementation:
% last cell connected to bulk air theta = 0 by:
%
%   resistance = conduction from cell center to surface + convection resistance
%
%   R_cond = (dr_last/2)/(k_g A_b)
%   R_conv = 1/(h A_b)
%
%   G_robin = 1/(R_cond + R_conv)
%
% If h = 0, G_robin = 0, insulated exterior.

Ab = 4*pi*b^2;

if h == 0
    Grobin = 0;
else
    Rcond_surf = (dr2(N2)/2) / (k2(N2) * Ab);
    Rconv = 1 / (h * Ab);
    Grobin = 1 / (Rcond_surf + Rconv);
end

K2(N2,N2) = K2(N2,N2) + Grobin;

Cmat2 = spdiags(Cdiag2, 0, N2, N2);

nmodes2 = min(N2-2, max(50, ceil(0.60 * N2)));

fprintf("Grid cells: water %d, glass %d, total %d\n", Nw2, Ng2, N2);
fprintf("Using %d eigenmodes\n", nmodes2);
fprintf("Robin h = %.3f W/m^2/K\n", h);
fprintf("Robin boundary conductance = %.6e W/K\n", Grobin);

[Vmode2, Lambda2] = eigs(K2, Cmat2, nmodes2, "sm", opts);

lambda2 = diag(Lambda2);
[lambda2, idx2] = sort(lambda2);
Vmode2 = Vmode2(:,idx2);

keep2 = find(lambda2 > 0);
lambda2 = lambda2(keep2);
Vmode2 = Vmode2(:,keep2);
nmodes2 = length(lambda2);

% C-orthonormalize
for n = 1:nmodes2
    normC = sqrt(Vmode2(:,n)' * (Cdiag2 .* Vmode2(:,n)));
    Vmode2(:,n) = Vmode2(:,n) / normC;
end

theta02 = zeros(N2,1);
theta02(rho2 < r) = T0 - Tf;

coef2 = zeros(nmodes2,1);
for n = 1:nmodes2
    coef2(n) = Vmode2(:,n)' * (Cdiag2 .* theta02);
end

tau1_robin = 1 / lambda2(1);

fprintf("Dominant eigenvalue = %.6e 1/s\n", lambda2(1));
fprintf("Dominant time scale = %.3f hours\n", tau1_robin/3600);

% ============================================================
% COMMON TIME GRID FOR COMPARISON
% ============================================================

% Use enough time to show both behaviors.
tmax = max(5*tau1_robin, 2.5*tau1_stag);

% Reasonable display bounds
tmax = max(tmax, 4*3600);       % at least 4 hours
tmax = min(tmax, 96*3600);      % at most 96 hours

time = linspace(0, tmax, nt)';

fprintf("\n============================================================\n");
fprintf("COMPARISON TIME HORIZON\n");
fprintf("============================================================\n");
fprintf("Final plotted time = %.3f hours\n", tmax/3600);

% ============================================================
% RECONSTRUCT BOTH SOLUTIONS
% ============================================================

water_cells1 = find(rho1 < r);
glass_cells1 = find(rho1 >= r & rho1 < b);
air_cells1   = find(rho1 >= b);

water_cells2 = find(rho2 < r);
glass_cells2 = find(rho2 >= r);

Vwater1 = sum(V1(water_cells1));
Vglass1 = sum(V1(glass_cells1));
Vair1   = sum(V1(air_cells1));

Vwater2 = sum(V2(water_cells2));
Vglass2 = sum(V2(glass_cells2));

Twatermean1 = zeros(nt,1);
Tcenter1    = zeros(nt,1);
Tsurface1   = zeros(nt,1);
Tglassmean1 = zeros(nt,1);
Tairmean1   = zeros(nt,1);

Twatermean2 = zeros(nt,1);
Tcenter2    = zeros(nt,1);
Tsurface2   = zeros(nt,1);
Tglassmean2 = zeros(nt,1);

% Surface cell approximations
surf_idx1 = max(find(rho1 < b));
surf_idx2 = N2;

for it = 1:nt
    t = time(it);

    decay1 = exp(-lambda1 * t);
    theta1 = Vmode1 * (coef1 .* decay1);
    T1 = Tf + theta1;

    decay2 = exp(-lambda2 * t);
    theta2 = Vmode2 * (coef2 .* decay2);
    T2 = Tf + theta2;

    Tcenter1(it) = T1(1);
    Tsurface1(it) = T1(surf_idx1);
    Twatermean1(it) = sum(T1(water_cells1) .* V1(water_cells1)) / Vwater1;
    Tglassmean1(it) = sum(T1(glass_cells1) .* V1(glass_cells1)) / Vglass1;
    Tairmean1(it)   = sum(T1(air_cells1)   .* V1(air_cells1))   / Vair1;

    Tcenter2(it) = T2(1);
    Tsurface2(it) = T2(surf_idx2);
    Twatermean2(it) = sum(T2(water_cells2) .* V2(water_cells2)) / Vwater2;
    Tglassmean2(it) = sum(T2(glass_cells2) .* V2(glass_cells2)) / Vglass2;
end

% Selected profile times based on Robin cooling, but within common horizon
profile_times = [0, 0.25, 0.50, 1.00, 2.00, 4.00] * tau1_robin;
profile_times = profile_times(profile_times <= tmax);

if length(profile_times) < 3
    profile_times = [0, 0.25, 0.50, 1.00] * tmax;
end

profiles1 = zeros(N1, length(profile_times));
profiles2 = zeros(N2, length(profile_times));

for ip = 1:length(profile_times)
    t = profile_times(ip);

    decay1 = exp(-lambda1 * t);
    theta1 = Vmode1 * (coef1 .* decay1);
    profiles1(:,ip) = Tf + theta1;

    decay2 = exp(-lambda2 * t);
    theta2 = Vmode2 * (coef2 .* decay2);
    profiles2(:,ip) = Tf + theta2;
end

% ============================================================
% PRINT DIAGNOSTIC SUMMARY
% ============================================================

fprintf("\n============================================================\n");
fprintf("SUMMARY\n");
fprintf("============================================================\n");

fprintf("Stagnant-air model tau_1: %.3f hours\n", tau1_stag/3600);
fprintf("Robin model tau_1:        %.3f hours\n", tau1_robin/3600);
fprintf("Ratio stagnant/Robin:     %.3f\n", tau1_stag/tau1_robin);

fprintf("\nFirst few stagnant-air eigenvalues:\n");
for n = 1:min(6,nmodes1)
    fprintf("  n=%3d lambda=%.6e tau=%.4f hours\n", ...
            n, lambda1(n), 1/lambda1(n)/3600);
end

fprintf("\nFirst few Robin eigenvalues:\n");
for n = 1:min(6,nmodes2)
    fprintf("  n=%3d lambda=%.6e tau=%.4f hours\n", ...
            n, lambda2(n), 1/lambda2(n)/3600);
end

% Initial reconstruction checks
theta01_recon = Vmode1 * coef1;
theta02_recon = Vmode2 * coef2;

err1 = norm(theta01_recon - theta01) / max(norm(theta01), eps);
err2 = norm(theta02_recon - theta02) / max(norm(theta02), eps);

fprintf("\nInitial condition reconstruction error:\n");
fprintf("  stagnant-air model: %.6e\n", err1);
fprintf("  Robin model:        %.6e\n", err2);

% ============================================================
% PLOT 1: MAIN COMPARISON, WATER MEAN AND CENTER
% ============================================================

figure;
hold on;
set(gca, "colororder", palette);

plot(time/3600, Twatermean1, "-",  "color", col_blue,   "linewidth", 2.4);
plot(time/3600, Twatermean2, "-",  "color", col_orange, "linewidth", 2.4);

plot(time/3600, Tcenter1,    "--", "color", col_blue,   "linewidth", 1.8);
plot(time/3600, Tcenter2,    "--", "color", col_orange, "linewidth", 1.8);

xl = xlim();
plot(xl, [Tf Tf], ":", "color", col_black, "linewidth", 1.4);

xlabel("time, hours");
ylabel("temperature, K");
title("Cooling comparison: stagnant air vs Robin convection");
legend("mean water: stagnant air PDE", ...
       sprintf("mean water: Robin h = %.1f", h), ...
       "center water: stagnant air PDE", ...
       sprintf("center water: Robin h = %.1f", h), ...
       "Tf", ...
       "location", "northeast");
grid on;

% ============================================================
% PLOT 2: SURFACE AND MEAN GLASS TEMPERATURES
% ============================================================

figure;
hold on;
set(gca, "colororder", palette);

plot(time/3600, Tsurface1,   "-",  "color", col_blue,   "linewidth", 2.2);
plot(time/3600, Tsurface2,   "-",  "color", col_orange, "linewidth", 2.2);

plot(time/3600, Tglassmean1, "--", "color", col_blue,   "linewidth", 1.7);
plot(time/3600, Tglassmean2, "--", "color", col_orange, "linewidth", 1.7);

xl = xlim();
plot(xl, [Tf Tf], ":", "color", col_black, "linewidth", 1.4);

xlabel("time, hours");
ylabel("temperature, K");
title("Glass/surface temperatures");
legend("outer glass: stagnant air PDE", ...
       sprintf("outer glass: Robin h = %.1f", h), ...
       "mean glass: stagnant air PDE", ...
       sprintf("mean glass: Robin h = %.1f", h), ...
       "Tf", ...
       "location", "northeast");
grid on;

% ============================================================
% PLOT 3: RADIAL PROFILES, STAGNANT-AIR MODEL
% ============================================================

figure;
hold on;

profile_cols = [col_black; col_blue; col_orange; col_purple; col_teal; col_verm];

for ip = 1:length(profile_times)
    cidx = 1 + mod(ip-1, rows(profile_cols));
    plot(rho1, profiles1(:,ip), "-", ...
         "color", profile_cols(cidx,:), "linewidth", 1.8);
end

yl = ylim();
plot([r r], yl, "--", "color", col_black, "linewidth", 1.2);
plot([b b], yl, "--", "color", col_black, "linewidth", 1.2);

xl = xlim();
plot(xl, [Tf Tf], ":", "color", col_black, "linewidth", 1.2);

xlabel("radius rho, meters");
ylabel("temperature, K");
title("Radial profiles: stagnant-air PDE model");
grid on;

leg = cell(length(profile_times),1);
for ip = 1:length(profile_times)
    leg{ip} = sprintf("t = %.2f h", profile_times(ip)/3600);
end
legend(leg, "location", "northeast");

% ============================================================
% PLOT 4: RADIAL PROFILES, ROBIN MODEL
% ============================================================

figure;
hold on;

for ip = 1:length(profile_times)
    cidx = 1 + mod(ip-1, rows(profile_cols));
    plot(rho2, profiles2(:,ip), "-", ...
         "color", profile_cols(cidx,:), "linewidth", 1.8);
end

yl = ylim();
plot([r r], yl, "--", "color", col_black, "linewidth", 1.2);
plot([b b], yl, "--", "color", col_black, "linewidth", 1.2);

xl = xlim();
plot(xl, [Tf Tf], ":", "color", col_black, "linewidth", 1.2);

xlabel("radius rho, meters");
ylabel("temperature, K");
title(sprintf("Radial profiles: Robin convection model, h = %.1f W/m^2/K", h));
grid on;

leg = cell(length(profile_times),1);
for ip = 1:length(profile_times)
    leg{ip} = sprintf("t = %.2f h", profile_times(ip)/3600);
end
legend(leg, "location", "northeast");

% ============================================================
% PLOT 5: FIRST FEW EIGENFUNCTIONS, BOTH MODELS
% ============================================================

figure;
hold on;

nplot = min(5, min(nmodes1, nmodes2));

for n = 1:nplot
    cidx = 1 + mod(n-1, rows(profile_cols));

    v1 = Vmode1(:,n);
    v1 = v1 / max(abs(v1));

    v2 = Vmode2(:,n);
    v2 = v2 / max(abs(v2));

    plot(rho1, v1, "-",  "color", profile_cols(cidx,:), "linewidth", 1.5);
    plot(rho2, v2, "--", "color", profile_cols(cidx,:), "linewidth", 1.8);
end

yl = ylim();
plot([r r], yl, ":", "color", col_black, "linewidth", 1.2);
plot([b b], yl, ":", "color", col_black, "linewidth", 1.2);

xlabel("radius rho, meters");
ylabel("scaled eigenfunction");
title("First few eigenfunctions: solid = stagnant air, dashed = Robin");
grid on;

leg = cell(2*nplot,1);
idx = 1;
for n = 1:nplot
    leg{idx} = sprintf("stagnant mode %d", n);
    idx = idx + 1;
    leg{idx} = sprintf("Robin mode %d", n);
    idx = idx + 1;
end
legend(leg, "location", "northeast");

% ============================================================
% OPTIONAL SENSITIVITY: ROBIN h VALUES
% ============================================================
%
% This section uses a cheap implicit time-step approximation by reusing the
% same water+glass grid and recomputing only the Robin boundary conductance.
% It compares mean water temperature for several h values using eigenmodes.
%
% Set do_h_sensitivity = false if you do not want this plot.

do_h_sensitivity = true;

if do_h_sensitivity

    h_values = [2, 5, 8, 12, 20];

    figure;
    hold on;

    for ih = 1:length(h_values)

        hh = h_values(ih);

        % Rebuild K for water+glass interior using same grid.
        Kh = sparse(N2,N2);

        for j = 1:N2-1
            rf = edges2(j+1);
            Af = 4*pi*rf^2;

            Rleft  = (dr2(j)/2)   / (k2(j)   * Af);
            Rright = (dr2(j+1)/2) / (k2(j+1) * Af);

            G = 1 / (Rleft + Rright);

            Kh(j,j)       = Kh(j,j)       + G;
            Kh(j+1,j+1)   = Kh(j+1,j+1)   + G;
            Kh(j,j+1)     = Kh(j,j+1)     - G;
            Kh(j+1,j)     = Kh(j+1,j)     - G;
        end

        if hh == 0
            Gh = 0;
        else
            Rcond_surf = (dr2(N2)/2) / (k2(N2) * Ab);
            Rconv = 1 / (hh * Ab);
            Gh = 1 / (Rcond_surf + Rconv);
        end

        Kh(N2,N2) = Kh(N2,N2) + Gh;

        nmh = nmodes2;

        [Vh, Lh] = eigs(Kh, Cmat2, nmh, "sm", opts);

        lambdah = diag(Lh);
        [lambdah, idxh] = sort(lambdah);
        Vh = Vh(:,idxh);

        keeph = find(lambdah > 0);
        lambdah = lambdah(keeph);
        Vh = Vh(:,keeph);
        nmh = length(lambdah);

        for n = 1:nmh
            normC = sqrt(Vh(:,n)' * (Cdiag2 .* Vh(:,n)));
            Vh(:,n) = Vh(:,n) / normC;
        end

        coefh = zeros(nmh,1);
        for n = 1:nmh
            coefh(n) = Vh(:,n)' * (Cdiag2 .* theta02);
        end

        Tmean_h = zeros(nt,1);

        for it = 1:nt
            t = time(it);
            decayh = exp(-lambdah * t);
            thetah = Vh * (coefh .* decayh);
            Th = Tf + thetah;

            Tmean_h(it) = sum(Th(water_cells2) .* V2(water_cells2)) / Vwater2;
        end

        cidx = 1 + mod(ih-1, rows(profile_cols));

        if ih == 1
            ls = "-";
        elseif ih == 2
            ls = "--";
        elseif ih == 3
            ls = "-.";
        elseif ih == 4
            ls = ":";
        else
            ls = "-";
        end

        plot(time/3600, Tmean_h, ls, ...
             "color", profile_cols(cidx,:), "linewidth", 2.2);
    end

    xl = xlim();
    plot(xl, [Tf Tf], ":", "color", col_black, "linewidth", 1.3);

    xlabel("time, hours");
    ylabel("mean water temperature, K");
    title("Robin model sensitivity to convective heat-transfer coefficient h");
    grid on;

    leg = cell(length(h_values)+1,1);
    for ih = 1:length(h_values)
        leg{ih} = sprintf("h = %.0f W/m^2/K", h_values(ih));
    end
    leg{end} = "Tf";
    legend(leg, "location", "northeast");
end

% ============================================================
% COMBINED SUMMARY FIGURE WITH SUBPLOTS
% Same y-axis limits in all four panels
% ============================================================

% Compute common y-axis limits across all plotted temperature series/profiles.
all_yvals = [
    Twatermean1(:);
    Twatermean2(:);
    Tcenter1(:);
    Tcenter2(:);
    Tsurface1(:);
    Tsurface2(:);
    Tglassmean1(:);
    Tglassmean2(:);
    profiles1(:);
    profiles2(:);
    Tf
];

ymin = min(all_yvals);
ymax = max(all_yvals);

% Add small padding.
ypad = 0.05 * (ymax - ymin);
if ypad == 0
    ypad = 1.0;
end

common_ylim = [ymin - ypad, ymax + ypad];

figure;

% ----------------------------
% Subplot 1: water cooling
% ----------------------------
subplot(2,2,1);
hold on;

plot(time/3600, Twatermean1, "-",  "color", col_blue,   "linewidth", 2.2);
plot(time/3600, Twatermean2, "-",  "color", col_orange, "linewidth", 2.2);

plot(time/3600, Tcenter1,    "--", "color", col_blue,   "linewidth", 1.6);
plot(time/3600, Tcenter2,    "--", "color", col_orange, "linewidth", 1.6);

ylim(common_ylim);
xl = xlim();
plot(xl, [Tf Tf], ":", "color", col_black, "linewidth", 1.2);

xlabel("time, hours");
ylabel("temperature, K");
title("Water cooling");
legend("mean stagnant", ...
       "mean Robin", ...
       "center stagnant", ...
       "center Robin", ...
       "Tf", ...
       "location", "northeast");
grid on;

% ----------------------------
% Subplot 2: glass/surface temperatures
% ----------------------------
subplot(2,2,2);
hold on;

plot(time/3600, Tsurface1,   "-",  "color", col_blue,   "linewidth", 2.0);
plot(time/3600, Tsurface2,   "-",  "color", col_orange, "linewidth", 2.0);

plot(time/3600, Tglassmean1, "--", "color", col_blue,   "linewidth", 1.5);
plot(time/3600, Tglassmean2, "--", "color", col_orange, "linewidth", 1.5);

ylim(common_ylim);
xl = xlim();
plot(xl, [Tf Tf], ":", "color", col_black, "linewidth", 1.2);

xlabel("time, hours");
ylabel("temperature, K");
title("Glass and surface");
legend("surface stagnant", ...
       "surface Robin", ...
       "mean glass stagnant", ...
       "mean glass Robin", ...
       "Tf", ...
       "location", "northeast");
grid on;

% ----------------------------
% Subplot 3: stagnant-air radial profiles
% ----------------------------
subplot(2,2,3);
hold on;

for ip = 1:length(profile_times)
    cidx = 1 + mod(ip-1, rows(profile_cols));
    plot(rho1, profiles1(:,ip), "-", ...
         "color", profile_cols(cidx,:), "linewidth", 1.6);
end

ylim(common_ylim);
yl = ylim();
plot([r r], yl, "--", "color", col_black, "linewidth", 1.0);
plot([b b], yl, "--", "color", col_black, "linewidth", 1.0);

xl = xlim();
plot(xl, [Tf Tf], ":", "color", col_black, "linewidth", 1.0);

xlabel("radius rho, meters");
ylabel("temperature, K");
title("Profiles: stagnant-air PDE");
grid on;

leg3 = cell(length(profile_times),1);
for ip = 1:length(profile_times)
    leg3{ip} = sprintf("%.2f h", profile_times(ip)/3600);
end
legend(leg3, "location", "northeast");

% ----------------------------
% Subplot 4: Robin radial profiles
% ----------------------------
subplot(2,2,4);
hold on;

for ip = 1:length(profile_times)
    cidx = 1 + mod(ip-1, rows(profile_cols));
    plot(rho2, profiles2(:,ip), "-", ...
         "color", profile_cols(cidx,:), "linewidth", 1.6);
end

ylim(common_ylim);
yl = ylim();
plot([r r], yl, "--", "color", col_black, "linewidth", 1.0);
plot([b b], yl, "--", "color", col_black, "linewidth", 1.0);

xl = xlim();
plot(xl, [Tf Tf], ":", "color", col_black, "linewidth", 1.0);

xlabel("radius rho, meters");
ylabel("temperature, K");
title(sprintf("Profiles: Robin h = %.1f", h));
grid on;

leg4 = cell(length(profile_times),1);
for ip = 1:length(profile_times)
    leg4{ip} = sprintf("%.2f h", profile_times(ip)/3600);
end
legend(leg4, "location", "northeast");

% ----------------------------
% Overall figure title if available
% ----------------------------

if exist("sgtitle", "file") || exist("sgtitle", "builtin")
    sgtitle("Coupled spherical cooling: stagnant air vs Robin convection");
else
    annotation("textbox", [0 0.95 1 0.05], ...
               "string", "Coupled spherical cooling: stagnant air vs Robin convection", ...
               "edgecolor", "none", ...
               "horizontalalignment", "center", ...
               "fontsize", 14, ...
               "fontweight", "bold");
end
