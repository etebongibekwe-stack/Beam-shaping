%------------------------------------------------------------
% Plane Wave Expansion method for 2D Photonic Crystal
% Hexagonal (HCP) lattice
% TE polarization
% Based on lecturer's square-lattice code
%------------------------------------------------------------

clear all; clf;

%------------------------------------------------------------
% PHYSICAL PARAMETERS
%------------------------------------------------------------

% Lattice constant
a = 75e-9;              % meters

% Radius of air hole
r = 20e-9;              % meters

% Permittivities
eps1 = 9;               % background dielectric
eps2 = 1;               % air holes (voids)

%------------------------------------------------------------
% NUMERICAL PARAMETERS
%------------------------------------------------------------

precis = 6;             % number of k-points between symmetry points
nG = 4;                 % plane wave expansion order
precisStruct = 40;      % unit cell discretization

%------------------------------------------------------------
% REAL-SPACE LATTICE VECTORS (HEXAGONAL)
%------------------------------------------------------------

a1 = [a, 0];
a2 = [a/2, sqrt(3)*a/2];

% Area of the hexagonal unit cell
Acell = abs(det([a1; a2]));

%------------------------------------------------------------
% DISCRETIZATION OF THE UNIT CELL
%------------------------------------------------------------

nx = 1;
for u = 0:precisStruct
    ny = 1;
    for v = 0:precisStruct

        % Coordinates inside hexagonal unit cell
        x = (u/precisStruct)*a1(1) + (v/precisStruct)*a2(1);
        y = (u/precisStruct)*a1(2) + (v/precisStruct)*a2(2);

        % Center of the hole
        xc = (a1(1)+a2(1))/3;
        yc = (a1(2)+a2(2))/3;

        % Definition of circular air hole
        if sqrt((x-xc)^2 + (y-yc)^2) < r
            struct(nx,ny) = 1/eps2;
        else
            struct(nx,ny) = 1/eps1;
        end

        xSet(nx) = x;
        ySet(ny) = y;

        ny = ny + 1;
    end
    nx = nx + 1;
end

% Mesh element area
dS = Acell / precisStruct^2;

% Mesh grids
xMesh = meshgrid(xSet(1:end-1));
yMesh = meshgrid(ySet(1:end-1))';

% Normalized inverse dielectric function
structMesh = struct(1:end-1,1:end-1) * dS / Acell;

%------------------------------------------------------------
% RECIPROCAL LATTICE VECTORS (HEXAGONAL)
%------------------------------------------------------------

b1 = (2*pi/a)*[1, -1/sqrt(3)];
b2 = (2*pi/a)*[0,  2/sqrt(3)];

numG = 1;
for m = -nG:nG
    for n = -nG:nG
        G(numG,:) = m*b1 + n*b2;
        numG = numG + 1;
    end
end
numG = numG - 1;

%------------------------------------------------------------
% k-PATH IN HEXAGONAL BRILLOUIN ZONE (Γ–M–K–Γ)
%------------------------------------------------------------

GAMMA = [0, 0];
M = [pi/a, pi/(sqrt(3)*a)];
K = [4*pi/(3*a), 0];

kpath = [];

for t = linspace(0,1,precis)
    kpath = [kpath; (1-t)*GAMMA + t*M];
end
for t = linspace(0,1,precis)
    kpath = [kpath; (1-t)*M + t*K];
end
for t = linspace(0,1,precis)
    kpath = [kpath; (1-t)*K + t*GAMMA];
end

kx = kpath(:,1)';
ky = kpath(:,2)';

%------------------------------------------------------------
% FOURIER COEFFICIENTS OF INVERSE DIELECTRIC FUNCTION
%------------------------------------------------------------

for p = 1:numG
    for q = 1:numG
        CN2D_N(p,q) = sum(sum( structMesh .* ...
            exp(1i*((G(p,1)-G(q,1))*xMesh + ...
                    (G(p,2)-G(q,2))*yMesh)) ));
    end
end

%------------------------------------------------------------
% MATRIX OPERATOR FOR TE POLARIZATION
%------------------------------------------------------------

for kidx = 1:length(kx)
    for p = 1:numG
        for q = 1:numG
            Mmat(kidx,p,q) = CN2D_N(p,q) * ...
                ((kx(kidx)+G(p,1))*(kx(kidx)+G(q,1)) + ...
                 (ky(kidx)+G(p,2))*(ky(kidx)+G(q,2)));
        end
    end
end

% ------------------------------------------------------------
% MATRIX OPERATOR FOR TM POLARIZATION
%------------------------------------------------------------

for kidx = 1:length(kx)
    for p = 1:numG
        for q = 1:numG
            if p == q
                Mmat_TM(kidx,p,q) = ...
                    (kx(kidx)+G(p,1))^2 + (ky(kidx)+G(p,2))^2;
            else
                Mmat_TM(kidx,p,q) = 0;
            end
        end
    end
end


%------------------------------------------------------------
% EIGENVALUE PROBLEM AND BAND STRUCTURE
%------------------------------------------------------------

for kidx = 1:length(kx)
    MM = squeeze(Mmat(kidx,:,:));
    [~,V] = eig(MM);
    dispe(:,kidx) = sqrt(real(diag(V))) * a / (2*pi);
end

%------------------------------------------------------------
% TM BAND STRUCTURE COMPUTATION
%------------------------------------------------------------

for kidx = 1:length(kx)
    MM_TM = squeeze(Mmat_TM(kidx,:,:));
    [~,V_TM] = eig(MM_TM);
    dispe_TM(:,kidx) = sqrt(real(diag(V_TM))) * a / (2*pi);
end

%% ---------------- TE bands ----------------
figure(1);
hold on; box on;

nbands = 8;

for n = 1:nbands
    plot(abs(dispe(n,:)),'r','LineWidth',2);

    % Detection of complete band gap
    if n < nbands
        if min(dispe(n+1,:)) > max(dispe(n,:))
            rectangle('Position',[1, max(dispe(n,:)), ...
                length(kx)-1, ...
                min(dispe(n+1,:))-max(dispe(n,:))], ...
                'FaceColor',[0.7 0.7 1], 'EdgeColor','none');
        end
    end
end

set(gca,'XTick',[1 precis 2*precis 3*precis]);
set(gca,'XTickLabel',{'\Gamma','M','K','\Gamma'});
xlabel('Wave vector','FontSize',14);
ylabel('Normalized frequency  \omega a / 2\pi c','FontSize',14);
title('2D Hexagonal Photonic Crystal – TE polarization');
xlim([1 length(kx)]);
grid on;


%% ---------------- TM bands ----------------
figure(2);
hold on; box on;

for n = 1:nbands
    plot(abs(dispe_TM(n,:)),'b','LineWidth',2);
end

set(gca,'XTick',[1 precis 2*precis 3*precis]);
set(gca,'XTickLabel',{'\Gamma','M','K','\Gamma'});
xlabel('Wave vector');
ylabel('\omega a / 2\pi c');
title('TM polarization – Hexagonal lattice');
grid on;
%%
figure(3);
hold on; box on;

for n = 1:6
    plot(abs(dispe(n,:)),'r','LineWidth',2);      % TE
    plot(abs(dispe_TM(n,:)),'b--','LineWidth',2); % TM
end

legend('TE','TM');
title('TE vs TM polarization comparison');


%% ============================================================
% STEP 3 : COMPLETE PHOTONIC BAND GAP MAP vs r/a (TE ONLY)
%============================================================

% Range of r/a values
r_ratio = 0.15 : 0.02 : 0.45;

gap_low  = NaN(size(r_ratio));
gap_high = NaN(size(r_ratio));

nbands_gap = 8;   % number of bands to check for gaps

for ir = 1:length(r_ratio)

    %--------------------------------------------------------
    % Update radius
    %--------------------------------------------------------
    r = r_ratio(ir) * a;

    %--------------------------------------------------------
    % Recompute unit cell (same as before)
    %--------------------------------------------------------
    nx = 1;
    for u = 0:precisStruct
        ny = 1;
        for v = 0:precisStruct

            x = (u/precisStruct)*a1(1) + (v/precisStruct)*a2(1);
            y = (u/precisStruct)*a1(2) + (v/precisStruct)*a2(2);

            xc = (a1(1)+a2(1))/3;
            yc = (a1(2)+a2(2))/3;

            if sqrt((x-xc)^2 + (y-yc)^2) < r
                struct(nx,ny) = 1/eps2;
            else
                struct(nx,ny) = 1/eps1;
            end

            xSet(nx) = x;
            ySet(ny) = y;

            ny = ny + 1;
        end
        nx = nx + 1;
    end

    structMesh = struct(1:end-1,1:end-1) * dS / Acell;

    %--------------------------------------------------------
    % Fourier coefficients
    %--------------------------------------------------------
    for p = 1:numG
        for q = 1:numG
            CN2D_N(p,q) = sum(sum( structMesh .* ...
                exp(1i*((G(p,1)-G(q,1))*xMesh + ...
                        (G(p,2)-G(q,2))*yMesh)) ));
        end
    end

    %--------------------------------------------------------
    % TE band structure
    %--------------------------------------------------------
   for kidx = 1:length(kx)

    % Force Mmat to be 2D (important!)
    Mmat2D = zeros(numG,numG);

    for p = 1:numG
        for q = 1:numG
            Mmat2D(p,q) = CN2D_N(p,q) * ...
                ((kx(kidx)+G(p,1))*(kx(kidx)+G(q,1)) + ...
                 (ky(kidx)+G(p,2))*(ky(kidx)+G(q,2)));
        end
    end

    [~,V] = eig(Mmat2D);
    omega(:,kidx) = sqrt(real(diag(V))) * a / (2*pi);
end


    %--------------------------------------------------------
    % Detection of COMPLETE band gaps
    %--------------------------------------------------------
    for n = 1:nbands_gap-1
        if min(omega(n+1,:)) > max(omega(n,:))
            gap_low(ir)  = max(omega(n,:));
            gap_high(ir) = min(omega(n+1,:));
            break;
        end
    end

end

%------------------------------------------------------------
% PLOT : Complete PBG vs r/a
%------------------------------------------------------------

figure;
hold on; box on;

for i = 1:length(r_ratio)
    if ~isnan(gap_low(i))
        plot([r_ratio(i) r_ratio(i)], ...
             [gap_low(i) gap_high(i)], ...
             'r','LineWidth',4);
    end
end

xlabel('r / a','FontSize',14);
ylabel('Normalized frequency  \omega a / 2\pi c','FontSize',14);
title('Complete TE Photonic Band Gap vs r/a (Hexagonal lattice)');
grid on;
