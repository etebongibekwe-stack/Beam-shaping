%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  BPM_Full_Gaussian_Lens_Cone_fixed_STATIC_COLORBAR.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% ========================================================================
% PHYSICAL PARAMETERS
% ========================================================================
lambda = 1.03;      
k0 = 2*pi/lambda;

n_air = 1.00;
n_Si  = 3.5651;

omega0 = 50;  % This beam diameter is maintained. Changing it affect the nature of the observed rings
f = 10000;          

%% ========================================================================
% GEOMETRY PARAMETERS
% ========================================================================
z_lens   = 9990; % propagation distance after the lens
z_Si     = 10;   % propagation distance in the slilcon before the cone. The sum of the distances matches the focal length of the lens
h_cone   = 10;
well_D   = 10;
cone_D0  = 5;
cone_D1  = 1;

z_final_above_cone = 200;

%% ========================================================================
% SIMULATION GRID
% ========================================================================
Nx = 512; Ny = 512;
Lx = 200; Ly = Lx;

dx = Lx/Nx; dy = Ly/Ny;
x = (-Nx/2:Nx/2-1)*dx;
y = (-Ny/2:Ny/2-1)*dy;
[X,Y] = meshgrid(x,y);

Zmax = z_Si + h_cone + z_final_above_cone;
dz = 0.5;
z = 0:dz:Zmax;
Nz = numel(z);

%% ========================================================================
% SPECTRAL GRID
% ========================================================================
fx = (-Nx/2:Nx/2-1)/Lx;
fy = (-Ny/2:Ny/2-1)/Ly;
[KX,KY] = meshgrid(2*pi*fx,2*pi*fy);
H_half = exp(1i*(KX.^2+KY.^2)/(2*k0*n_air)*(dz/2));

%% ========================================================================
% INITIAL FIELD
% ========================================================================
E = exp(-(X.^2+Y.^2)/omega0^2);
E = E / sqrt(sum(abs(E(:)).^2)*dx*dy);

%% ========================================================================
% LENS
% ========================================================================
t_lens = exp(1i*k0/(2*f)*(X.^2+Y.^2));
[~,iz_lens] = min(abs(z-z_lens));

%% ========================================================================
% REFRACTIVE INDEX PARAMETERS
% ========================================================================
par.n_air = n_air;
par.n_Si  = n_Si;
par.z_Si  = z_Si;
par.h_cone = h_cone;
par.well_D = well_D;
par.cone_D0 = cone_D0;
par.cone_D1 = cone_D1;
par.cx = 0; par.cy = 0;

%% ========================================================================
% STORAGE
% ========================================================================
cy = round(Ny/2);
I_xz = zeros(Nx,Nz);

save_every = max(1,round(Nz/200));
I_stack = zeros(Ny,Nx,ceil(Nz/save_every),'single');
z_store = zeros(1,size(I_stack,3));
frame_idx = 0;

%z_targets = [z_Si-1, z_Si+1, z_Si+h_cone/2, z_Si+h_cone+1];
z_targets = [25 50 75 100 125 150 175 200];
E_at_plane = cell(size(z_targets));
plane_found = false(size(z_targets));

%% ========================================================================
% PROPAGATION
% ========================================================================
for iz = 1:Nz
    zc = z(iz);

    Ek = fftshift(fft2(ifftshift(E)));
    Ek = Ek .* H_half;
    E  = fftshift(ifft2(ifftshift(Ek)));

    nxy = index_map(X,Y,zc,par);
    E = E .* exp(1i*k0*(nxy-n_air)*dz);

    if iz==iz_lens
        E = E .* t_lens;
    end

    Ek = fftshift(fft2(ifftshift(E)));
    Ek = Ek .* H_half;
    E  = fftshift(ifft2(ifftshift(Ek)));

    I_xz(:,iz) = abs(E(cy,:)).^2;

    if mod(iz,save_every)==0 || iz==1 || iz==Nz
        frame_idx = frame_idx+1;
        I_stack(:,:,frame_idx) = abs(E).^2;
        z_store(frame_idx) = zc;
    end

    for m=1:length(z_targets)
        if ~plane_found(m) && abs(zc-z_targets(m))<dz/2
            E_at_plane{m} = E;
            plane_found(m)=true;
        end
    end
end

%% ========================================================================
% GLOBAL COLOR SCALE
% ========================================================================
Imin = 0;
Imax = prctile(I_stack(:),99);
%% ========================================================================
% 2D INTENSITY MAPS AT FIXED z AFTER CONE
% ========================================================================

z_tip = z_Si + h_cone;
z_after = z_tip + [25 50 75 100 125 150 200];

E_after_planes = cell(size(z_after));
plane_found = false(size(z_after));

% --- Extract closest stored planes ---
for iz = 1:length(z_store)
    zc = z_store(iz);
    for m = 1:length(z_after)
        if ~plane_found(m) && abs(zc - z_after(m)) < dz
            E_after_planes{m} = I_stack(:,:,iz);
            plane_found(m) = true;
        end
    end
end


%% ========================================================================
% 2D INTENSITY AT SELECTED Z-PLANES  ✅
% ========================================================================
figure('Name','2D Intensity at Selected z-Planes');
for m=1:length(z_targets)
    subplot(2,4,m)
    imagesc(x,y,abs(E_at_plane{m}).^2)
    axis image; set(gca,'YDir','normal');
    caxis([Imin Imax]);
    colorbar;
    cb = colorbar;
    cb.Label.String = 'Intensity (arb. units)';
    ylabel('x (\mum)')
    xlabel('y (\mum)')
    title(sprintf('z = %.1f µm',z_targets(m)))
end

%% ========================================================================
% LONGITUDINAL x-z INTENSITY
% ========================================================================
figure('Name','Longitudinal x-z Intensity');
imagesc(z,x,I_xz)
set(gca,'YDir','normal')
caxis([Imin Imax])
cb = colorbar
cb = colorbar;
cb.Label.String = 'Intensity (arb. units)';
xlabel('z (µm)'); ylabel('x (µm)')
title('Longitudinal Intensity')

%% ========================================================================
% FINAL PLANE
% ========================================================================
figure('Name','Final Intensity');
imagesc(x,y,I_stack(:,:,end))
axis image; set(gca,'YDir','normal')
caxis([Imin Imax])
colorbar
cb = colorbar;
cb.Label.String = 'Intensity (arb. units)';
ylabel('x (\mum)')
xlabel('y (\mum)')
title(sprintf('z = %.1f µm',z(end)))

%% ========================================================================
% ANIMATION
% ========================================================================
figure('Name','Propagation Animation');
for k=1:size(I_stack,3)
    imagesc(x,y,I_stack(:,:,k))
    axis image; set(gca,'YDir','normal')
    caxis([Imin Imax])
    colorbar
    cb = colorbar;
    cb.Label.String = 'Intensity (arb. units)';
    ylabel('x (\mum)')
    xlabel('y (\mum)')
    ylabel('x(µm)')
    title(sprintf('z = %.1f µm',z_store(k)))
    drawnow
end

%% ========================================================================
% VIDEO EXPORT
% ========================================================================
v = VideoWriter('beam_propagation.mp4','MPEG-4');
v.FrameRate = 20;
open(v)

figure;
for k=1:size(I_stack,3);
    imagesc(x,y,I_stack(:,:,k))
    axis image; set(gca,'YDir','normal')
    caxis([Imin Imax])
    colorbar
    cb = colorbar;
    cb.Label.String = 'Intensity (arb. units)';
    ylabel('x (\mum)')
    xlabel('y (\mum)')
    title(sprintf('z = %.1f µm',z_store(k)))
    writeVideo(v,getframe(gcf))
end
close(v)

fprintf('Video saved with static colorbar.\n')

%% ========================================================================
% INDEX MAP FUNCTION
% ========================================================================
function nxy = index_map(X,Y,zc,par)

nxy = par.n_air*ones(size(X));

if zc < par.z_Si
    nxy(:) = par.n_Si;
    return
end

if zc <= par.z_Si+par.h_cone
    Rwell = par.well_D/2;
    h_rel = (zc-par.z_Si)/par.h_cone;
    Rcone = (par.cone_D0*(1-h_rel)+par.cone_D1*h_rel)/2;

    r = sqrt((X-par.cx).^2+(Y-par.cy).^2);
    nxy(r<=Rwell) = par.n_Si;
    nxy(r<=Rcone) = par.n_Si;
end
end
