% =========================================================================
% --- FDTD 2D Nanophotonic Project: FULL PHOTONIC NANOJET + REFLECTION ---
% --- Adaptive visualization + FULL analysis restored
% =========================================================================
clear; close all; clc;

%% =========================================================================
% 0. CONSTANTS AND GRID
% =========================================================================
c0   = 3e8;
eps0 = 8.854187817e-12;
mu0  = 4*pi*1e-7;
eta0 = sqrt(mu0/eps0);

lambda = 800e-9;
omega  = 2*pi*c0/lambda;

dx = lambda/30;
dy = dx;
Nx = 610; Ny = 310;
S  = 1/sqrt(2);
dt = S*dx/c0;
Niter = 1200;
npml = 70;

%% =========================================================================
% 1. SOURCE
% =========================================================================
tau = 10e-15;
t0  = 4*tau;
src_x = round(0.08*Nx);
src_amp = 0.2

wy = 5e-6;
ypos = ((1:Ny)-0.5*Ny)*dy;
spatial_profile = src_amp*exp(-(ypos/wy).^2);

%% =========================================================================
% 2. MATERIAL
% =========================================================================
n_bg = 1.0;
sigma_particle = 10;

R_values = [1 2 2.5 3]*1e-6;
n_values = [1.6];

%% =========================================================================
% 3. SCENARIOS
% =========================================================================
Scenarios(1).TwoSpheres = 0;
Scenarios(1).Title = 'Single Sphere';
Scenarios(2).TwoSpheres = 1;
Scenarios(2).Title = 'Two Interacting Spheres';

%% =========================================================================
% 4. PML
% =========================================================================
m = 5;
sigma_max = (m+1)/(npml*dx*eta0);
sigma_x = zeros(Nx,1); sigma_y = zeros(1,Ny);
for i=1:npml
    s = sigma_max*((npml-i+0.5)/npml)^m;
    sigma_x(i)=s; sigma_x(Nx-i+1)=s;
    sigma_y(i)=s; sigma_y(Ny-i+1)=s;
end

sigma_e_pml = repmat(sigma_x,1,Ny) + repmat(sigma_y,Nx,1);
sigma_m_x = (mu0/eps0)*repmat(sigma_x,1,Ny);
sigma_m_y = (mu0/eps0)*repmat(sigma_y,Nx,1);

aHx = (1 - dt*sigma_m_y/(2*mu0)) ./ (1 + dt*sigma_m_y/(2*mu0));
bHx = (dt/mu0) ./ (1 + dt*sigma_m_y/(2*mu0));
aHy = (1 - dt*sigma_m_x/(2*mu0)) ./ (1 + dt*sigma_m_x/(2*mu0));
bHy = (dt/mu0) ./ (1 + dt*sigma_m_x/(2*mu0));

%% =========================================================================
% 5. MAIN LOOP
% =========================================================================
snap_interval = 60;
vid_fps = 5;

for run_idx = 1:length(Scenarios)

    TwoSpheres = Scenarios(run_idx).TwoSpheres;
    run_title  = Scenarios(run_idx).Title;

    % Summary arrays
    FEF_all = zeros(length(R_values),1);
    FWHM_all = zeros(length(R_values),1);
    NanojetL_all = zeros(length(R_values),1);
    Joule_all = zeros(length(R_values),1);
    Reflection_all = zeros(length(R_values),1);

    for r_idx = 1:length(R_values)

        R = R_values(r_idx);

        %% --- Field initialization
        Ez = zeros(Nx,Ny); Hx = zeros(Nx,Ny); Hy = zeros(Nx,Ny);
        eps_r = ones(Nx,Ny)*n_bg^2;
        sigma_mat = zeros(Nx,Ny);

        cx = round(0.45*Nx); 
        cy = round(0.5*Ny);
        gap_cells = round(1e-6/dx);

        for i=1:Nx
            for j=1:Ny
                x=(i-cx)*dx; y=(j-cy)*dy;
                if TwoSpheres
                    if sqrt(x^2+(y-gap_cells*dy-R)^2)<=R || ...
                       sqrt(x^2+(y+gap_cells*dy+R)^2)<=R
                        eps_r(i,j)=n_values(1)^2;
                        sigma_mat(i,j)=sigma_particle;
                    end
                else
                    if sqrt(x^2+y^2)<=R
                        eps_r(i,j)=n_values(1)^2;
                        sigma_mat(i,j)=sigma_particle;
                    end
                end
            end
        end

        vis_mask = eps_r > 1.05;

        eps_map = eps0*eps_r;
        sigma_total = sigma_e_pml + sigma_mat;

        aEz = (1 - dt*sigma_total./(2*eps_map)) ./ ...
              (1 + dt*sigma_total./(2*eps_map));
        bEz = (dt./eps_map) ./ ...
              (1 + dt*sigma_total./(2*eps_map));

        %% --- Diagnostics
        W_acc = zeros(Nx,Ny);
        U_Joule = zeros(Nx,Ny);
        E_inc = zeros(Niter,1);
        E_src_total = zeros(Niter,1);

        %% --- Video
        vidEz = VideoWriter(sprintf('Ez_%s_R%.1fum.mp4',run_title,R*1e6),'MPEG-4');
        vidNano = VideoWriter(sprintf('Nanojet_%s_R%.1fum.mp4',run_title,R*1e6),'MPEG-4');
        vidEz.FrameRate = vid_fps;
        vidNano.FrameRate = vid_fps;
        open(vidEz); open(vidNano);

        figEz = figure('Color','w');
        figNano = figure('Color','w');

        Imax_smooth = 0;
        Wmax_smooth = 0;

        %% --- TIME LOOP
        for n=1:Niter
            t = n*dt;

            Hx(:,1:end-1)=aHx(:,1:end-1).*Hx(:,1:end-1) ...
                - bHx(:,1:end-1).*(Ez(:,2:end)-Ez(:,1:end-1))/dy;
            Hy(1:end-1,:)=aHy(1:end-1,:).*Hy(1:end-1,:) ...
                + bHy(1:end-1,:).*(Ez(2:end,:)-Ez(1:end-1,:))/dx;

            curlH=(Hy-[zeros(1,Ny);Hy(1:end-1,:)])/dx ...
                 -(Hx-[zeros(Nx,1),Hx(:,1:end-1)])/dy;

            Ez = aEz.*Ez + bEz.*curlH;

            src = exp(-((t-t0)/tau)^2).*sin(omega*t);
            Ez(src_x,:) = Ez(src_x,:) + src*spatial_profile;

            % Diagnostics accumulation
            W_acc = W_acc + abs(Ez).^2*dt;
            U_Joule = U_Joule + 0.5*sigma_mat.*abs(Ez).^2*dt;
            E_inc(n) = Ez(src_x + round(1e-6/dx), cy);
            E_src_total(n) = Ez(src_x, cy);

            % Instantaneous field
            if mod(n,snap_interval)==0
                Imax_now = max(abs(Ez(:)).^2);
                Imax_smooth = 0.95*Imax_smooth + 0.05*Imax_now;

                figure(figEz); clf;
                imagesc(abs(Ez').^2);
                axis image; colormap hot; colorbar;
                xlabel('x (µm)')
                ylabel('y (µm)')
                caxis([0 Imax_smooth+eps]);
                hold on; contour(vis_mask',[1 1],'c','LineWidth',1.2); hold off;
                title(sprintf('%s | |Ez|^2 | t=%.1f fs',run_title,t*1e15));
                drawnow;
                writeVideo(vidEz,getframe(figEz));
            end

            % Nanojet accumulation
            if mod(n,snap_interval)==0 && n > round(1.5*t0/dt)
                Wmax_now = max(W_acc(:));
                Wmax_smooth = 0.95*Wmax_smooth + 0.05*Wmax_now;

                figure(figNano); clf;
                imagesc(W_acc');
                axis image; colormap jet; colorbar;
                xlabel('x (µm)')
                ylabel('y (µm)')
                caxis([0 Wmax_smooth+eps]);
                hold on; contour(vis_mask',[1 1],'w','LineWidth',1.2); hold off;
                title(sprintf('%s | Nanojet Build-Up',run_title));
                drawnow;
                writeVideo(vidNano,getframe(figNano));
            end
        end

        close(vidEz); close(vidNano);

        %% =========================================================================
        % ANALYSIS (RESTORED)
        % =========================================================================
        Ez_final = Ez;

        E_inc_peak = max(abs(E_inc));
        E_local_max = max(abs(Ez_final(:)));
        FEF = E_local_max / E_inc_peak;

        [~, x_peak_idx] = max(max(abs(Ez_final).^2,[],2));
        y_profile = abs(Ez_final(x_peak_idx,:)).^2;
        y_profile = y_profile / max(y_profile);
        idx = find(y_profile >= 0.5);
        if isempty(idx)
            FWHM_um = NaN;
        else
            FWHM_um = (idx(end)-idx(1)+1)*dy*1e6;
        end

        centerline = abs(Ez_final(:,cy)).^2;
        idx = find(centerline >= 0.5*max(centerline));
        if isempty(idx)
            nanojet_length = NaN;
        else
            nanojet_length = (idx(end)-idx(1))*dx*1e6;
        end

        Total_Joule = sum(U_Joule(:))*dx*dy;

        E_ref = E_src_total - E_inc;
        R_coeff = max(abs(E_ref).^2)/max(abs(E_inc).^2);

        Ez_centerline = squeeze(Ez_final(:,cy));
        fft_spectrum = abs(fft(Ez_centerline));
        fft_axis = (0:Nx-1)/(Nx*dx);

       %% === Fourier spectrum along nanojet axis (x) ===
        Ez_centerline = squeeze(Ez_final(:,cy));
        Ez_centerline = Ez_centerline - mean(Ez_centerline); % remove DC
        
        fft_spectrum = abs(fftshift(fft(Ez_centerline)));
        k_axis = (-Nx/2:Nx/2-1)/(Nx*dx);   % spatial frequency (1/m)
        
        % === CHECK PLOT ===
        x_axis = ((1:Nx)-0.5*Nx)*dx*1e6;   % µm
        
        figure('Color','w');
        
        subplot(2,1,1)
        plot(x_axis, real(Ez_centerline),'LineWidth',1.2)
        xlabel('x (µm)')
        ylabel('Ez')
        title('Electric field Ez along nanojet axis')
        grid on
        
        subplot(2,1,2)
        plot(k_axis*1e-6, fft_spectrum/max(fft_spectrum),'LineWidth',1.2)
        xlabel('Spatial frequency k_x (µm^{-1})')
        ylabel('Normalized |FFT(Ez)|')
        title('Spatial Fourier spectrum along x')
        grid on


        %% Store
        Results(run_idx,r_idx).Radius = R;
        Results(run_idx,r_idx).FEF = FEF;
        Results(run_idx,r_idx).FWHM_um = FWHM_um;
        Results(run_idx,r_idx).NanojetLength_um = nanojet_length;
        Results(run_idx,r_idx).Total_Joule = Total_Joule;
        Results(run_idx,r_idx).Reflection = R_coeff;
        Results(run_idx,r_idx).FFT_axis = fft_axis;
        Results(run_idx,r_idx).FFT_spectrum = fft_spectrum;

        FEF_all(r_idx) = FEF;
        FWHM_all(r_idx) = FWHM_um;
        NanojetL_all(r_idx) = nanojet_length;
        Joule_all(r_idx) = Total_Joule;
        Reflection_all(r_idx) = R_coeff;

        % Final centerline plot
        x_axis = ((1:Nx)-0.5*Nx)*dx*1e6;
        figFinal = figure('Color','w');
        plot(x_axis, abs(Ez_final(:,cy)).^2,'LineWidth',1.4);
        xlabel('x (µm)'); ylabel('|Ez|^2');
        title(sprintf('%s | R=%.2f µm | FEF=%.2f',run_title,R*1e6,FEF));
        grid on;
        saveas(figFinal, sprintf('Final_%s_R%.1fum.png', ...
            strrep(run_title,' ','_'), R*1e6));
    end

    %% --- Summary plots
    figure('Color','w');
    subplot(2,2,1);
    plot(R_values*1e6, FEF_all,'-o','LineWidth',1.5);
    xlabel('R (µm)'); ylabel('FEF'); grid on;

    subplot(2,2,2);
    plot(R_values*1e6, FWHM_all,'-o','LineWidth',1.5);
    xlabel('R (µm)'); ylabel('FWHM (µm)'); grid on;

    subplot(2,2,3);
    plot(R_values*1e6, NanojetL_all,'-o','LineWidth',1.5);
    xlabel('R (µm)'); ylabel('Nanojet Length (µm)'); grid on;

    subplot(2,2,4);
    plot(R_values*1e6, Reflection_all,'-o','LineWidth',1.5);
    xlabel('R (µm)'); ylabel('Reflection'); grid on;

    sgtitle(sprintf('%s Summary', run_title));
end
%%
k0 = 2*pi/lambda;

Ex = zeros(size(Ez_final));
Ex(2:end-1,:) = -(1/(1i*k0)) * ...
    (Ez_final(3:end,:) - Ez_final(1:end-2,:)) / (2*dx);

%%
figure;
imagesc(Ez');
axis image; colorbar;
xlabel('x (µm)')
ylabel('y (µm)')
title('Electric field Ez');
%%
figure;
imagesc(Hx');
axis image; colorbar;
xlabel('x (µm)')
ylabel('y (µm)')
title('Magnetic field Hx');
%%
figure;
imagesc(Hy');
axis image; colorbar;
xlabel('x (µm)')
ylabel('y (µm)')
title('Magnetic field Hy');
%%
figure;
imagesc(abs(Ex').^2);
axis image; colorbar;
xlabel('x (µm)')
ylabel('y (µm)')
title('Longitudinal electric field |Ex|^2');
