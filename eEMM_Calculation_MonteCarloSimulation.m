clear all
close all
clc

load('simulation_Q_matrices.mat')
load('simulation_B1data.mat')
load('experimental_B1data.mat')


%% User Input

B1_perCh_Exp = experimental_B1data;      %per-channel experimental b1 data; 4D complex data with the format of [x,y,z,channels]
B1_perCh_Sim = simulation_B1data;        %per-channel simulation b1 data; 4D complex data with the format of [x,y,z,channels]
Q = simulation_Q_matrices;               %10g-averaged local SAR data;  3D complex data with the size of [# of channels x # of channels x # of voxels]



voxSiz_fixed = [2 2 20];                 % Set the spatial resolution (mm) of the experimental data; Example here: 2mm in-plane res and 4mm slice res with 20mm slice distance
%Voxel size (mm) in directions of 1st, 2nd, and 3rd indices of the 3D experimental B1 data
voxSiz_moving = [2 2 2];                 %Set the spatial resolution (mm) of the simulation data; Example here: 2mm isotropic res  
%Voxel size (mm) in directions of 1st, 2nd, and 3rd indices of the 3D Simulation B1 data

Nc = 16;                                 %The Number of Tx Channels; Example here: 16 channel

V_exc(1:Nc,1) = exp(-1j*pi/180*(0:360/Nc:360*(1-1/Nc)).');    %Set desired sets of excitation vectors to generate the Mode of interest (Step #1); Example in this code: CP and zero-phase modes. Each excitation vector should be in the form of NofChannelsx1
V_exc(1:Nc,2) = exp(-1j*pi/180*(zeros(16,1)));           

N_monte_carlo = 1e5;                     %Number of iterations in Monte-Carlo simulations (see the manuscript for more information)

Perct = 99.9;                            %Set the percentile of error region for eEMM calculation; Example in this code: 99.9%

%% Per-channel NRMSE Calculation (Step #0: Preparation)
Nc = size(B1_perCh_Exp,4);   %The Number of Tx Channels
perChannel_NRMSE = zeros(Nc,1);
clear B1_perCh_Exp_reg B1_perCh_Sim_reg
for jj = 1:Nc  %Registering simulation per-channlel B1 maps to experimental data
    disp(['Registering Ch #',num2str(jj)])

    Vfixed  = abs(squeeze(B1_perCh_Exp(:,:,:,jj)));
    Vmoving = abs(squeeze(B1_perCh_Sim(:,:,:,jj)));
    Vmoving_phase = angle(squeeze(B1_perCh_Sim(:,:,:,jj)));
    [Vfixed_crop, Vmoving_reg, VMoving_phase_reg, Vmoving_reg_crop] = RegisterMaps(Vfixed,voxSiz_fixed,Vmoving,voxSiz_moving,Vmoving_phase);

    B1_perCh_Exp_reg(:,:,:,jj) = Vfixed;
    B1_perCh_Sim_reg(:,:,:,jj) = Vmoving_reg.*exp(1j*VMoving_phase_reg);

    Mask_reg = logical(abs(Vmoving_reg_crop));
    error_perChannel = abs(Vmoving_reg_crop - Vfixed_crop);
    perChannel_NRMSE(jj,1) = sqrt(mean(error_perChannel(Mask_reg).^2)) / mean(Vfixed_crop(Mask_reg));

    figure
    subplot(1,2,1), imagesc(MosaicGen(Vfixed_crop.*Mask_reg))
    axis equal
    axis off
    colorbar
    clim([0 4*mean(Vfixed_crop(Mask_reg))])
    colormap jet
    title('Experiment')
    subplot(1,2,2), imagesc(MosaicGen(Vmoving_reg_crop.*Mask_reg))
    axis equal
    axis off
    colorbar
    clim([0 4*mean(Vfixed_crop(Mask_reg))])
    colormap jet
    title('Simulation')
    sgtitle(['NRMSE = ',num2str(round(perChannel_NRMSE(jj)*100)),'% (Ch#',num2str(jj),')'])
    print('-djpeg','-r600',['perChannel_Registration_data Ch',num2str(jj)]);
end

Reg_chk = input('Do these regsitrations and per-channel comparisons look reasonable? Please enter y or n: ', 's');
if Reg_chk == "n"
    clc
    disp('Please check you validation steps and verify the units and scaling factors, then rerun this code')
    return
end

%% Generate the Mode of Interest (Step #1)
N_modes = size(V_exc,2);
disp(['Number of Modes of interest to be tested: ',num2str(N_modes)])
NRMSE_MOI0 = 0;
for nn = 1:N_modes
    disp(['Generating Mode #',num2str(nn)])
    pattern_exp = 0;
    pattern_sim = 0;
    for jj = 1:Nc
        pattern_exp = pattern_exp + V_exc(jj,nn) * B1_perCh_Exp(:,:,:,jj);
        pattern_sim = pattern_sim +  V_exc(jj,nn) * B1_perCh_Sim(:,:,:,jj);
    end
    Vfixed  = abs(squeeze(pattern_exp));
    Vmoving = abs(squeeze(pattern_sim));
    [Vfixed_crop, ~ , Vmoving_reg_crop] = RegisterMaps(Vfixed,voxSiz_fixed,Vmoving,voxSiz_moving);

    Mask_reg = logical(abs(Vmoving_reg_crop));
    error_MOI0 = abs(Vmoving_reg_crop - Vfixed_crop);
    NRMSE_MOI0(nn,1) = sqrt(mean(error_MOI0(Mask_reg).^2)) / mean(Vfixed_crop(Mask_reg));

    figure
    subplot(1,2,1), imagesc(MosaicGen(Vfixed_crop.*Mask_reg))
    axis equal
    axis off
    colorbar
    clim([0 2*mean(Vfixed_crop(Mask_reg))])
    colormap jet
    title('Experiment')
    subplot(1,2,2), imagesc(MosaicGen(Vmoving_reg_crop.*Mask_reg))
    axis equal
    axis off
    colorbar
    clim([0 2*mean(Vfixed_crop(Mask_reg))])
    colormap jet
    title('Simulation')
    sgtitle(['NRMSE = ',num2str(round(NRMSE_MOI0(nn)*100)),'% (MOI#',num2str(nn),')'])
    print('-djpeg','-r600',['ModeOfInterest_Registration_data_mode#',num2str(nn)]);
end

Reg_chk = input('Do these regsitrations and B1 patterns comparisons look reasonable? Please enter y or n: ', 's');
if Reg_chk == "n"
    clc
    disp('Please check you validation steps and verify the units and scaling factors, then rerun this code')
    return
end

%% Monte-Carlo Simulation and Error Spaces (Steps #2)
Perturbations = perChannel_NRMSE;
ind_Q = find(squeeze(abs(Q(1,1,:))));
Q_reduced = Q(:,:,ind_Q);
B1_Sim = reshape(B1_perCh_Sim_reg, [prod(size(B1_perCh_Sim_reg,(1:3))) , Nc]);
ind = find(abs(B1_Sim(:,1)));
NRMSE_sim = zeros(N_monte_carlo,Shim_No);
SAR_error = zeros(N_monte_carlo,Shim_No);
for Shim_No = 1:N_modes    % Monte-Carlo simulations
    X0 = V_exc(:,Shim_No);
    X0 = X0 / sqrt(X0' * X0);
    temp = B1_Sim * X0;
    B0 = temp;
    PSAR0 = pSARcalc(X0,Q_reduced);

    disp(['Shimming_NO: ',num2str(Shim_No)])
    tic
    parfor nn = 1:N_monte_carlo
        Xper = X0 + (-1).^round(rand(Nc,1)).*Perturbations.*rand(Nc,1).*real(X0) + 1j*(-1).^round(rand(Nc,1)).*Perturbations.*rand(Nc,1).*imag(X0);
        Xper = Xper / sqrt(Xper' * Xper);
        temp = B1_Sim * Xper;
        Bper = temp;
        error = abs(B0(ind)) - abs(Bper(ind));
        RMSE = sqrt( mean( (error).^2 ) );
        NRMSE_sim(nn,Shim_No) = RMSE / mean(abs(B0(ind))) * 100;

        PSARper = pSARcalc(Xper,Q_reduced);
        SAR_error(nn,Shim_No) = (PSARper - PSAR0)/PSAR0 * 100;
    end
    toc
end

%% Propagate B1 Error Region to pSAR Error Region (Steps #3)
SAR_error_region = NaN(N_monte_carlo,N_modes);
clear temp1 SAR_error_region
for ss = 1:N_modes
    temp1(:,ss) = NRMSE_sim(:,ss) < NRMSE_MOI0(ss)*100;
    temp2 = SAR_error(temp1(:,ss),ss);
    SAR_error_region(1:length(temp2),ss) = temp2;
end

%% Percentile of pSAR Error Region (Step #4)
close all
clear SAR_Modeling_error
for ss = 1:N_modes
    SAR_Modeling_error(ss) = prctile(SAR_error_region(:,ss),Perct);
    figure
    h = histogram(SAR_error_region(:,ss));
    hold on
    plot ([SAR_Modeling_error(ss) SAR_Modeling_error(ss)],[0 1000],'--r','linewidth',2)
    text(SAR_Modeling_error(ss)+5,100,[num2str(round(SAR_Modeling_error(ss))),' %'],'fontweight','bold','FontSize',18)
    set(gca,'fontsize',16)
    title(['pSAR_1_0_g EM Modeling Error (Mode #',num2str(ss),')'])
    xlabel('pSAR_1_0_g Deviation (%)')
    ylabel('# of drive vectors')
    xLabP = (round(SAR_Modeling_error(ss)/50)+1)*50;
    xLabN = (round(prctile(SAR_error_region(:,ss),0.1)/50)-1)*50;
    yLab = ceil(max(h.BinCounts) / (10^round(log10(max(h.BinCounts))))) * (10^round(log10(max(h.BinCounts))));
    xlim([xLabN xLabP])
    set(gca,'XTick',xLabN:50:xLabP);
    set(gca,'XTickLabel', xLabN:50:xLabP);
    ylim([0 yLab])
    grid on
    print('-djpeg','-r600',['pSAR_1_0_g EM Modeling Error (Mode#',num2str(ss),')']);
end

%% eEMM Calculation (Step #5)
[eEMM,pp] = max(SAR_Modeling_error);
fprintf('pSAR10g EM modeling error for tested %1d modes of interest is %1d%%, corresponding to the Mode #%1d\n\n',N_modes,round(eEMM),pp)





%% Utility Functions

function [pSAR,Pin,SARdist] = pSARcalc(Vexc,Q)
% Peak Local SAR corresponding to an excitation vector given Q-matrix/VOPs
% [pSAR,Pin] = pSAR_calc(Vexc,Qmatrix)
%
% Inputs:
%
% Vexc: Excitation Vector (peak voltage); Format: Nchx1, Complex (V)
% Q: NgAveSARmatrices; Format: NchxNchxNvox, Complex (1/kg/ohm)
%
% Outputs:
%
% pSAR: Peak Local SAR corresponding to Vexc; Format: 1x1, Positive Scalar (W/kg)
% Pin: Total input power corresponding to Vexc; Format: 1x1, Positive Scalar (W)
% SARdist: Masked SAR distribution; Format: 1xNvox, positive Scalar (W/kg)
%
% Note:
%
% To make a 3D SARdist -> SARdist3D = double(Mask); SARdist3D(Mask)= SARdist;
% where Mask is a 3D logical mask of the sample

siz = size(Vexc);
if length(Vexc)>1 && siz(1) == 1
    Vexc = permute(Vexc,[2,1]);
end

if nargout==1
    pSAR = norm(Vexc.' * reshape(Vexc' * reshape(Q,[length(Vexc),length(Vexc)*length(Q)]) , [length(Vexc),length(Q)]), 'inf');
elseif nargout==2
    pSAR = norm(Vexc.' * reshape(Vexc' * reshape(Q,[length(Vexc),length(Vexc)*length(Q)]) , [length(Vexc),length(Q)]), 'inf');
    Pin = Vexc'*Vexc/2/50;
elseif nargout==3
    pSAR = norm(Vexc.' * reshape(Vexc' * reshape(Q,[length(Vexc),length(Vexc)*length(Q)]) , [length(Vexc),length(Q)]), 'inf');
    Pin = Vexc'*Vexc/2/50;
    SARdist = real(Vexc.' * reshape(Vexc' * reshape(Q,[length(Vexc),length(Vexc)*length(Q)]) , [length(Vexc),length(Q)]));
end

end

function [Vfixed_crop, Vmoving_reg, VMoving_phase_reg, Vmoving_reg_crop] = RegisterMaps(Vfixed,voxSiz_fixed,Vmoving,voxSiz_moving,Vmoving_phase)
%Inputs: 
%
%Vfixed: Fixed refernce volume (3D Matrix)
%voxSiz_fixed: A vector containing 3 scalar numbers (Voxel size (mm) in directions of 1st, 2nd , and 3rd indices of matrix Vfixed)
%Vmoving: To-be-registred moving volume (3D Matrix)
%voxSiz_moving: A vector containing 3 scalar numbers (Voxel size (mm) in directions of 1st, 2nd , and 3rd indices of matrix Vfixed)
%Vmoving_phase: Phase of the to-be-registred moving volume (3D Matrix)
%
%Output: 
%
%Vfixed_crop: Fixed volume cropped to match the size of the moving volume
%Vmoving_reg: Moving volume registred to the fixed refeeence volume 
%VMoving_phase_reg: Phase of the moving volume registred to the fixed refeeence volume
%Vmoving_reg_crop: Moving volume registred to the fixed refeeence volume and cropped to match the size 

%Registration with different resolutions
dx_fixed = voxSiz_fixed(2);  dy_fixed = voxSiz_fixed(1);  dz_fixed = voxSiz_fixed(3);
dx_mov   = voxSiz_moving(2);  dy_mov   = voxSiz_moving(1);  dz_mov   = voxSiz_moving(3);
[sYf, sXf, sZf] = size(Vfixed);
[sYm, sXm, sZm] = size(Vmoving);
Rfixed  = imref3d(size(Vfixed), ...
    [0, dx_fixed * sXf], ...
    [0, dy_fixed * sYf], ...
    [0, dz_fixed * sZf]);
Rmoving = imref3d(size(Vmoving), ...
    [0, dx_mov * sXm], ...
    [0, dy_mov * sYm], ...
    [0, dz_mov * sZm]);
[optimizer, metric] = imregconfig('multimodal'); % or 'monomodal'
optimizer.MaximumIterations = 300;
if min(size(Vfixed))<16 || min(size(Vmoving))<16
tform = imregtform( ...
    Vmoving, Rmoving, ...
    Vfixed,  Rfixed, ...
    'rigid', optimizer, metric, ...
    'PyramidLevels', 2);
else
    tform = imregtform( ...
    Vmoving, Rmoving, ...
    Vfixed,  Rfixed, ...
    'rigid', optimizer, metric);
end

%Reslice moving volume onto fixed grid
fillVal = 0;  % or some value you know is "background"
Vmoving_reg = imwarp( ...
    Vmoving, Rmoving, ...
    tform, ...
    'OutputView', Rfixed, ...
    'FillValues', fillVal);
VMoving_phase_reg = imwarp( ...
    Vmoving_phase, Rmoving, ...
    tform, ...
    'OutputView', Rfixed, ...
    'InterpolationMethod', 'nearest', ...  % important!
    'FillValues', 0);

% Compute bounding box of the overlapping FOV and crop both
mask = Vmoving_reg ~= fillVal;
if ~any(mask(:))
    error('No overlap between fixed and moving volumes after registration!');
end
[rows, cols, slices] = ind2sub(size(mask), find(mask));
rmin = min(rows); rmax = max(rows);
cmin = min(cols); cmax = max(cols);
zmin = min(slices); zmax = max(slices);

% Crop both volumes to the same limited FOV
Vfixed_crop      = Vfixed(     rmin:rmax, cmin:cmax, zmin:zmax);
Vmoving_reg_crop = Vmoving_reg(rmin:rmax, cmin:cmax, zmin:zmax);
end

