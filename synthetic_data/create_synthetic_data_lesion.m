%% Simulate dynamic deuterium MRSI data, based on in vivo maps of CSF, WM and GM %%
% Written by Anna Duguid, Jun 2025

clearvars;
addpath(genpath('../MRSI_functions'))

%% Setting parameters for simulation

load(strcat(pwd,'/LesionPars.mat'))
temp_data_hr.RecoPar = LesionPars;
temp_data_hr.Par = LesionPars;
temp_data_lr.RecoPar = LesionPars;
temp_data_lr.Par = LesionPars;

Settings.Phase_Offsets_hr_flag = true;
Settings.B0_Inhomogeneity_flag = true; % Apply B0 map to high resolution data
Settings.lambda = 1; % Scaling applied to B0 map before application to high resolution spatiospectral phantom
Settings.Simulate_Masks_lr_flag = true; % Generate downsampled GM, WM, CSF maps, and low res mask
Settings.ApplyFreqAlign_lr_flag = true; % Perform frequency alignment for downsampled data following frequency shifts of high res data
Settings.CorrectT2_flag = true; % Apply tissue specific T2 times to metabolite signals and multiply with corresponding tissue probability maps
Settings.Hamming_flag = false; % Apply hamming filter to truncated k-space data

%% Fixing parameters and loading data

if Settings.ApplyFreqAlign_lr_flag; Settings.Simulate_Masks_lr_flag = true; end

%% Reading in tissue maps, B0 map and brain masks

load(strcat(pwd,'/maps/mask.mat'));
load(strcat(pwd,'/maps/GM.mat'));
load(strcat(pwd,'/maps/WM.mat'));
load(strcat(pwd,'/maps/CSF.mat'));

for ii=1:size(mask_sparse, 2)
    mask_hr(:,:,ii)=full(mask_sparse{ii});
    GM_hr(:,:,ii)=full(GM_sparse{ii});
    WM_hr(:,:,ii)=full(WM_sparse{ii});
    CSF_hr(:,:,ii)=full(CSF_sparse{ii});
end

clear mask_hr_sparse GM_hr_sparse WM_hr_sparse CSF_hr_sparse

mask_hr = logical(mask_hr);
GM_hr = single(GM_hr);
WM_hr = single(WM_hr);
CSF_hr = single(CSF_hr);

% Defining spherical lesion mask in high resolution
Lesion_centre = [110 130 165];
Lesion_radius = 18;
[x, y, z] = ndgrid(1:size(mask_hr, 1), 1:size(mask_hr, 2), 1:size(mask_hr, 3));
Lesion_hr = sqrt((x - Lesion_centre(1)).^2 + (y - Lesion_centre(2)).^2 + (z - Lesion_centre(3)).^2) <= Lesion_radius;

GM_hr(Lesion_hr==1) = 0;
WM_hr(Lesion_hr==1) = 0;
CSF_hr(Lesion_hr==1) = 0;

% Loading B0 map and converting from proton to deuteron frequencies

if(Settings.B0_Inhomogeneity_flag)
    
    load(strcat(pwd,'/maps/B0_ASPIRE.mat'));
    for ii=1:size(B0_ASPIRE_sparse, 2)
        B0_Map(:,:,ii)=full(B0_ASPIRE_sparse{ii});
    end
    
    clear B0_ASPIRE_sparse
    
    B0_Map = single(B0_Map);
    B0_Map = B0_Map * 6.536 / 42.576 * Settings.lambda; % Accounting for gyromagnetic ratio of deuterium and scaling by lambda parameter
    
end

if(Settings.Phase_Offsets_hr_flag)
    
    load(strcat(pwd,'/maps/phase.mat'));
    for ii=1:size(phase_sparse, 2)
        phase_map_hr(:,:,ii)=full(phase_sparse{ii});
    end
    
    clear phase_sparse
    
    phase_map_hr = single(phase_map_hr);
    phase_map_hr = phase_map_hr * 6.536 / 42.576; % Scaling based on gyromagnetic ratio
    phase_map_hr = exp(1i*phase_map_hr).* mask_hr; 
    phase_map_hr(isnan(phase_map_hr)) = 0;
    phase_map_hr = repmat(phase_map_hr, [1 1 1 LesionPars.vecSize]);

end

%% Indices for for downsampling by k-space truncation

start_index = floor((LesionPars.ImageSize_hr - LesionPars.TargetImageSize) / 2) + 1;
end_index = start_index + LesionPars.TargetImageSize - 1;

%% Generating low resolution maps and brain mask
if Settings.Simulate_Masks_lr_flag

    % Generating lr brain mask

    mask_lr = imresize3(single(mask_hr), [22,22,21], 'nearest');
    mask_lr = logical(mask_lr);
    
    % Make sum map and calculate in HR what the 'air' component is
    Sum_Map_hr = WM_hr + GM_hr + CSF_hr + Lesion_hr;
    Non_MR_hr = 1-Sum_Map_hr;

    % 3D Fourier transform to k-space
    
    GM_lr = fftshift(fft(ifftshift(GM_hr,1), [], 1), 1);
    GM_lr = fftshift(fft(ifftshift(GM_lr,2), [], 2), 2);
    GM_lr = fftshift(fft(ifftshift(GM_lr,3), [], 3), 3);

    WM_lr = fftshift(fft(ifftshift(WM_hr,1), [], 1), 1);
    WM_lr = fftshift(fft(ifftshift(WM_lr,2), [], 2), 2);
    WM_lr = fftshift(fft(ifftshift(WM_lr,3), [], 3), 3);
    
    CSF_lr = fftshift(fft(ifftshift(CSF_hr,1), [], 1), 1);
    CSF_lr = fftshift(fft(ifftshift(CSF_lr,2), [], 2), 2);
    CSF_lr = fftshift(fft(ifftshift(CSF_lr,3), [], 3), 3);
    
    Non_lr = fftshift(fft(ifftshift(Non_MR_hr,1), [], 1), 1);
    Non_lr = fftshift(fft(ifftshift(Non_lr,2), [], 2), 2);
    Non_lr = fftshift(fft(ifftshift(Non_lr,3), [], 3), 3);
    
    Lesion_lr = fftshift(fft(ifftshift(Lesion_hr,1), [], 1), 1);
    Lesion_lr = fftshift(fft(ifftshift(Lesion_lr,2), [], 2), 2);
    Lesion_lr = fftshift(fft(ifftshift(Lesion_lr,3), [], 3), 3);

    % Cutting out desired resolution

    WM_lr = WM_lr(start_index(1):end_index(1),start_index(2):end_index(2),start_index(3):end_index(3)); % Check these numbers are right
    GM_lr = GM_lr(start_index(1):end_index(1),start_index(2):end_index(2),start_index(3):end_index(3)); % Check these numbers are right
    CSF_lr = CSF_lr(start_index(1):end_index(1),start_index(2):end_index(2),start_index(3):end_index(3)); % Check these numbers are right
    Non_lr = Non_lr(start_index(1):end_index(1),start_index(2):end_index(2),start_index(3):end_index(3)); % Check these numbers are right
    Lesion_lr = Lesion_lr(start_index(1):end_index(1),start_index(2):end_index(2),start_index(3):end_index(3));
    
    GM_lr=HammingFilter(GM_lr,[1 2 3],100,'OuterProduct',1);
    WM_lr=HammingFilter(WM_lr,[1 2 3],100,'OuterProduct',1);
    CSF_lr=HammingFilter(CSF_lr,[1 2 3],100,'OuterProduct',1);
    Lesion_lr=HammingFilter(Lesion_lr,[1 2 3],100,'OuterProduct',1);
    Non_lr=HammingFilter(Non_lr,[1 2 3],100,'OuterProduct',1);

    % FoV shift correction

    for x=1:size(GM_lr,1)
        GM_lr(x,:,:)=GM_lr(x,:,:)*exp(1i*pi*x/size(GM_lr,1));
        WM_lr(x,:,:)=WM_lr(x,:,:)*exp(1i*pi*x/size(WM_lr,1));
        CSF_lr(x,:,:)=CSF_lr(x,:,:)*exp(1i*pi*x/size(CSF_lr,1));
        Lesion_lr(x,:,:)=Lesion_lr(x,:,:)*exp(1i*pi*x/size(Lesion_lr,1));
        Non_lr(x,:,:)=Non_lr(x,:,:)*exp(1i*pi*x/size(Non_lr,1));
    end

    for y=1:size(GM_lr,2)
        GM_lr(:,y,:)=GM_lr(:,y,:)*exp(1i*pi*y/size(GM_lr,2));
        WM_lr(:,y,:)=WM_lr(:,y,:)*exp(1i*pi*y/size(WM_lr,2));
        CSF_lr(:,y,:)=CSF_lr(:,y,:)*exp(1i*pi*y/size(CSF_lr,2));
        Lesion_lr(:,y,:)=Lesion_lr(:,y,:)*exp(1i*pi*y/size(Lesion_lr,2));
        Non_lr(:,y,:)=Non_lr(:,y,:)*exp(1i*pi*y/size(Non_lr,2));
    end

    % k-space to image space

    GM_lr = fftshift(ifft(ifftshift(GM_lr,1), [], 1), 1);
    GM_lr = fftshift(ifft(ifftshift(GM_lr,2), [], 2), 2);
    GM_lr = abs(fftshift(ifft(ifftshift(GM_lr,3), [], 3), 3));

    WM_lr = fftshift(ifft(ifftshift(WM_lr,1), [], 1), 1);
    WM_lr = fftshift(ifft(ifftshift(WM_lr,2), [], 2), 2);
    WM_lr = abs(fftshift(ifft(ifftshift(WM_lr,3), [], 3), 3));

    CSF_lr = fftshift(ifft(ifftshift(CSF_lr,1), [], 1), 1);
    CSF_lr = fftshift(ifft(ifftshift(CSF_lr,2), [], 2), 2);
    CSF_lr = abs(fftshift(ifft(ifftshift(CSF_lr,3), [], 3), 3));
    
    Lesion_lr = fftshift(ifft(ifftshift(Lesion_lr,1), [], 1), 1);
    Lesion_lr = fftshift(ifft(ifftshift(Lesion_lr,2), [], 2), 2);
    Lesion_lr = abs(fftshift(ifft(ifftshift(Lesion_lr,3), [], 3), 3));
    
    Non_lr = fftshift(ifft(ifftshift(Non_lr,1), [], 1), 1);
    Non_lr = fftshift(ifft(ifftshift(Non_lr,2), [], 2), 2);
    Non_lr = abs(fftshift(ifft(ifftshift(Non_lr,3), [], 3), 3));
    
    GM_lr(mask_lr==0) = 0;
    WM_lr(mask_lr==0) = 0;
    CSF_lr(mask_lr==0) = 0;
    Lesion_lr(mask_lr==0) = 0;
    Non_lr(Non_lr==0) = 0;
    summed_Map_lr = GM_lr + WM_lr + CSF_lr + Lesion_lr + Non_lr;

    % Normalising low rank masks

    GM_lr = (GM_lr ./ summed_Map_lr); GM_lr(mask_lr==0) = 0;
    WM_lr = (WM_lr ./ summed_Map_lr); WM_lr(mask_lr==0) = 0;
    CSF_lr = (CSF_lr ./ summed_Map_lr); CSF_lr(mask_lr==0) = 0;
    Lesion_lr = (Lesion_lr ./ summed_Map_lr); Lesion_lr(mask_lr==0) = 0;

    % Saving low resolution maps and brain mask as nifti files
    niftiwrite(single(mask_lr), strcat(pwd, '/maps/mask_lr.nii'));
    niftiwrite(single(GM_lr), strcat(pwd, '/maps/GM_lr.nii'));
    niftiwrite(single(WM_lr), strcat(pwd, '/maps/WM_lr.nii'));
    niftiwrite(single(CSF_lr), strcat(pwd, '/maps/CSF_lr.nii'));
    niftiwrite(single(Lesion_lr), strcat(pwd, '/maps/Lesion_lr.nii'));

end

%% Simulating characteristic FID and temporal changes for each tissue type (GM, WM, and CSF)

% Loading metabolite signals

load(strcat(pwd, '/metabolites/glc_lesion.mat'));
load(strcat(pwd, '/metabolites/glx_lesion.mat'));
load(strcat(pwd, '/metabolites/lac_lesion.mat'));
load(strcat(pwd, '/metabolites/water_lesion.mat'));

if Settings.CorrectT2_flag
    
    % Removing basis function T2 decay
    Ts = 1:1:size(water,2);
    fitt=polyfit(Ts, log(water),1);  
    CorrectionFactors = exp(Ts*(-fitt(1)));
    glc = glc.*CorrectionFactors;
    glx = glx.*CorrectionFactors;
    lac = lac.*CorrectionFactors;
    water = water.*CorrectionFactors;
    
    % Bader 2024: T2 times (ms) for GM/WM
    % Cocking et al: T2 CSF/GM/WM = 90/32/30 - assuming a T2 ratio of 90/31 to calculate expected CSF T2 times...
        
    t2_glc_GM = 36;
    t2_glc_WM = 35;
    t2_glc_lesion = (t2_glc_GM+t2_glc_WM)/2;
    t2_glc_CSF = 103.1;

    t2_glx_GM = 37;
    t2_glx_WM = 34;
    t2_glx_lesion = (t2_glx_GM+t2_glx_WM)/2;
    t2_glx_CSF = 103.1;

    t2_lac_lesion = t2_glc_lesion*61/32; % Using ratio at 11.7T between glucose and lactate t2 times

    t2_water_GM = 36;
    t2_water_WM = 33;
    t2_water_lesion = (t2_water_GM+t2_water_WM)/2;
    t2_water_CSF = 100.2;

    % FID acquisition times in ms, NOTE: basis set already includes a 2ms acquisition delay so these are actually the FID_times - 2ms
    add_DelayTime_ms = 0;
    FID_times = add_DelayTime_ms + (0:1:LesionPars.vecSize-1) * LesionPars.Dwelltimes * 1e-6;

    % Applying t2 decay - separation for GM / WM / CSF
    glc_GM = reshape(glc.*exp(-FID_times/t2_glc_GM), [1 1 1 LesionPars.vecSize]);
    glc_WM = reshape(glc.*exp(-FID_times/t2_glc_WM), [1 1 1 LesionPars.vecSize]);
    glc_lesion = reshape(glc.*exp(-FID_times/t2_glc_lesion), [1 1 1 LesionPars.vecSize]);
    glc_CSF = reshape(glc.*exp(-FID_times/t2_glc_CSF), [1 1 1 LesionPars.vecSize]);
    
    glx_GM = reshape(glx.*exp(-FID_times/t2_glx_GM), [1 1 1 LesionPars.vecSize]);
    glx_WM = reshape(glx.*exp(-FID_times/t2_glx_WM), [1 1 1 LesionPars.vecSize]);
    glx_lesion = reshape(glx.*exp(-FID_times/t2_glx_lesion), [1 1 1 LesionPars.vecSize]);
    glx_CSF = reshape(glx.*exp(-FID_times/t2_glx_CSF), [1 1 1 LesionPars.vecSize]);

    lac_lesion = reshape(lac.*exp(-FID_times/t2_lac_lesion), [1 1 1 LesionPars.vecSize]);
    
    water_GM = reshape(water.*exp(-FID_times/t2_water_GM), [1 1 1 LesionPars.vecSize]);
    water_WM = reshape(water.*exp(-FID_times/t2_water_WM), [1 1 1 LesionPars.vecSize]);
    water_lesion = reshape(water.*exp(-FID_times/t2_water_lesion), [1 1 1 LesionPars.vecSize]);
    water_CSF = reshape(water.*exp(-FID_times/t2_water_CSF), [1 1 1 LesionPars.vecSize]);

else
    
    add_DelayTime_ms = 0;
    
    glc_GM = reshape(glc, [1 1 1 LesionPars.vecSize]);
    glc_WM = reshape(glc, [1 1 1 LesionPars.vecSize]);
    glc_lesion = reshape(glc, [1 1 1 LesionPars.vecSize]);
    glc_CSF = reshape(glc, [1 1 1 LesionPars.vecSize]);
    
    glx_GM = reshape(glx, [1 1 1 LesionPars.vecSize]);
    glx_WM = reshape(glx, [1 1 1 LesionPars.vecSize]);
    glx_lesion = reshape(glx, [1 1 1 LesionPars.vecSize]);
    glx_CSF = reshape(glx, [1 1 1 LesionPars.vecSize]);
    
    lac_lesion = reshape(lac, [1 1 1 LesionPars.vecSize]);
    
    water_GM = reshape(water, [1 1 1 LesionPars.vecSize]);
    water_WM = reshape(water, [1 1 1 LesionPars.vecSize]);
    water_lesion = reshape(water, [1 1 1 LesionPars.vecSize]);
    water_CSF = reshape(water, [1 1 1 LesionPars.vecSize]);
    
end


% Temporal behaviour of metabolites over 14-63 min (based on in vivo observation)

glx_dyn_WM = single(linspace(0.25,1,LesionPars.nRep)); 
glx_dyn_CSF = single(linspace(0.25,1,LesionPars.nRep));
glx_dyn_lesion = single(linspace(0.25,1,LesionPars.nRep));
lac_dyn_lesion = single(linspace(0.25,1,LesionPars.nRep)); 
glx_dyn_GM = single(linspace(0.1,1.3,LesionPars.nRep));
glc_dyn = single(1-exp(-(0:LesionPars.nRep-1:55)/20)+0.0381);
water_dyn = single(linspace(0.9,1,LesionPars.nRep));

%% Simulating csi data

% Generating 4D data seperately for each repetition (memory limit) and then cutting out required k-space and saving low resolution data

csi_data_lr.Data = zeros([LesionPars.TargetImageSize LesionPars.vecSize LesionPars.nRep]);
temp_data_hr.Data = zeros([LesionPars.ImageSize_hr LesionPars.vecSize]);

for Rep = 1:LesionPars.nRep
    
    fprintf('\n Generating 2H-MRSI data repetition %d of 8 \n', Rep);
        
    temp_data_hr.Data = 0.7 * (glc_lesion .* Lesion_hr + glc_GM .* GM_hr + 0.7 * glc_WM .* WM_hr + glc_CSF .* CSF_hr) * glc_dyn(Rep); 
        
    temp_data_hr.Data = temp_data_hr.Data + 0.8 * (0.1 * glx_lesion .* Lesion_hr + 1 * glx_GM .* GM_hr * glx_dyn_GM(Rep) + 0.25 * glx_WM .* WM_hr * glx_dyn_WM(Rep) +  0.125 * glx_CSF .* CSF_hr * glx_dyn_CSF(Rep)); 
    
    temp_data_hr.Data = temp_data_hr.Data + 2 * (1.25 * water_lesion .* Lesion_hr + 1.25 * water_GM .* GM_hr + water_WM .* WM_hr + 1.5 * water_CSF .* CSF_hr) * water_dyn(Rep);
    
    temp_data_hr.Data = temp_data_hr.Data + lac_lesion .* Lesion_hr * lac_dyn_lesion(Rep);
    
    % Add BO inhomogeneity here
    
    if Settings.B0_Inhomogeneity_flag
        AdditionalInput.B0Map_Hz = B0_Map;
        AdditionalInput.Mask = mask_hr;
        Settings.Printing = 1;
        temp_data_hr.Mask = mask_hr;
        Settings.ApplyAlongDim = 4;
        Settings.MaxShifts_ppm = [-0.3 0.3];
        [temp_data_hr,~] = op_FrequencyAlign2_AD(temp_data_hr,Settings,AdditionalInput);
        clear AdditionalInput;
    end
    
    if Settings.Phase_Offsets_hr_flag
        temp_data_hr.Data = temp_data_hr.Data .* phase_map_hr;
    end
    
    % image to k-space
    temp_data_hr.Data = fftshift(fft(ifftshift(temp_data_hr.Data,1), [], 1), 1);
    temp_data_hr.Data = fftshift(fft(ifftshift(temp_data_hr.Data,2), [], 2), 2);
    temp_data_hr.Data = fftshift(fft(ifftshift(temp_data_hr.Data,3), [], 3), 3);

    % Downsampling: cut desired k-space out
    
    temp_data_lr.Data = temp_data_hr.Data(start_index(1):end_index(1),start_index(2):end_index(2),start_index(3):end_index(3),:); % Check these numbers are right

    for x=1:size(temp_data_lr.Data,1)
        temp_data_lr.Data(x,:,:,:)=temp_data_lr.Data(x,:,:,:)*exp(1i*pi*x/size(temp_data_lr.Data,1));   
    end

    for y=1:size(temp_data_lr.Data,2)
        temp_data_lr.Data(:,y,:,:)=temp_data_lr.Data(:,y,:,:)*exp(1i*pi*y/size(temp_data_lr.Data,2));   
    end

    if Settings.Hamming_flag
        temp_data_lr.Data=HammingFilter(temp_data_lr.Data,[1 2 3],100,'OuterProduct',1);
    end

    % k-space to image space
    temp_data_lr.Data = fftshift(ifft(ifftshift(temp_data_lr.Data,1), [], 1), 1);
    temp_data_lr.Data = fftshift(ifft(ifftshift(temp_data_lr.Data,2), [], 2), 2);
    temp_data_lr.Data = fftshift(ifft(ifftshift(temp_data_lr.Data,3), [], 3), 3);

    % Frequency offset correction for low resolution data
    if Settings.ApplyFreqAlign_lr_flag
        Settings.ApplyAlongDim = 4;
        Settings.ZerofillingFactor = 10;
        Settings.UseSVDForRef_flag = 1; % Using SVD for ref spec - not sure if this is right??
        Settings.UseSVDForRef_NoSingVals = 1;
        Settings.PeakSearchPPM = 4.65; % Is this right??
        Settings.Printing = 1;
        Settings.MaxShifts_ppm = [-0.5 0.5];
        AdditionalInput.Mask = mask_lr;
        [temp_data_lr,AddOut] = op_FrequencyAlign2_AD(temp_data_lr, Settings, AdditionalInput);
    end

    csi_data_lr.Data(:,:,:,:,Rep) = temp_data_lr.Data;

end

%% Saving Synthetic Data

csi_data_lr.Pars = LesionPars;
csi_data_lr.Settings = Settings;
save('deut_lesion_sim.mat', 'csi_data_lr')
clearvars -except csi_data_lr
