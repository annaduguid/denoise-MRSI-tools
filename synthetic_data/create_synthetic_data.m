%% Simulate dynamic deuterium MRSI data, based on in vivo maps of CSF, WM and GM %%
% Written by Anna Duguid, Jun 2025

clearvars;
addpath(genpath('../MRSI_functions'))

%% Setting parameters for simulation

load(strcat(pwd,'/Pars.mat'))
temp_data_hr.RecoPar = Pars;
temp_data_hr.Par = Pars;
temp_data_lr.RecoPar = Pars;
temp_data_lr.Par = Pars;

Settings.Phase_Offsets_hr_flag = true;
Settings.B0_Inhomogeneity_flag = true; % Apply B0 map to high resolution data
Settings.lambda = 1; % Scaling applied to B0 map before application to high resolution spatiospectral phantom
Settings.Simulate_Masks_lr_flag = true; % Generate downsampled GM, WM, CSF maps, and low res mask
Settings.ApplyFreqAlign_lr_flag = true; % Perform frequency alignment for downsampled data following frequency shifts of high res data
Settings.CorrectT2_flag = true; % Apply tissue specific T2 times to metabolite signals and multiply with corresponding tissue probability maps
Settings.Hamming_flag = false; % Apply hamming filter to truncated k-space data

%% Fixing parameters and loading data

if(Settings.ApplyFreqAlign_lr_flag); Settings.Simulate_Masks_lr_flag = true; end

%% Reading in sparse tissue maps, B0 map and brain masks

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
    phase_map_hr = repmat(phase_map_hr, [1 1 1 Pars.vecSize]);

end

%% Indices for for downsampling by k-space truncation

start_index = floor((Pars.ImageSize_hr - Pars.TargetImageSize) / 2) + 1;
end_index = start_index + Pars.TargetImageSize - 1;

%% Generating low resolution maps and brain mask

if(Settings.Simulate_Masks_lr_flag)

    % Generating lr brain mask

    mask_lr = imresize3(single(mask_hr), [22,22,21], 'nearest');
    mask_lr = logical(mask_lr);
    
    % Make sum map and calculate in HR what the 'air' component is
    Sum_Map_hr = WM_hr + GM_hr + CSF_hr;
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
    
    Non_MR_lr = fftshift(fft(ifftshift(Non_MR_hr,1), [], 1), 1);
    Non_MR_lr = fftshift(fft(ifftshift(Non_MR_lr,2), [], 2), 2);
    Non_MR_lr = fftshift(fft(ifftshift(Non_MR_lr,3), [], 3), 3);

    % Cutting out desired resolution

    WM_lr = WM_lr(start_index(1):end_index(1),start_index(2):end_index(2),start_index(3):end_index(3)); % Check these numbers are right
    GM_lr = GM_lr(start_index(1):end_index(1),start_index(2):end_index(2),start_index(3):end_index(3)); % Check these numbers are right
    CSF_lr = CSF_lr(start_index(1):end_index(1),start_index(2):end_index(2),start_index(3):end_index(3)); % Check these numbers are right
    Non_MR_lr = Non_MR_lr(start_index(1):end_index(1),start_index(2):end_index(2),start_index(3):end_index(3)); % Check these numbers are right
    
    GM_lr=HammingFilter(GM_lr,[1 2 3],100,'OuterProduct',1);
    WM_lr=HammingFilter(WM_lr,[1 2 3],100,'OuterProduct',1);
    CSF_lr=HammingFilter(CSF_lr,[1 2 3],100,'OuterProduct',1);
    Non_MR_lr=HammingFilter(Non_MR_lr,[1 2 3],100,'OuterProduct',1);

    % FoV shift correction
    for x=1:size(GM_lr,1)
        GM_lr(x,:,:)=GM_lr(x,:,:)*exp(1i*pi*x/size(GM_lr,1));
        WM_lr(x,:,:)=WM_lr(x,:,:)*exp(1i*pi*x/size(WM_lr,1));
        CSF_lr(x,:,:)=CSF_lr(x,:,:)*exp(1i*pi*x/size(CSF_lr,1));
        Non_MR_lr(x,:,:)=Non_MR_lr(x,:,:)*exp(1i*pi*x/size(Non_MR_lr,1));
    end

    for y=1:size(GM_lr,2)
        GM_lr(:,y,:)=GM_lr(:,y,:)*exp(1i*pi*y/size(GM_lr,2));
        WM_lr(:,y,:)=WM_lr(:,y,:)*exp(1i*pi*y/size(WM_lr,2));
        CSF_lr(:,y,:)=CSF_lr(:,y,:)*exp(1i*pi*y/size(CSF_lr,2));
        Non_MR_lr(x,:,:)=Non_MR_lr(x,:,:)*exp(1i*pi*x/size(Non_MR_lr,1));
    end

    % Fourier transform back to image space

    GM_lr = fftshift(ifft(ifftshift(GM_lr,1), [], 1), 1);
    GM_lr = fftshift(ifft(ifftshift(GM_lr,2), [], 2), 2);
    GM_lr = abs(fftshift(ifft(ifftshift(GM_lr,3), [], 3), 3));

    WM_lr = fftshift(ifft(ifftshift(WM_lr,1), [], 1), 1);
    WM_lr = fftshift(ifft(ifftshift(WM_lr,2), [], 2), 2);
    WM_lr = abs(fftshift(ifft(ifftshift(WM_lr,3), [], 3), 3));

    CSF_lr = fftshift(ifft(ifftshift(CSF_lr,1), [], 1), 1);
    CSF_lr = fftshift(ifft(ifftshift(CSF_lr,2), [], 2), 2);
    CSF_lr = abs(fftshift(ifft(ifftshift(CSF_lr,3), [], 3), 3));
    
    Non_MR_lr = fftshift(ifft(ifftshift(Non_MR_lr,1), [], 1), 1);
    Non_MR_lr = fftshift(ifft(ifftshift(Non_MR_lr,2), [], 2), 2);
    Non_MR_lr = abs(fftshift(ifft(ifftshift(Non_MR_lr,3), [], 3), 3));
    
    GM_lr(mask_lr==0) = 0;
    WM_lr(mask_lr==0) = 0;
    CSF_lr(mask_lr==0) = 0;
    Non_MR_lr(mask_lr==0) = 0;

    summed_Map_lr = GM_lr + WM_lr + CSF_lr + Non_MR_lr;

    % Normalising low rank masks

    GM_lr = (GM_lr ./ summed_Map_lr); GM_lr(mask_lr==0) = 0;
    WM_lr = (WM_lr ./ summed_Map_lr); WM_lr(mask_lr==0) = 0;
    CSF_lr = (CSF_lr ./ summed_Map_lr); CSF_lr(mask_lr==0) = 0;

    % Saving low resolution maps and brain mask as nifti files
    niftiwrite(single(mask_lr), strcat(pwd, '/maps/mask_lr.nii'));
    niftiwrite(single(GM_lr), strcat(pwd, '/maps/GM_lr.nii'));
    niftiwrite(single(WM_lr), strcat(pwd, '/maps/WM_lr.nii'));
    niftiwrite(single(CSF_lr), strcat(pwd, '/maps/CSF_lr.nii'));

end

%% Simulating characteristic FID and temporal changes for each tissue type (GM, WM, and CSF)

% Loading metabolite signals

glc = load(strcat(pwd, '/metabolites/glc.mat'));
glc = single(conj(glc.data));
glx = load(strcat(pwd, '/metabolites/glx.mat'));
glx = single(conj(glx.data));
water = load(strcat(pwd, '/metabolites/water.mat'));
water = single(conj(water.data));

if(Settings.CorrectT2_flag)
    % Removing basis function T2 decay
    Ts = 1:1:96;
    fitt=polyfit(Ts, log(water),1);  
    CorrectionFactors = exp(Ts*(-fitt(1)));
    glc = glc.*CorrectionFactors;
    glx = glx.*CorrectionFactors;
    water = water.*CorrectionFactors;

    % Bader 2024: T2 times (ms) for GM/WM
    % Cocking et al: T2 CSF/GM/WM = 90/32/30 - assuming a T2 ratio of 90/31 to calculate expected CSF T2 times...
       
    t2_glc_GM = 36;
    t2_glc_WM = 35;
    t2_glc_CSF = 103.1;

    t2_glx_GM = 37;
    t2_glx_WM = 34;
    t2_glx_CSF = 103.1;

    t2_water_GM = 36;
    t2_water_WM = 33;
    t2_water_CSF = 100.2;
        

    % FID acquisition times in ms, NOTE: basis set already includes a 2ms acquisition delay so these are actually the FID_times - 2ms
    add_DelayTime_ms = 0;
    FID_times = add_DelayTime_ms + (0:1:Pars.vecSize-1) * Pars.Dwelltimes * 1e-6;

    % Applying t2 decay - separation for GM / WM / CSF
    glc_GM = reshape(glc.*exp(-FID_times/t2_glc_GM), [1 1 1 Pars.vecSize]);
    glc_WM = reshape(glc.*exp(-FID_times/t2_glc_WM), [1 1 1 Pars.vecSize]);
    glc_CSF = reshape(glc.*exp(-FID_times/t2_glc_CSF), [1 1 1 Pars.vecSize]);
    
    glx_GM = reshape(glx.*exp(-FID_times/t2_glx_GM), [1 1 1 Pars.vecSize]);
    glx_WM = reshape(glx.*exp(-FID_times/t2_glx_WM), [1 1 1 Pars.vecSize]);
    glx_CSF = reshape(glx.*exp(-FID_times/t2_glx_CSF), [1 1 1 Pars.vecSize]);
    
    water_GM = reshape(water.*exp(-FID_times/t2_water_GM), [1 1 1 Pars.vecSize]);
    water_WM = reshape(water.*exp(-FID_times/t2_water_WM), [1 1 1 Pars.vecSize]);
    water_CSF = reshape(water.*exp(-FID_times/t2_water_CSF), [1 1 1 Pars.vecSize]);

else
    
    add_DelayTime_ms = 0;
    
    glc_GM = reshape(glc, [1 1 1 Pars.vecSize]);
    glc_WM = reshape(glc, [1 1 1 Pars.vecSize]);
    glc_CSF = reshape(glc, [1 1 1 Pars.vecSize]);
    
    glx_GM = reshape(glx, [1 1 1 Pars.vecSize]);
    glx_WM = reshape(glx, [1 1 1 Pars.vecSize]);
    glx_CSF = reshape(glx, [1 1 1 Pars.vecSize]);
    
    water_GM = reshape(water, [1 1 1 Pars.vecSize]);
    water_WM = reshape(water, [1 1 1 Pars.vecSize]);
    water_CSF = reshape(water, [1 1 1 Pars.vecSize]);
    
end

% Dynamic behaviour of metabolites over 14-63 min (based on in vivo observation)

glx_dyn_WM = linspace(0.25,1,Pars.nRep); 
glx_dyn_CSF = linspace(0.25,1,Pars.nRep); 
glx_dyn_GM = linspace(0.1,1.3,Pars.nRep);
glc_dyn_increase = 1-exp(-(0:Pars.nRep-1:55)/20)+0.0381;
water_dyn_increase = linspace(0.9,1,Pars.nRep);

%% Simulating csi data

% Generating 4D data seperately for each repetition (memory limit) and then cutting out required k-space and saving low resolution data

csi_data_lr.Data = zeros([Pars.TargetImageSize Pars.vecSize Pars.nRep]);
temp_data_hr.Data = zeros([Pars.ImageSize_hr Pars.vecSize]);

% Plotting one voxel as a check of shift applications

for Rep = 1:Pars.nRep
    
    % Additional correction due to different higher resolution -> low resolution scaling compared to previous synthetic data
    
    fprintf('\n Generating 2H-MRSI data repetition %d of 8 \n', Rep);
    
    temp_data_hr.Data = ( 0.7 * (glc_GM .* GM_hr + 0.7 * glc_WM .* WM_hr + glc_CSF .* CSF_hr) * glc_dyn_increase(Rep)); 
        
    temp_data_hr.Data = temp_data_hr.Data + (0.8 * (1 * glx_GM .* GM_hr * glx_dyn_GM(Rep) + 0.25 * glx_WM .* WM_hr * glx_dyn_WM(Rep) +  0.125 * glx_CSF .* CSF_hr * glx_dyn_CSF(Rep))); 
    
    temp_data_hr.Data = temp_data_hr.Data + (2 * (1.25 * water_GM .* GM_hr + water_WM .* WM_hr + 1.5 * water_CSF .* CSF_hr) * water_dyn_increase(Rep));
    
    % Add BO inhomogeneity here
    
    if Settings.B0_Inhomogeneity_flag
        AdditionalInput.B0Map_Hz = B0_Map;
        AdditionalInput.Mask = mask_hr;
        Settings.Printing = 1;
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

    % cut desired k-space
    
    temp_data_lr.Data = temp_data_hr.Data(start_index(1):end_index(1),start_index(2):end_index(2),start_index(3):end_index(3),:); % Check these numbers are right

    for x=1:size(temp_data_lr.Data,1)
        temp_data_lr.Data(x,:,:,:)=temp_data_lr.Data(x,:,:,:)*exp(1i*pi*x/size(temp_data_lr.Data,1));   
    end

    for y=1:size(temp_data_lr.Data,2)
        temp_data_lr.Data(:,y,:,:)=temp_data_lr.Data(:,y,:,:)*exp(1i*pi*y/size(temp_data_lr.Data,2));   
    end

    % Since no density weighting in k-space acquisition
    if Settings.Hamming_flag
        temp_data_lr.Data=HammingFilter(temp_data_lr.Data,[1 2 3],100,'OuterProduct',1);
    end

    % k-space to image space
    temp_data_lr.Data = fftshift(ifft(ifftshift(temp_data_lr.Data,1), [], 1), 1);
    temp_data_lr.Data = fftshift(ifft(ifftshift(temp_data_lr.Data,2), [], 2), 2);
    temp_data_lr.Data = fftshift(ifft(ifftshift(temp_data_lr.Data,3), [], 3), 3);

    % Realigning frequency for low resolution data
   
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

%% Saving

csi_data_lr.Pars = Pars;
csi_data_lr.Settings = Settings;
save('deut_sim.mat', 'csi_data_lr')
clearvars -except csi_data_lr

