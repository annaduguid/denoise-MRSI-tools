%% Denoise array with different methods (Global PS, Local PS, Stacked PS, SPIN-SVD)

clearvars;

% Functions and toolboxes
% addpath(genpath('/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/aduguid/Matlab/MatlabFunctions'));
addpath('/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/aduguid/Matlab/LowRankDenoising/WIP_ForGithub/Functions/');

% Loading lesion simulation
load(strcat(pwd,'/deut_lesion_sim.mat'));

%% Add noise in  k-space then hamming filter...

csi_data_lr.Data = fftshift(fft(ifftshift(csi_data_lr.Data,1), [], 1), 1);
csi_data_lr.Data = fftshift(fft(ifftshift(csi_data_lr.Data,2), [], 2), 2);
csi_data_lr.Data = fftshift(fft(ifftshift(csi_data_lr.Data,3), [], 3), 3);

rng(1);
Noisesigma = 1;

noise=(randn(size(csi_data_lr.Data))+1i*randn(size(csi_data_lr.Data)))*Noisesigma;

csi_data_lr.Data=csi_data_lr.Data + noise;

csi_data_lr.Data=HammingFilter(csi_data_lr.Data,[1 2 3],100,'OuterProduct',1);

csi_data_lr.Data = fftshift(ifft(ifftshift(csi_data_lr.Data,1), [], 1), 1);
csi_data_lr.Data = fftshift(ifft(ifftshift(csi_data_lr.Data,2), [], 2), 2);
csi_data_lr.Data = fftshift(ifft(ifftshift(csi_data_lr.Data,3), [], 3), 3);
    

%% Definitions

% MetaboInfo.Water.PPMRange = [3.85 5.7];%4.8
% MetaboInfo.Glc.PPMRange = [2.8 3.85];%3.9
% MetaboInfo.Glx.PPMRange = [1.6 2.4];%2.4
% MetaboInfo.Lac.PPMRange = [0.2 1.6];%1.4

MetaboInfo.Water.PPMRange = [3.85 5.7];%4.8
MetaboInfo.Glc.PPMRange = [2.8 3.85];%3.9
MetaboInfo.Glx.PPMRange = [2.0 2.8];%2.4
MetaboInfo.Lac.PPMRange = [0.9 1.5];%1.3

%% Show non-denoised results

MetaboMap = CalcMetaboMap(csi_data_lr,[],MetaboInfo);
PlotMetaboMaps(MetaboMap,15,8) % slice 15, Time Pt 8

%% LocalPS

csi_dn = csi_data_lr;

% With default rank, patch, step
csi_dn.Data = denoise_array(csi_data_lr.Data,'LocalPS');

MetaboMap = CalcMetaboMap(csi_dn,[],MetaboInfo);
PlotMetaboMaps(MetaboMap,15,8) % plotting slice 15, 2H-MRSI repetition 8

% With smaller patch
csi_dn.Data = denoise_array(csi_data_lr.Data,'LocalPS','patch',[4 4 4]);

MetaboMap = CalcMetaboMap(csi_dn,[],MetaboInfo);
PlotMetaboMaps(MetaboMap,15,8) % slice 15, Time Pt 8

% With even smaller patch
csi_dn.Data = denoise_array(csi_data_lr.Data,'LocalPS','patch',[3 3 3]);

MetaboMap = CalcMetaboMap(csi_dn,[],MetaboInfo);
PlotMetaboMaps(MetaboMap,15,8) % slice 15, Time Pt 8


%% GlobalPS

csi_dn = csi;

% With default rank
csi_dn.Data = denoise_array(csi_data_lr.Data,'GlobalPS');

MetaboMap = CalcMetaboMap(csi_dn,[],MetaboInfo);
PlotMetaboMaps(MetaboMap,15,8) % slice 15, Time Pt 8

% With higher rank
csi_dn.Data = denoise_array(csi_data_lr.Data,'GlobalPS','rank',40);

MetaboMap = CalcMetaboMap(csi_dn,[],MetaboInfo);
PlotMetaboMaps(MetaboMap,15,8) % slice 15, Time Pt 8

%% StackedPS

% With default rank
csi_dn.Data = denoise_array(csi_data_lr.Data,'StackedPS');

MetaboMap = CalcMetaboMap(csi_dn,[],MetaboInfo);
PlotMetaboMaps(MetaboMap,15,8) % slice 15, Time Pt 8

% With higher rank
csi_dn.Data = denoise_array(csi_data_lr.Data,'StackedPS','rank',40);

MetaboMap = CalcMetaboMap(csi_dn,[],MetaboInfo);
PlotMetaboMaps(MetaboMap,15,8) % slice 15, Time Pt 8

%% SpinSVD

% With default settings
csi_dn.Data = denoise_array(csi_data_lr.Data,'SpinSVD');

MetaboMap = CalcMetaboMap(csi_dn,[],MetaboInfo);
PlotMetaboMaps(MetaboMap,15,8) % slice 15, Time Pt 8

% With higher rank
csi_dn.Data = denoise_array(csi_data_lr.Data,'SpinSVD','rank',40);

MetaboMap = CalcMetaboMap(csi_dn,[],MetaboInfo);
PlotMetaboMaps(MetaboMap,15,8) % slice 15, Time Pt 8

% With higher FIDptsperCol
csi_dn.Data = denoise_array(csi_data_lr.Data,'SpinSVD','FIDptsPerCol',24);

MetaboMap = CalcMetaboMap(csi_dn,[],MetaboInfo);
PlotMetaboMaps(MetaboMap,15,8) % slice 15, Time Pt 8


%% Function for calculating metabolic maps from signal

function MetaboMap = CalcMetaboMap(Signal,Mask,MetaboInfo)

    if(~exist('Mask','var') || isempty(Mask))
        Sizze = size(Signal.Data); % Size function to determine dimentionality of csi data
        Mask = ones([Sizze(1:3) 1 Sizze(5:end)]); % Mask is a matrix of ones size 22 x 22 x 21 x 1 x 8
    
    end
    ppmvec = compute_chemshift_vector(Signal.Par);
    Signal_Spec = fftshift(fft(Signal.Data,[],4),4); % Fourier transform performed along the fourth dimension. For a vector fft shift swaps the left and right halevs of x.
    Fields = fieldnames(MetaboInfo); % Lists all of the defined metabolites e.g. water, glc, glx
    
    for ii = 1:numel(Fields) % For each metabolite. Numel counts the number of array elements.
        CurField = Fields{ii}; % The name of the current metabolite
        Tmp = FindClosestIndex(ppmvec,MetaboInfo.(CurField).PPMRange(1)); SumFromTo(1) = Tmp{1}; %finding the index of ppmvec that is closest to the ppm range
        Tmp = FindClosestIndex(ppmvec,MetaboInfo.(CurField).PPMRange(2)); SumFromTo(2) = Tmp{1}; % " "
        SumFromTo = sort(SumFromTo);
        MetaboMap.(CurField) = squeeze_single_dim(sum(abs(Signal_Spec(:,:,:,SumFromTo(1):SumFromTo(2),:,:)),4) .* Mask,4); % ???
    end
    
end

%% Function for plotting metabolic maps from signal

function [] = PlotMetaboMaps(MetaboMap,PlotSlc,PlotTimePt) % Function plots metabolic map for given slice (PlotSlc) and time point 

    Fields = fieldnames(MetaboMap);
    if(~exist('PlotSlc','var')) % Tilde is a logical negotiation so the if will be true if 
        PlotSlc = round(size(MetaboMap.(Fields{1}),3)/2); % Plots middle slice if no slice has been selected
    end
    if(~exist('PlotTimePt','var'))
        PlotTimePt = round(size(MetaboMap.(Fields{1}),4)/2); % Plots middle time point if none has been selected
    end    
    
    NoOfSubplots = ceil(sqrt(numel(Fields))); % ceil(X) rounds each element of X to the nearest integer greater than or equal to that element.
    figure;
    for ii = 1:numel(Fields)
        CurField = Fields{ii};
        subplot(NoOfSubplots,NoOfSubplots,ii);
        imagesc(MetaboMap.(CurField)(:,:,PlotSlc,PlotTimePt))
        axis square
        title(CurField)
    end
    
end