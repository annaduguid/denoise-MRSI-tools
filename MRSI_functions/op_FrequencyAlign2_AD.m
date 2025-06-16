function [MRStruct,AdditionalOut] = op_FrequencyAlign2(MRStruct,Settings,AdditionalIn)
% FrequencyAlignment Align frequencies of csi spectra.
%
% This function was written by Bernhard Strasser, April 2015.
%
%
% The function masks the data in k-space, so that k-space values outside of an ellipsoid-shaped mask are set to zero. The mask can be a
% 3d-ellipsoid, or an 2d-ellipse. The equation for the mask is
% mask = {(x,y,z) E R³ | (x/a)² + (y/b)² + (z/c)² <= R²}
% a, b, c, and R can be chosen by the user.
%
%
% [MRStruct,mask] = FrequencyAlignment(MRStruct,ApplyAlongDims,EllipsoidCoefficients,PerformFFT_flag)
%
% Input: 
% -     MRStruct                     ...    Input array to which the filter should be applied
% -     ApplyAlongDims               ...    Along these dimensions the filter is applied. If this vector has two elements, a two dimensional 
%                                          Filter is applied. Otherwise, a 3d filter is used.
% -     EllipsoidCoefficients        ...    The values for [a b c R], which determine the shape and size of the ellipsoid. For two dimensional
%                                          Filter, set c = 1;
% -     PerformFFT_flag              ...    If it is 0, the image gets Fourier transformed to k-space before applying the filter, 
%                                          and transformed back to image domain afterwards
%
% Output:
% -     MRStruct                     ...     The filtered/masked output array
% -     mask                         ...     The mask of the filter
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: myrepmat_1_0

% Further remarks: 



%% 0. Declarations, Preparations, Definitions


% 0.1 Preparations
if(nargin < 1)
	MRStruct = 0;
	return
end
% if(nargin < 2)
%     return
% end
if(~isfield(MRStruct,'Par'))
    MRStruct.Par = struct;
    if(isfield(MRStruct,'Data_file'))
        MRStruct.Par = read_ascconv(MRStruct.Data_file); 
    end
end
if(~isfield(MRStruct,'RecoPar'))
    MRStruct.RecoPar = MRStruct.Par; % Mine wont have this .RecoPar so will be assigned as same as .Par
end

size_InArray = size(MRStruct.Data); % Size of csi data

if(~exist('Settings','var'))
    Settings = struct;
end
if(~isfield(Settings,'ApplyAlongDim'))
    Settings.ApplyAlongDim = find(size_InArray == max(size_InArray)); % This assumes that the maximum dimension is the spectral dimension so that we know which dimension to apply shifts along
end
if(Settings.ApplyAlongDim ~= 4)
    error('ApplyAlongDim has to be 4 currently, but it is %d',Settings.ApplyAlongDim); % But dimension must always be for (3 spatial dimensions, one spectral)
end
if(~isfield(Settings,'ZerofillingFactor'))
    Settings.ZerofillingFactor = 2; % We dont want zero filling so should set this to 1 (I dont think it will do anything anyway)
end
if(~isfield(Settings,'PeakSearchPPM'))
    Settings.PeakSearchPPM = 3.02;  % Default: Cr . We are not using peak searching so should be fine to leave this
end
if(~isfield(Settings,'PeakSearchRangePPM')) % Not used when only applying shift
    Settings.PeakSearchRangePPM = 0.2;
end
if(~isfield(Settings,'UseSVDForRef_flag')) % Not used when only applying shift
    Settings.UseSVDForRef_flag = true;
end
if(~isfield(Settings,'UseSVDForRef_NoSingVals')) % Not used when only applying shift
    Settings.UseSVDForRef_NoSingVals = 1; 
end
if(~isfield(Settings,'MaxShifts_ppm')) % Not used when only applying shift
    Settings.MaxShifts_ppm = [-0.1 0.1];     % to left and right direction
end
if(~isfield(Settings,'AlignRefPeak_flag')) % Not used when only applying shift
    Settings.AlignRefPeak_flag = 1;
end
if(~isfield(Settings,'PeakDetection')) % Not used when only applying shift
    Settings.PeakDetection = struct;
end
if(~isfield(Settings.PeakDetection,'MinSNR')) % Not used when only applying shift
    Settings.PeakDetection.MinSNR = 2;
end
if(~isfield(Settings.PeakDetection,'MinDotProd')) % Not used when only applying shift
    Settings.PeakDetection.MinDotProd = 0.55;
end
if(~isfield(Settings,'FindShiftAlgorithm'))
    Settings.FindShiftAlgorithm = 'DotProduct'; % Not used when only applying shift
end
if(~isfield(Settings,'UseThisInStructMask'))
    Settings.UseThisInStructMask = 'BrainMask'; % Not compulsory for masking can use other arguments
end
if(~exist('AdditionalIn','var'))
    AdditionalIn = struct();
end
if(isfield(Settings,'Printing')) % Not used when only applying shift
    Settings.Printing_flag = true;
else
    Settings.Printing_flag = 0;
end

if(~exist('AdditionalIn','var') || ~isfield(AdditionalIn,'RefSpecIn')) % This will assign RefSpec as empty array because I dont provide one
	RefSpec = [];
else
    RefSpec = reshape(AdditionalIn.RefSpecIn,[ones([1 Settings.ApplyAlongDim-1]) numel(AdditionalIn.RefSpecIn)]);
    RefVox = [-1 -1 -1];
end

% Handle mask
MaskWasProvided = true; % Can either provide in MRStruct as 'Mask' or 'BrainMask' or AdditionalIn as 'Mask'
if(isfield(AdditionalIn,'Mask'))
    mask = AdditionalIn.Mask;
elseif(isfield(Settings,'UseThisInStructMask') && isfield(MRStruct,(Settings.UseThisInStructMask))) % If don't provide argument in here then it will skip
    mask = MRStruct.(Settings.UseThisInStructMask);
elseif(isfield(MRStruct,'BrainMask'))
	mask = MRStruct.BrainMask;
elseif(isfield(MRStruct,'Mask')) % Make sure I've capitalised Mask so that I'm only applying shift within the mask - much less computationally expensive
	mask = MRStruct.Mask;    
else
	mask = ones(size_InArray(1:3));     % The mask is always only spatial
	MaskWasProvided = false;	
	error('Mask not read in - would have to apply shifts to whole array and takes ages');
end


% Assume B0Map is in Hz
if(exist('AdditionalIn','var') && isfield(AdditionalIn,'B0Map'))
    AdditionalIn.B0Map_Hz = AdditionalIn.B0Map; 
    AdditionalIn = rmfield(AdditionalIn,'B0Map');
end
if(exist('AdditionalIn','var') && isfield(AdditionalIn,'B0Map_Hz'))
    HzPerPt = 1E9/MRStruct.RecoPar.Dwelltimes(1)/MRStruct.RecoPar.vecSize; 
    AdditionalIn.ShiftMap = AdditionalIn.B0Map_Hz/ HzPerPt; % Why do we divide by Hz per point?? Do I want to do this???
end

if(exist('AdditionalIn','var') && isfield(AdditionalIn,'ShiftMap'))
    AdditionalOut.ShiftMap = AdditionalIn.ShiftMap;
	OnlyApplyShiftMap_flag = true;
    RefSpecWasProvided = false;
    Settings.AlignRefPeak_flag = false;
else
	OnlyApplyShiftMap_flag = false;
    disp('Applying frequency alignment')
% 	error('Running extra code rather than only applying shift map');

end


%% Define Reference Spectrum if necessary

if(~OnlyApplyShiftMap_flag)
    RefSpecWasProvided = true;
    if(isempty(RefSpec))
        RefSpecWasProvided = false;
        
        if(Settings.UseSVDForRef_flag)
            if(Settings.Printing_flag);disp('Using SVD for Ref Spec');end
            CurData = MRStruct.Data .* mask;    % We don't want to take all the lipids etc.
            [~, ~, V] = svd(reshape(CurData,[numel(CurData)/size(CurData,4) size(CurData,4)]),'econ');
            V = V';
            RefSpec = reshape(sum(V(1:Settings.UseSVDForRef_NoSingVals,:),1),[ones([1 Settings.ApplyAlongDim-1]) size(V,1)]); RefVox = [0 0 0];
            clear CurData
            if(isempty(RefSpec))
                disp('WARNING: RefSpec Empty')
            end
            if(~isempty(RefSpec))
                if(Settings.Printing_flag);disp('Created Ref Spec using SVD');end
            end
        else
        
            if(isfield(Settings,'RefVox') && numel(Settings.RefVox) > 1)
                RefVox = Settings.RefVox;
            elseif(MaskWasProvided)
                % Define Reference Voxel as Center of Mass Voxel
                RefVox = regionprops(mask, 'centroid');
                RefVox = round(RefVox.Centroid);
            else
                RefVox = floor(size_InArray/2)+1;
                fprintf('\n\nWARNING: No mask input for FrequencyAlignment. Reference voxel will be set as (%d,%d,%d). Might be wrong!',RefVox(1),RefVox(2),RefVox(3))
            end

            if(numel(RefVox) < 3)
                RefVox(3) = 1;
            end

            % Set RefSpec
            RefSpec = MRStruct.Data(RefVox(1),RefVox(2),RefVox(3),:);
        end
    end
end





% 0.2 Declarations


% 0.3 Definitions
    
size_SearchArray = size_InArray; size_SearchArray(Settings.ApplyAlongDim) = size_SearchArray(Settings.ApplyAlongDim)*Settings.ZerofillingFactor;
SearchArray = MRStruct.Data;

if(~OnlyApplyShiftMap_flag)
    CS_vec_zf = compute_chemshift_vector(MRStruct.RecoPar.LarmorFreq,MRStruct.RecoPar.Dwelltimes(1)/10^9,MRStruct.RecoPar.vecSize*Settings.ZerofillingFactor);
    SearchForPeak_Center_Pts = find(min(abs(CS_vec_zf - Settings.PeakSearchPPM)) == abs(CS_vec_zf - Settings.PeakSearchPPM));
  	SearchForPeak_LeftPt_Pts = find(min(abs(CS_vec_zf - Settings.PeakSearchPPM - Settings.PeakSearchRangePPM)) == abs(CS_vec_zf - Settings.PeakSearchPPM - Settings.PeakSearchRangePPM));
	SearchForPeak_RightPt_Pts = find(min(abs(CS_vec_zf - Settings.PeakSearchPPM + Settings.PeakSearchRangePPM)) == abs(CS_vec_zf - Settings.PeakSearchPPM + Settings.PeakSearchRangePPM));
end


%% 1. Zerofilling & FFT SearchArray

if(~OnlyApplyShiftMap_flag)
    SearchArray = Zerofilling_Spectral(SearchArray,size_SearchArray,0);
    SearchArray = fftshift(fft(SearchArray,[],Settings.ApplyAlongDim),Settings.ApplyAlongDim);
    RefSpec2 = Zerofilling_Spectral(RefSpec,[ones([1 numel(size_SearchArray)-1]) size_SearchArray(Settings.ApplyAlongDim)],0);
    RefSpec2 = fftshift(fft(RefSpec2,[],Settings.ApplyAlongDim),Settings.ApplyAlongDim); % Time to frequency domain
%     figure;plot(abs(squeeze(RefSpec2(1,1,1,:))));title('RefSpec freqdomain')
    if(isempty(RefSpec2))
        disp('WARNING: RefSpec2 Empty')
    end
end


%% Align Ref Spectrum


if(Settings.AlignRefPeak_flag)
    % Find closest local maximum to Settings.PeakSearchPPM
    
    % Find all peaks with prominence higher than 1
    [PeakHght,PeakPos] = findpeaks(abs(squeeze(RefSpec2))./max(abs(RefSpec2(SearchForPeak_LeftPt_Pts:SearchForPeak_RightPt_Pts))),'MinPeakProminence',0.3);
%     figure; plot(abs(squeeze(RefSpec2))); hold on, scatter(PeakPos,PeakHght,'r'), hold off        % (Debug)
    if(Settings.Printing_flag);fprintf('PeakHeightRefSpec %f \n\n\n', PeakHght);end
    if(Settings.Printing_flag);fprintf('PeakPosnRefSpec  %f \n\n\n', PeakPos);end
    % Restrict peaks to those inside of the Search-range
    DeletePeaks = PeakPos < SearchForPeak_LeftPt_Pts | PeakPos > SearchForPeak_RightPt_Pts;
    PeakPos(DeletePeaks) = []; PeakHght(DeletePeaks) = []; 
    
    % Use the one that is closest to SearchForPeak_Center_Pts
    tmp = min(abs(PeakPos - SearchForPeak_Center_Pts)) == abs(PeakPos - SearchForPeak_Center_Pts);
    PeakPos = PeakPos(tmp); PeakHght = PeakHght(tmp);
    
    % Circshift RefSpectrum
    RefSpec2 = circshift(RefSpec2,SearchForPeak_Center_Pts-PeakPos,Settings.ApplyAlongDim);
    
    % Debug:
%     figure; plot(abs(squeeze(RefSpec2))); hold on, scatter(PeakPos,PeakHght,'r'), hold off        % (Debug)
    clear tmp PeakPos PeakHght DeletePeaks
end



%% 2. Create Matrix of Differently Shifted RefSpecs
% Instead of shifting the individual spectra of all voxels N times, computing the DotProduct and finding the maximum,
% do instead: Shift the RefSpec once N times, save that, and calculate the DotProducts for all shifts of the RefSpec.
% In the fist case, we would need to do NumbOfVox x N shifts, whereas in the latter we only need N shifts.
if(~OnlyApplyShiftMap_flag)
	
	AdditionalOut.ShiftMap = zeros([size(MRStruct.Data,1) size(MRStruct.Data,2) size(MRStruct.Data,3)]);
    AdditionalOut.DotProdMap = AdditionalOut.ShiftMap;

    MaxShifts_Pts = round(Settings.MaxShifts_ppm * (MRStruct.RecoPar.LarmorFreq/1E6)* MRStruct.RecoPar.vecSize/(1E9/ MRStruct.RecoPar.Dwelltimes(1)));
    
    if(Settings.Printing_flag);fprintf('MaxShifts_Pts %f \n\n\n', MaxShifts_Pts);end
    
	ReferenceSpecMat_Spec = squeeze(RefSpec2(1,1,1,SearchForPeak_LeftPt_Pts:SearchForPeak_RightPt_Pts)); 
    ReferenceSpecMat_Spec = ReferenceSpecMat_Spec/norm(ReferenceSpecMat_Spec);
	ReferenceSpecMat = zeros(size(ReferenceSpecMat_Spec,1));
    
    if(Settings.Printing_flag);fprintf('SearchForPeak RightPt Pts %d \n\n\n', SearchForPeak_RightPt_Pts);end
    if(Settings.Printing_flag);fprintf('SearchForPeak LeftPt Pts %d \n\n\n', SearchForPeak_LeftPt_Pts);end

	CircShiftVec = SearchForPeak_LeftPt_Pts-SearchForPeak_Center_Pts :1: SearchForPeak_RightPt_Pts-SearchForPeak_Center_Pts;
    if(Settings.Printing_flag);str = sprintf('%g ',CircShiftVec);str = strtrim(str);fprintf('CircShiftVec: %s\n\n\n', str);end
    
	for i = 1:abs(SearchForPeak_RightPt_Pts-SearchForPeak_LeftPt_Pts+1)
		ReferenceSpecMat(i,:) = circshift(ReferenceSpecMat_Spec,CircShiftVec(i)); % Shifting normalised ref spec mat spec by different amounts
	end
	
end


%% 3. Calculate & Apply ShiftMap
ShiftsIncluded = [];
% MRStruct.Data = fftshift(fft(MRStruct.Data,[],Settings.ApplyAlongDim),Settings.ApplyAlongDim);
if(~OnlyApplyShiftMap_flag)
    for x = 1:size(MRStruct.Data,1)
        for y = 1:size(MRStruct.Data,2)
            for z = 1:size(MRStruct.Data,3)

                if(mask(x,y,z) == 0 || (RefSpecWasProvided && x==RefVox(1) && y == RefVox(2) && z == RefVox(3)))			% Dont process if mask=0 or reference voxel
                    continue
                end

    %             if(x==23 && y == 24)
    %                 sadf = 1; 
    %             end

                TmpSpec = squeeze(SearchArray(x,y,z,SearchForPeak_LeftPt_Pts:SearchForPeak_RightPt_Pts)); % This is the spectrum for voxel x,y,z within the peak search range
    %             TmpSpec = TmpSpec - fftshift(fft(exp(-transpose(0:numel(TmpSpec)-1)/1) .* ifft(ifftshift(TmpSpec))));
    %             TmpSpec = abs(TmpSpec);
                TmpSpec = conj(TmpSpec) / norm(TmpSpec);
                DotProd = abs(ReferenceSpecMat * TmpSpec);
                AdditionalOut.DotProdMap(x,y,z) = max(DotProd);
                if(max(abs(TmpSpec)) > Settings.PeakDetection.MinSNR*std(TmpSpec)) % If peak signal is greater than min signal from SNR and std dev noise
                    % Calculate ShiftMap
                    if(strcmpi(Settings.FindShiftAlgorithm,'DotProduct'))
                        MaxDotProd =  max(DotProd);
                        MaxDotProdMatch = find(DotProd == MaxDotProd); MaxDotProdMatch = MaxDotProdMatch(1); % finding location of maximum dot product
                        if x==10 && y ==10 && z == 10
                                if(Settings.Printing_flag);str = sprintf('%g ', DotProd);str = strtrim(str);fprintf('DotProd: %s\n\n\n', str);end
                                if(Settings.Printing_flag);fprintf('MaxDotProd %g \n\n\n', MaxDotProd);end
                                if(Settings.Printing_flag);fprintf('MaxDotProdMatch %g \n\n\n', MaxDotProdMatch);end
                        end
                        ShiftPoint = -( CircShiftVec(MaxDotProdMatch) / Settings.ZerofillingFactor);
                    else
                        [~, MaxPos] = max(abs(abs(TmpSpec))); 
                        ShiftPoint = (ceil(numel(TmpSpec)/2) - MaxPos )/ Settings.ZerofillingFactor; % Max should be always in center of TmpSpec!
                        MaxDotProd = 1E9;    % So that we have no condition on DotProd
                    end
                    if(MaxDotProd > Settings.PeakDetection.MinDotProd && ShiftPoint > min(MaxShifts_Pts) && ShiftPoint < max(MaxShifts_Pts)) % 
                        AdditionalOut.ShiftMap(x,y,z) = ShiftPoint; % - because we shifted the reference, but below we want to shift the other spectra
                        ShiftsIncluded = [ShiftsIncluded, ShiftPoint];
                        
                    else
                        AdditionalOut.ShiftMap(x,y,z) = NaN; continue;
                    end
                end

                % Apply ShiftMap
    % 			MRStruct.Data(x,y,z,:) = circshift(squeeze(MRStruct.Data(x,y,z,:)),[AdditionalOut.ShiftMap(x,y,z) 1]);

            end
        end
    end
end

% if(Settings.Printing_flag);fprintf('Shifts included %d \n\n\n', ShiftsIncluded);end

% MRStruct.Data = ifft(fftshift(MRStruct.Data,Settings.ApplyAlongDim),[],Settings.ApplyAlongDim);
CurMap = AdditionalOut.ShiftMap; CurMap(isnan(CurMap)) = 0;         % Can't multiply with NaN, otherwise data gets NaNed
t = (0:(MRStruct.RecoPar.vecSize-1))/MRStruct.RecoPar.vecSize;

% Adding acquisition delay...
if(isfield(AdditionalIn,'B0Map_Hz'))
    t = (t * MRStruct.RecoPar.Dwelltimes(1)*1E-9 * MRStruct.RecoPar.vecSize)+ 2E-3;
    CurMap = CurMap*HzPerPt;
end

MRStruct.Data = MRStruct.Data .* exp(1i*2*pi*CurMap.*reshape(t,[ones([1 Settings.ApplyAlongDim-1]) numel(t) ones([1 ndims(MRStruct.Data)-Settings.ApplyAlongDim])]));

if(nargout > 1)
    HzPerPt = 1E9/MRStruct.RecoPar.Dwelltimes(1)/MRStruct.RecoPar.vecSize;
    AdditionalOut.B0Map_Hz = AdditionalOut.ShiftMap * HzPerPt;
end

%% 5. Postparations

supp_UpdateRecoSteps(MRStruct,Settings);




