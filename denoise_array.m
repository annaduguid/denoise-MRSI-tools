function X = denoise_array(X,method,varargin)

% Input:
% - X       ...     Array to be denoised
% - method	...     Chosen denoising method, options: 'GlobalPS'; 'LocalPS'; 'StackedPS'; 'SpinSVD'
% Optional Name-Value Inputs:
% - 'rank'         ... Number of components to retain (default: 6 for Global/LocalPS, 8 otherwise)
% - 'FIDptsPerCol' ... Number of FID points per column for SpinSVD (default: 8)
% - 'PatchSize'    ... Patch size for Local PS (default: [5 5 5])
% - 'StepSize'     ... Step size for Local PS (default: 1)
%
% Output:
% - X       ...     Denoised Array

% Reading optional inputs
p = inputParser;
addParameter(p, 'rank', []);
addParameter(p, 'FIDptsPerCol', 8);
addParameter(p, 'patch', [5 5 5]);
addParameter(p, 'step', 1);
parse(p, varargin{:});
rank = p.Results.rank;
FIDptsPerCol = p.Results.FIDptsPerCol;
patch = p.Results.patch;
step = p.Results.step;

X = squeeze(X);
data_size = size(X);

if strcmp(method,'GlobalPS')
    
    if isempty(rank)
        rank = 6;
    end
    
    for rep = 1:data_size(5) % erforming denoising of each repetition separately
        Input = reshape(X(:,:,:,:,rep),[prod(data_size(1:3)) data_size(4)]); % Forming Carorati matrix
        [U_MM S_MM V_MM] = svd(Input,'econ'); % SVD of Casorati matrix         
        X(:,:,:,:,rep) = reshape(U_MM(:,1:rank) * (S_MM(1:rank,1:rank) * transpose(conj(V_MM(:,1:rank)))),data_size(1:4)); % Low rank approximation
    end
    
elseif strcmp(method,'LocalPS')
    
    if isempty(rank)
        rank = 6;
    end
    
    M = false(data_size(1:3));
    M(1:step:end, 1:step:end, 1:step:end) = true;
    [ii, jj, kk] = ind2sub(size(M), find(M));
    
    idx_i = 0:min(patch(1)-1, data_size(1)-1);
    idx_j = 0:min(patch(2)-1, data_size(2)-1);
    idx_k = 0:min(patch(3)-1, data_size(3)-1);
    
    navg = zeros([length(ii), data_size(1:3)]);
    
    for rep = 1:data_size(5) % Performing denoising of each repetition separately
        X_tmp = X(:,:,:,:,rep); % current repetition
        X_dn = zeros(data_size(1:4));
        
        for count=1:length(ii)
            idx = ii(count); jdx = jj(count); kdx = kk(count);
            i = mod(idx + idx_i - 1, data_size(1)) + 1;
            j = mod(jdx + idx_j - 1, data_size(2)) + 1;
            k = mod(kdx + idx_k - 1, data_size(3)) + 1;
            patch_data = X_tmp(i,j,k,:);
            [U_MM S_MM V_MM] = svd(reshape(patch_data,[prod(patch) data_size(4)]),'econ'); % SVD of Casorati matrix         
            patch_dn = reshape(U_MM(:,1:rank) * (S_MM(1:rank,1:rank) * transpose(conj(V_MM(:,1:rank)))),[patch data_size(4)]); % Low rank approximation
            X_dn(i, j, k, :) =  X_dn(i, j, k, :) + patch_dn;
            if rep==1
                navg(count, i, j, k) = 1;
            end
        end
        
        if rep==1
            navg_all = squeeze(sum(navg, 1));
        end
        
        X(:,:,:,:,rep)= X_dn./navg_all;
        
    end
    
elseif strcmp(method,'StackedPS')
    
    if isempty(rank)
        rank = 8;
    end
    
    
    X = permute(X, [1 2 3 5 4]);
    Input = reshape(X,[data_size(1)*data_size(2)*data_size(3)*data_size(5) data_size(4)]); % Casorati matrices for all reps stacked along voxels dimension
    [U_MM, S_MM, V_MM] = svd(Input,'econ');
    X = U_MM(:,1:rank) * (S_MM(1:rank,1:rank) * transpose(conj(V_MM(:,1:rank))));
    X = permute(reshape(X,[data_size(1:3) data_size(5) data_size(4)]), [1 2 3 5 4]);
    
elseif strcmp(method,'SpinSVD')
    
    if isempty(rank)
        rank = 8;
    end
    
    if mod(data_size(4),data_size(5)) ~= 0
        error('Input FID size (%d) must be divisible by nReps (%d).',data_size(4),data_size(5))
    end
    
    Input = reshape(X,[data_size(1)*data_size(2)*data_size(3)*FIDptsPerCol data_size(4)*data_size(5)/FIDptsPerCol]);
    [U_MM, S_MM, V_MM] = svd(Input,'econ');
    X = U_MM(:,1:rank) * (S_MM(1:rank,1:rank) * transpose(conj(V_MM(:,1:rank))));
    X = reshape(X,data_size);
    
end
    
end
