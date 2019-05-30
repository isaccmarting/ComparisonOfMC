% -----------------------------------------------------------------------     
%Inputs:
% noise_image:  the noisy image whose noise level requires to be estimated
% patch_size: the predefined size of patches
% 
%Outputs:
% estsigma: Estimated result given by our method
% -----------------------------------------------------------------------
function estimate_noise_std = EstimateNoise(noise_image, patch_size)

%p_out = image2cols(noise_image, patch_size, 3);
% function res = image2cols(noise_image, pSz, step)
[height, width, channel] = size(noise_image); 
numRowPatches = height - patch_size + 1; 
numColPatches = width - patch_size + 1; 
step = 3; 
patches_row_seq = 1:step:numRowPatches; 
patches_col_seq = 1:step:numColPatches; 
%channel = size(noise_image,3); 
if patches_row_seq(end) ~= numRowPatches
    patches_row_seq = [patches_row_seq numRowPatches]; 
end
if patches_col_seq(end) ~= numColPatches
    patches_col_seq = [patches_col_seq numColPatches]; 
end
%sz = length(range_y)*length(range_x);

numPatches = length(patches_row_seq)*length(patches_col_seq); 
patches = zeros(patch_size*patch_size*channel, numPatches); 
iter = 0; 
for i = patches_row_seq
    for j = patches_col_seq
        iter = iter + 1; 
        %p = noise_image(y:y+pSz-1,x:x+pSz-1,:);
        %tmp(:,idx) = p(:);
        image_patch = noise_image(i:i+patch_size-1, j:j+patch_size-1); 
        patches(:, iter) = image_patch(:); 
    end
end

centered_patches = patches - repmat(mean(patches, 2), [1, numPatches]); 
eigen_value = sort(eig(centered_patches * centered_patches' / numPatches), 'ascend'); 
for endNum = size(eigen_value, 1):-1:1
    meanEig = mean(eigen_value(1:endNum)); 
    if(sum(eigen_value(1:endNum)>meanEig) == sum(eigen_value(1:endNum)<meanEig))
        break; 
    end
end
estimate_noise_std = sqrt(meanEig); 
