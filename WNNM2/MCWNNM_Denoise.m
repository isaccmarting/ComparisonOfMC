function [estimate_image, Parameters]   =  MCWNNM_Denoise(noise_image, ref_image, Parameters)
[height, width, channel] = size(noise_image); 

% Precompute the all the patch indexes in the Searching window
numRowPatches = height - Parameters.patch_size + 1; 
numColPatches = width - Parameters.patch_size + 1; 
keyPatches_row_seq = 1:Parameters.step:numRowPatches; 
keyPatches_row_seq = [keyPatches_row_seq keyPatches_row_seq(end)+1:numRowPatches]; 
keyPatches_col_seq =  1:Parameters.step:numColPatches; 
keyPatches_col_seq =  [keyPatches_col_seq keyPatches_col_seq(end)+1:numColPatches]; 
% Total number of patches in the test image
numPatches = numRowPatches * numColPatches; 
% Total number of seed patches being processed
numKeyPatches = length(keyPatches_row_seq) * length(keyPatches_col_seq); 
% index of each patch in image 
patches_index = reshape(1:numPatches, [numRowPatches numColPatches]); 
% preset variables for all the patch indexes in the Searching window 
keyPatches_neighborsIndex = int32(zeros((2*Parameters.window_size+1)^2, numKeyPatches)); 
numNeighbors = int32(zeros(1, numKeyPatches)); 
keyPatchIndex2PatchIndex = int32(zeros(1, numKeyPatches)); 
for i = 1:length(keyPatches_row_seq)
    for j = 1:length(keyPatches_col_seq)
        patch_row_index = keyPatches_row_seq(i); 
        patch_col_index = keyPatches_col_seq(j); 
        patch_index = (patch_col_index-1) * numRowPatches + patch_row_index; 
        keyPatch_index = (j-1) * length(keyPatches_row_seq) + i; 
        
        top = max(patch_row_index - Parameters.window_size, 1); 
        bottom = min(patch_row_index + Parameters.window_size, numRowPatches); 
        left = max(patch_col_index - Parameters.window_size, 1); 
        right = min(patch_col_index + Parameters.window_size, numColPatches); 
        neigborsIndex = patches_index(top:bottom, left:right); 
        neigborsIndex = neigborsIndex(:); 

        numNeighbors(keyPatch_index) = length(neigborsIndex); 
        keyPatches_neighborsIndex(1:numNeighbors(keyPatch_index), keyPatch_index) = neigborsIndex; 
        keyPatchIndex2PatchIndex(keyPatch_index) = patch_index; 
    end
end

% noisy image to patch
noise_patches =	Image2Patches(noise_image, Parameters.patch_size, channel, numPatches); 
patch_size_square = Parameters.patch_size * Parameters.patch_size; 
tempW = zeros(channel*patch_size_square, numPatches);
numSimilarPatches = Parameters.ori_numSimilarPatches;   % Initial Non-local Patch number 
estimate_image = noise_image;   % Estimated Image 
for k = 1 : Parameters.K2
    % iterative regularization
    estimate_image = estimate_image + Parameters.delta * (noise_image - estimate_image); 
    % image to patch
    estimate_patches = Image2Patches(estimate_image, Parameters.patch_size, channel, numPatches); 
    % estimate local noise variance
    for ch = 1:channel
        if(k == 1)
            sigma_channel = repmat(Parameters.channels_noise_std(ch), [1 size(estimate_patches, 2)]); 
        else
            sigma_channel = Parameters.lambda*sqrt(max(0, repmat(Parameters.channels_noise_std(ch)^2, 1, size(estimate_patches, 2)) - mean((noise_patches((ch-1)*patch_size_square+1:ch*patch_size_square, :) - estimate_patches((ch-1)*patch_size_square+1:ch*patch_size_square, :)).^2))); 
        end
        tempW((ch-1)*patch_size_square+1:ch*patch_size_square, :) = repmat(sigma_channel, [patch_size_square, 1]);
    end
    if (mod(k-1, Parameters.reBlkMacLoop) == 0)
        numSimilarPatches = numSimilarPatches - 10;  % Lower Noise level, less NL patches
        % record the indexs of patches similar to the seed patch
	    keyPatches_SimilarPatches = zeros(numSimilarPatches, numKeyPatches, 'single'); 
        for keyPatchIndex = 1:numKeyPatches
            keyPatch = estimate_patches(:, keyPatchIndex2PatchIndex(keyPatchIndex)); 
            neighbors = estimate_patches(:, keyPatches_neighborsIndex(1:numNeighbors(keyPatchIndex), keyPatchIndex)); 
            distance = sum(bsxfun(@minus, neighbors, keyPatch).^2, 1); 
            [~, sort_indices] = sort(distance); 
            keyPatches_SimilarPatches(:, keyPatchIndex) = keyPatches_neighborsIndex(sort_indices(1:numSimilarPatches), keyPatchIndex); 
        end
    end
    % Denoising by MCWNNM
    estimate_Xs = zeros(size(estimate_patches)); 
    patches_weight = zeros(size(estimate_patches)); 
    for keyPatchIndex = 1:numKeyPatches % For each keypatch group
        Y = estimate_patches(:, keyPatches_SimilarPatches(1:numSimilarPatches, keyPatchIndex)); % Non-local similar patches to the keypatch
        meanY = repmat(mean(Y, 2), 1, numSimilarPatches); 
        Y = Y - meanY; 

        temp = tempW(:, keyPatchIndex2PatchIndex(keyPatchIndex)); 
        minTemp = min(temp); 
        W = (minTemp+eps) ./ (temp+eps); 
        C = (2 * sqrt(2*numSimilarPatches) * minTemp^2); 
        X = MCWNNM_ADMM(Y, W, C, Parameters.rho, Parameters.mu, Parameters.K1) + meanY; % WNNM Estimation
        estimate_Xs(:, keyPatches_SimilarPatches(1:numSimilarPatches, keyPatchIndex)) = estimate_Xs(:, keyPatches_SimilarPatches(1:numSimilarPatches, keyPatchIndex)) + X; 
        patches_weight(:, keyPatches_SimilarPatches(1:numSimilarPatches, keyPatchIndex)) = patches_weight(:, keyPatches_SimilarPatches(1:numSimilarPatches, keyPatchIndex)) + ones(channel*patch_size_square, numSimilarPatches); 
    end
    estimate_image = Patches2Image(estimate_Xs, patches_weight, height, width, channel, Parameters.patch_size); 
    estimate_image(estimate_image > 255) = 255; 
    estimate_image(estimate_image < 0) = 0; 
    if Parameters.bEvaluate
        Parameters.PSNR(k, Parameters.image_index) = psnr_MultiCh(estimate_image, ref_image); 
        Parameters.SSIM(k, Parameters.image_index) = ssim(estimate_image, ref_image); 
        fprintf('Iter = %2.3f, PSNR = %2.4f, SSIM = %2.4f \n', k, Parameters.PSNR(k, Parameters.image_index), Parameters.SSIM(k, Parameters.image_index)); 
    end
end
