%-------------------------------------------------------------------------------------------------------------
% This is an implementation of the MCWNNM algorithm for real color image denoising. 
%-------------------------------------------------------------------------------------------------------------

clear;
image_dir = 'kodak_color/'; 
%image_dir = 'Real_ccnoise_denoised_part/'; 
%image_dir = 'NoiseClinicImages/'; 
% denoising_type: CleanImage, NoiseImageWithReference, NoiseImageWithoutReference 
denoising_type = 'CleanImage'; 
% Changing to 'NoiseImageWithReference', you need change Parameters and sigma_channel. 
%denoising_type = 'NoiseImageWithReference'; 
% Changing to 'NoiseImageWithReference', you need change Parameters and sigma_channel. 
%denoising_type = 'NoiseImageWithoutReference'; 
[Parameters, ref_images, noise_images] = Initialize(denoising_type, image_dir); 

if Parameters.output_noise_image
    noise_image_dir = ['noise_' image_dir]; 
    if ~mkdir(noise_image_dir)
        fprintf(['Cannot make directory ' noise_image_dir '! \n']); 
        return; 
    end
end
if Parameters.output_denoised_image
    denoised_image_dir = ['denoised_' image_dir]; 
    if ~mkdir(denoised_image_dir)
        fprintf(['Cannot make directory ' denoised_image_dir '! \n']); 
        return; 
    end
end

if strcmp(denoising_type, 'CleanImage')
    numImages = length(ref_images); 
else
    numImages = length(noise_images); 
end
% record all the results in each iteration 
Parameters.PSNR = zeros(Parameters.K2, numImages, 'single'); 
Parameters.SSIM = zeros(Parameters.K2, numImages, 'single'); 
for i = 1:numImages
    if strcmp(denoising_type, 'CleanImage')
        fprintf('%s :\n', ref_images(i).name); 
        ref_image = double(imread(fullfile(['../' image_dir], ref_images(i).name))); 
        [height, width, channel] = size(ref_image); 
        noise_image = zeros([height, width, channel]); 
        for ch = 1:channel
            randn('seed', 0); 
            noise_image(:, :, ch) = ref_image(:, :, ch) + Parameters.channels_noise_std(ch) * randn([height, width]); 
        end
        if Parameters.output_noise_image
            imwrite(noise_image/255, [noise_image_dir ref_images(i).name]); 
            fprintf(['The noise image is written to file ' noise_image_dir ref_images(i).name '. \n']); 
        end
    else
        fprintf('%s :\n', noise_images(i).name); 
        noise_image = double(imread(fullfile(['../' image_dir], noise_images(i).name))); 
        [height, width, channel] = size(noise_image); 
        Parameters.channels_noise_std = zeros(1, channel);
        for ch = 1:channel
            Parameters.channels_noise_std(ch) = EstimateNoise(noise_image(:, :, ch), Parameters.patch_size); 
        end
        fprintf('The noise levels are %2.2f, %2.2f, %2.2f. \n', Parameters.channels_noise_std(1), Parameters.channels_noise_std(2), Parameters.channels_noise_std(3)); 
        
        if strcmp(denoising_type, 'NoiseImageWithReference')
            ref_image = double(imread(fullfile(['../' image_dir], ref_images(i).name))); 
        else
            ref_image = noise_image; 
        end
    end
    
    if Parameters.bEvaluate
        fprintf('The initial PSNR = %2.4f, SSIM = %2.4f \n', psnr_MultiCh(noise_image, ref_image), ssim(noise_image, ref_image)); 
    end
    Parameters.image_index = i; 
    % WNNM3 
    Parameters.channels_noise_std = repmat(sqrt(mean(Parameters.channels_noise_std.^2)), size(Parameters.channels_noise_std)); 
    [denoised_image, Parameters] = MCWNNM_Denoise(noise_image, ref_image, Parameters); 
    if Parameters.output_denoised_image
        imwrite(denoised_image/255, [denoised_image_dir ref_images(i).name]); 
        fprintf(['The denoised image is written to file ' denoised_image_dir ref_images(i).name '. \n']); 
    end
end

if Parameters.bEvaluate
    [mPSNR, idx] = max(mean(Parameters.PSNR, 2));
    mSSIM = mean(Parameters.SSIM(idx,:));
    fprintf('The best PSNR result is at %d iteration and its average PSNR = %2.4f, SSIM = %2.4f. \n', idx, mPSNR, mSSIM); 
end
fprintf('Finish! \n'); 
