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
denoised_image_dir = ['denoised_' image_dir]; 
if ~mkdir(denoised_image_dir)
    fprintf(['Cannot make directory ' denoised_image_dir '! \n']); 
    return; 
end
diff_image_dir = ['diff_' image_dir]; 
if ~mkdir(diff_image_dir)
    fprintf(['Cannot make directory ' diff_image_dir '! \n']); 
    return; 
end

if strcmp(denoising_type, 'CleanImage')
    numImages = length(ref_images); 
else
    numImages = length(noise_images); 
end
% record all the results in each iteration 
Parameters.PSNR = zeros(numImages, 1, 'single'); 
Parameters.SSIM = zeros(numImages, 1, 'single'); 
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
    % NoiseClinic 
    if Parameters.output_noise_image
        system(['./NoiseClinic ' noise_image_dir ref_images(i).name ' ' denoised_image_dir ref_images(i).name ' ' diff_image_dir ref_images(i).name ' 2 1.0 0']); % The last one is verbose=0. 
    else
        system(['./NoiseClinic ' '../' image_dir noise_images(i).name ' ' denoised_image_dir ref_images(i).name ' ' diff_image_dir ref_images(i).name ' 2 1.0 0']); % The last one is verbose=0. 
    end
    fprintf(['The denoised image is written to file ' denoised_image_dir ref_images(i).name '. \n']); 
    if Parameters.bEvaluate
        denoised_image = double(imread(fullfile(denoised_image_dir, ref_images(i).name))); 
        Parameters.PSNR(i) = psnr_MultiCh(denoised_image, ref_image); 
        Parameters.SSIM(i) = ssim(denoised_image, ref_image); 
        fprintf('PSNR = %2.4f, SSIM = %2.4f. \n', Parameters.PSNR(i), Parameters.SSIM(i)); 
    end
end

if Parameters.bEvaluate
    fprintf('The number of images is %d and its average PSNR = %2.4f, SSIM = %2.4f. \n', numImages, mean(Parameters.PSNR), mean(Parameters.SSIM)); 
end
fprintf('Finish! \n'); 
