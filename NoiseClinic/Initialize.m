function [Parameters, ref_images, noise_images] = Initialize(denoising_type, image_dir)
% set parameters 

Parameters.patch_size = 6;                                 % Patch size, larger values would get better performance, but will be slower
%arameters.output_denoised_image = true; 

if strcmp(denoising_type, 'CleanImage')
    Parameters.channels_noise_std = [40 20 30]; 
    
    Parameters.bEvaluate = true; 
    Parameters.output_noise_image = true; 
    ref_images = dir(fullfile(['../' image_dir], '*.png')); 
    noise_images = ref_images; 
elseif strcmp(denoising_type, 'NoiseImageWithReference')
    Parameters.bEvaluate = true; 
    Parameters.output_noise_image = false; 
    ref_images = dir(fullfile(['../' image_dir], '*mean.png')); 
    noise_images = dir(fullfile(['../' image_dir], '*real.png')); 
else % if strcmp(denoising_type, 'NoiseImageWithoutReference')
    Parameters.bEvaluate = false; 
    Parameters.output_noise_image = false; 
    noise_images = dir(fullfile(['../' image_dir], '*.png')); 
    ref_images = noise_images; 
end
