function [Parameters, ref_images, noise_images] = Initialize(denoising_type, image_dir)
% set parameters 

% parameters for denoising
Parameters.window_size = 20;                                   % Non-local patch searching window
Parameters.reBlkMacLoop = 2;                            % InnerLoop Num of between re-blockmatching
Parameters.patch_size = 6;                                 % Patch size, larger values would get better performance, but will be slower
Parameters.step = 5;                             % The step between neighbor patches, smaller values would get better performance, but will be slower
Parameters.ori_numSimilarPatches = 70; 
% WNNM3 
Parameters.K2 = 10;                            % total iter numbers
Parameters.delta = 0.1;                 % iterative regularization parameter
Parameters.lambda = 0.6;                   % for different noise levels, this parameter should be tuned to achieve better performance

% parameters for ADMM
Parameters.K1 = 10;
Parameters.rho = 3; 
Parameters.mu = 1.001; 

Parameters.output_denoised_image = true; 

if strcmp(denoising_type, 'CleanImage')
    Parameters.channels_noise_std = [40 20 30]; 
    
    Parameters.bEvaluate = true; 
    Parameters.output_noise_image = true; 
    ref_images = dir(fullfile(['../' image_dir], '*.png')); 
    noise_images = ref_images; 
elseif strcmp(denoising_type, 'NoiseImageWithReference')
    Parameters.K2 = 2; 
    Parameters.delta = 0; 
    Parameters.lambda = 4; 
    Parameters.rho = 6; 
    Parameters.mu = 1; 

    Parameters.bEvaluate = true; 
    Parameters.output_noise_image = false; 
    ref_images = dir(fullfile(['../' image_dir], '*mean.png')); 
    noise_images = dir(fullfile(['../' image_dir], '*real.png')); 
else % if strcmp(denoising_type, 'NoiseImageWithoutReference')
    Parameters.K2 = 2; 
    Parameters.delta = 0; 
    Parameters.lambda = 1.5; 
    Parameters.rho = 6; 
    Parameters.mu = 1; 

    Parameters.bEvaluate = false; 
    Parameters.output_noise_image = false; 
    noise_images = dir(fullfile(['../' image_dir], '*.png')); 
    ref_images = noise_images; 
end