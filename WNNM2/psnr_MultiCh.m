function res = psnr_MultiCh(target_image, ref_image)

image_error = target_image - ref_image; 
[height, width, channel] = size(target_image); 
mse = 0; 
for ch = 1:channel
    mse = mse + mean(mean(image_error(:, :, ch).^2)); 
end
mse = mse / channel; 
res = 10 * log10(255^2/mse); 
