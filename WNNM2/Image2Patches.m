function patches = Image2Patches(target_image, patch_size, channel, numPatches)
% record the non-local patch set and the index of each patch in
% of seed patches in image
patches = zeros(channel*patch_size*patch_size, numPatches, 'double'); 
iter = 0; 
for ch = 1:channel
    for i = 1:patch_size
        for j = 1:patch_size
            iter = iter + 1; 
            image_patch  = target_image(i:end-patch_size+i, j:end-patch_size+j, ch); 
            patches(iter, :) = image_patch(:)'; 
        end
    end
end
