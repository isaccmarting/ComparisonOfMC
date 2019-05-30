function output_image = Patches2Image(patches, patches_weight, height, width, channel, patch_size)
% Reconstruction
output_image = zeros(height, width, channel); 
image_weight = zeros(height, width, channel); 
numRowPatches = height-patch_size+1; 
numColPatches = width-patch_size+1; 
row_index_seq = 1:numRowPatches; 
col_index_seq = 1:numColPatches; 

iter = 0; 
for ch = 1:channel
    for i = 1:patch_size
        for j = 1:patch_size
            iter = iter + 1; 
            output_image(row_index_seq-1+i, col_index_seq-1+j, ch) = output_image(row_index_seq-1+i, col_index_seq-1+j, ch) + reshape(patches(iter, :)', [numRowPatches numColPatches]); 
            image_weight(row_index_seq-1+i, col_index_seq-1+j, ch) = image_weight(row_index_seq-1+i, col_index_seq-1+j, ch) + reshape(patches_weight(iter, :)', [numRowPatches numColPatches]); 
        end
    end
end
output_image = output_image ./ image_weight; 
