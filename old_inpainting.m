close all;
clear; clc;

dataset_path = 'cnn-results/DFNet_rectangle';
output_path = 'output';
masks_path = 'masks';
output_flist_name = 'flist';
masks_flist_name = 'flist';
result_path = 'classic_inpainting_2';
extension = '.jpg';

image_height = 256;
image_width = 256;
dataset_size = 1228;

N = 92;
sigma = floor(N/4);

disp('Candidate names reading...');
candidates = read_flist([dataset_path, '/', output_path, '/', output_flist_name], dataset_size);
disp('   success');

disp('Candidate names reading...');
masks = read_flist([dataset_path, '/', masks_path, '/', masks_flist_name], dataset_size);
disp('   success');

for index=1:dataset_size
    
    candidate_full_path = [dataset_path, '/', output_path, '/', char(candidates(index))];
    mask_full_path = [dataset_path, '/', masks_path, '/', char(masks(index))];
    try
        disp('Inpaint ' + candidates(index) + ' ...');
        
        candidate = imread(candidate_full_path);
        result = candidate;
        
        %// Magic for the best border detection
        mask = imread(mask_full_path);
        gray_mask = rgb2gray(mask);
        cont_mask = double(gray_mask);
        cont_mask = cont_mask / max(cont_mask(:));
        cont_mask(cont_mask < 0.8) = 0;
        cont_mask(cont_mask >= 0.8) = 1;
        
        inpaintedR = uint8(inpaint(double(candidate(:,:,1)), double(cont_mask), N, sigma));
        inpaintedG = uint8(inpaint(double(candidate(:,:,2)), double(cont_mask), N, sigma));
        inpaintedB = uint8(inpaint(double(candidate(:,:,3)), double(cont_mask), N, sigma));
        
        res = cat(3, inpaintedR, inpaintedG, inpaintedB);
        full_result_path = [result_path, '/', char(candidates(index))];
        imwrite(res, full_result_path);
    catch ME
        disp([' ERROR! name: ', candidate_full_path, ' error: ', char(ME.identifier)]);
    end
end

function [inpainted] = inpaint(aImage, aMask, aMaxHoleSize, aFilterSigma)
[imageW, imageH] = size(aImage);
mult = -0.5/ (aFilterSigma * aFilterSigma);
inpainted = zeros(imageW,imageH);

for i=1:imageW
    for j=1:imageH
        if (aMask(i,j) == 1)
            inpainted(i,j) = aImage(i,j);
            continue;
        end
        
        
        sum = 0.0;
        sumFilt = 0.0;
        
        for v=-aMaxHoleSize:aMaxHoleSize
            for h=-aMaxHoleSize:aMaxHoleSize
                
                if (aMask(i,j) == 0)
                    valFilt = exp(mult * (v * v + h * h));
                    
                    indexI = i + v;
                    if(indexI < 1)
                       indexI = 1; 
                    end
                    if(indexI > 256)
                       indexI = 256;
                    end
                    
                    indexJ = j + h;
                    if(indexJ < 1)
                       indexJ = 1;
                    end
                    if(indexJ > 256)
                       indexJ = 256;
                    end
                    
                    sumFilt = sumFilt + valFilt;
                    sum = sum + valFilt * aImage(indexI, indexJ);
                end
                
            end
        end
        
        inpainted(i,j) = sum / sumFilt;
    end
end
end

function [names] = read_flist(path, size)
names = strings(size,1);
fid = fopen(path);
tline = fgetl(fid);
counter = 1;
while ischar(tline)
    %disp(tline)
    names(counter) = tline;
    counter = counter + 1;
    tline = fgetl(fid);
end
fclose(fid);
end