close all;
clear; clc;

dataset_path = 'classic_gauss';
output_path = 'output';
masks_path = 'masks';
output_flist_name = 'flist_hole';
masks_flist_name = 'flist_mask';
result_path = 'result';
extension = '.jpg';

image_height = 256;
image_width = 256;
dataset_size = 1228;
sigma_max = 10;

N = 92;
R = 15;

disp('Candidate names reading...');
candidates = read_flist([dataset_path, '/', output_path, '/', output_flist_name], dataset_size);
disp('   success');

disp('Candidate names reading...');
masks = read_flist([dataset_path, '/', masks_path, '/', masks_flist_name], dataset_size);
disp('   success');

for sigma=1:sigma_max
    for index=1:dataset_size
        
        candidate_full_path = [dataset_path, '/', output_path, '/', char(candidates(index)), extension];
        mask_full_path = [dataset_path, '/', masks_path, '/', char(masks(index)), extension];
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
            
            [startI, endI, startJ, endJ] = get_mask_coord(cont_mask);
            
            inpaintedR = inpaint(double(candidate(:,:,1)), double(cont_mask), N, sigma, R, startI, endI, startJ, endJ);
            inpaintedG = inpaint(double(candidate(:,:,2)), double(cont_mask), N, sigma, R, startI, endI, startJ, endJ);
            inpaintedB = inpaint(double(candidate(:,:,3)), double(cont_mask), N, sigma, R, startI, endI, startJ, endJ);
           
            res = cat(3, inpaintedR, inpaintedG, inpaintedB);
            res = uint8(res);
            
            full_result_path = [result_path, '/', char(candidates(index)), extension];
            imwrite(res, full_result_path);
        catch ME
            disp([' ERROR! name: ', candidate_full_path, ' error: ', char(ME.identifier)]);
        end
    end
end

function [startI, endI, startJ, endJ] = get_mask_coord(aMask)
[imageW, imageH] = size(aMask);

startI = -1;
startJ = -1;

for i=1:imageW
    for j=1:imageH
        if (aMask(i,j) == 0)
            
            if(startI == -1)
                startI = i;
            end
            
            if(startJ == -1)
                startJ = j;
            end
            
            endI = i;
            endJ = j;
            
        end
    end
end
end

function inpainted = inpaint(aImage, aMask, aMaxHoleSize, aFilterSigma, R, startI, endI, startJ, endJ)
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
            
            indexI = i + v;
            for h=-aMaxHoleSize:aMaxHoleSize
                indexJ = j + h;
                if ((indexI < 1) || (indexI > imageH) || (indexJ < 1)|| (indexJ > imageW))
                    continue;
                end
                
                if (aMask(indexI,indexJ) == 1 &&...
                        indexI > (startI - R) && indexI < (endI + R) &&...
                        indexJ > (startJ - R) && indexJ < (endJ + R))
                   
                    valFilt = exp(mult * (v * v + h * h));
 
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