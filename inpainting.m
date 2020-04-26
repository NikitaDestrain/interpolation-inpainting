close all;
clear; clc;

dataset_path = 'places_for3';
output_path = 'output';
masks_path = 'masks';
output_flist_name = 'flist_hole';
masks_flist_name = 'flist_mask';
result_path = 'classic_inpainting';
extension = '.jpg';

image_height = 256;
image_width = 256;
dataset_size = 1002;

N = 91;
sigma = 20;

shift = uint8(floor(N / 2));
quad_N = N * N;

disp('Gaussian kernel generating...')
h_special = fspecial('gaussian', [N N], sigma);
disp('   success');

disp('Candidate names reading...');
candidates = read_flist([dataset_path, '/', output_path, '/', output_flist_name], dataset_size);
disp('   success');

disp('Candidate names reading...');
masks = read_flist([dataset_path, '/', masks_path, '/', masks_flist_name], dataset_size);
disp('   success');

%candidates(1) = candidate_full_path;
%masks(1) = mask_full_path;


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
        
        %// Try to find mask coord
        start_x = -1;
        start_y = -1;
        end_x = -1;
        end_y = -1;
        
        for i=1:image_height
            for j=1:image_width
                pixel = cont_mask(i,j);
                if (pixel == 0)
                    if (start_x == -1)
                        start_x = j;
                    end
                    if (start_y == -1)
                        start_y = i;
                    end
                    
                    end_x = j;
                    end_y = i;
                end
            end
        end
        
        hole_size_x = 1 + end_x - start_x;
        hole_size_y = 1 + end_y - start_y;
        
        %// Cut part for conv
        %// Warning swap x and y cuz we have fault when we calculate size
        for runner_x=0:hole_size_x
            for runner_y=0:hole_size_y
                tmp_x = start_x + runner_x;
                tmp_y = start_y + runner_y;
                
                left_y = tmp_y - shift;
                %// The best padding ever for optimize time
                if (left_y < 1)
                    left_y = 1;
                end
                
                rigth_y = tmp_y + shift;
                %// The best padding ever for optimize time
                if (rigth_y > image_height)
                    rigth_y = image_height;
                end
                
                left_x = tmp_x - shift;
                %// The best padding ever for optimize time
                if (left_x < 1)
                    left_x = 1;
                end
                
                rigth_x = tmp_x + shift;
                %// The best padding ever for optimize time
                if (rigth_x > image_width)
                    rigth_x = image_width;
                end
                
                %// Thanks for 3 color channel
                cc_r = candidate(left_y:rigth_y, left_x:rigth_x, 1);
                cc_g = candidate(left_y:rigth_y, left_x:rigth_x, 2);
                cc_b = candidate(left_y:rigth_y, left_x:rigth_x, 3);
                
                cc = cat(3, cc_r, cc_g, cc_b);
                
                %// Predict
                tmp = imfilter(cc, h_special);
                
                delimeter = quad_N - (runner_x * runner_y);
                
                norm_cc_r = sum(sum(tmp(:,:,1))) / delimeter;
                norm_cc_g = sum(sum(tmp(:,:,2))) / delimeter;
                norm_cc_b = sum(sum(tmp(:,:,3))) / delimeter;
                
                result(tmp_y, tmp_x, 1) = norm_cc_r;
                result(tmp_y, tmp_x, 2) = norm_cc_g;
                result(tmp_y, tmp_x, 3) = norm_cc_b;
            end
        end
        
        %figure(); imshow(result);
        
        full_result_path = [result_path, '/', char(candidates(index)), '_inpaint', extension];
        disp(['   save result to ', full_result_path]);
        imwrite(result, full_result_path);
    catch ME
        disp([' ERROR! name: ', candidate_full_path, ' error: ', char(ME.identifier)]);
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