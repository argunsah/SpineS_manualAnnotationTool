function [bilFilt] = bil_wiener_filt(I,weinerSize,repetition)

I_temp = double(I);

for z = 1:size(I,3)
    for rep = 1:repetition      
        I_temp(:,:,z) = wiener2(I_temp(:,:,z),[weinerSize weinerSize]);
        I_temp(:,:,z) = bilateralFilter(I_temp(:,:,z));
    end
end

bilFilt = I_temp;
     