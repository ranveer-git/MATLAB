clc
close all
clear all
D = [1.24 1.5 1.7 1.9 2.1 2.25 2.3 2.48 2.5 2.67]

D1 = 1.2:0.2:2.6

transformationMatrix = zeros(2,length(D));
for dIndex=1:length(D);
    dCurrent = D(dIndex);
    for bucketIndex = 1 : length(D1);
       currentBucket = D1(bucketIndex);
       if(currentBucket > dCurrent)
        break
       end
       if(currentBucket <= dCurrent)
          transformationMatrix(1,dIndex) = dIndex;
          transformationMatrix(2,dIndex) = bucketIndex;
       end
    end
end

for i = 1:
    

