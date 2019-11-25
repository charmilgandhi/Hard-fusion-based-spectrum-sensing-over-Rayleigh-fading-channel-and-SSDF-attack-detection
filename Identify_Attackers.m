%Identify Honest and False Secondary Users (SUVs)

clc
clear all;

h = 30;
D = randi([0 1],h,h);

trust_value = zeros(h);
for k=1:h
honesty = randi([2000,4000], 1,1);
false = 4000 - honesty;
trust_value(k) = (1 + honesty)/(2 + (honesty + false));
end

x = zeros(h,h);
for ii = 1:h
    for jj = 1:h
        xor_matrix(:,jj) = xor(D(:,ii),D(:,jj));
    end
    temp_sum = zeros(1,h);
    for k = 1:h
        temp_sum = temp_sum + (k*xor_matrix(k,:));
    end
    x(ii,:) = temp_sum;
end

XD_normalize = zeros(h,h);

for i=1:h
     max1 = max(x(i,:));
     for j=1:h
        XD_normalize(i,j) = x(i,j)/max1; 
     end
end

threshold = 0.8;
for i=1:h
     for j=1:h
       if(trust_value(j) > threshold)
           omega(j) = XD_normalize(i,j);
       else
            omega(j) = NaN; 
       end
     end
     [M,L] = min(omega);
     Y1(i) = L; 
     [M,I] = max(omega);
     Y2(i) = I;
 end
 Y1
 Y2