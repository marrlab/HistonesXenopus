clearvars;
clc;

addpath(genpath('../HistonesXenopus')) 
H4K20_import;
% H4K20dummy_import;

count = 1;
for i = 1:size(DA(1).y,1)
   for j =  1:size(DA(1).y,2)
       clearvars a b
       for irep = 1:length(DA)
           a(irep) = exp(DA(irep).y(j,i));
           b(irep) = exp(DB(irep).y(j,i));
       end
       [h,p] = ttest2(a,b);
       H(count) = h;
       P(count) = p;
       count = count+1;
   end
end

length(find(P<0.05))