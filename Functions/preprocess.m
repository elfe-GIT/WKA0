function sys = preprocess(fileName,sheetName)
%UNTITLED3 Summary of this function goes here
%   read parameters and compose system matrix A and right-hand-side b
%   return matrices as Matlab-container (hashmap)

p = readtable(fileName,'Sheet',sheetName,'Range','A1:D29','ReadRowNames',true);
data = p(:,{'number','unit'});

% convert "data" to class of type planets
sys = model(data);

end

