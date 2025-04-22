disp('Included!')

msg = string('\n If you see this message that means this file included into PATH \n') + ...
     string(mfilename('fullpath'))
disp(msg)

clear msg

%
% Чтобы запустить - добавь к себе в пути
%
% [path2module, ~, ~] = fileparts(mfilename('fullpath'));
% cd(path2module)
% addpath(genpath( '/home/ivan/Documents/MATLAB/action_functional_modules' ))
% clear path2module