function displayProgress(msg,iter,total)
% display message for updating progress

% setup for first iteration
if(iter == 1)
    fprintf(1,[msg ' ']); % an extra space takes care of the first backspace
end

% update through the iterations
fprintf(1,[repmat('\b',1,numel(num2str(iter-1))) '%d'],iter);

% a new line for the end of the loop
if(iter == total)
    fprintf('\n');
end