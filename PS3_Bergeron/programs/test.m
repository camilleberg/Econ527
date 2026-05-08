% Example matrix and threshold
A = [10, 20, 5, 2; 5, 15, 25, 30; 8, 12, 18, 22];
threshold = 12;

% Returns row and column vectors
[B, C] = find(A < threshold);
disp('Row and column indices of elements less than threshold:');
disp(B);
disp(C);

newValue = 0; % Value to replace elements below threshold

A(A < threshold) = newValue;
disp('Updated matrix A:');
disp(A);

