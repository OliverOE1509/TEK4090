% Generate a vector of 500 normally distributed random numbers with mean =
% 10 and standard deviation =5
foo = 10 + 5*randn(500, 1);
% Use a normplot to see whether the sample appears to be normally
% distributed
normplot(foo)

return
figure
% Generate a vector of 500 random numbers drawn from a gamma distribution
% with aplha = 5 and beta = 3
bar = gamrnd(5,3, 500,1)
% Use a normplot to see whether the sample appears to be normally
% distributed
normplot(bar)