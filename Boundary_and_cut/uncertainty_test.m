mu = Pd(2:end);
sigma = 0.1*mu;
data = [];
for i = 1:32
    data(i,:) = normrnd(mu(i),sigma(i),[1,5000]);
end
C = cov(data');
upper = chol(C,'upper');
