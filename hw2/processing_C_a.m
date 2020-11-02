load('data_backup.mat')
max_error = zeros(size(Ms));
for k = 1:length(Ms)
    max_error(k) = max(abs(errors(k,:)));
end
plot(log(hs),log(max_error),'b')
coeffs=polyfit(log(hs),log(max_error),1)