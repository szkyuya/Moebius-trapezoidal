%% Mobius transform
phi= @(z) -cot(z./2); % mobius transform
phi_=@(z) 1./(2.*(sin(z./2)).^2); % Jacobian
M=16; %number of points; calculated up until 2^M points
err=zeros(2,M);
%% Define your weight
rho= @(x) exp(-(x.^2)./2)./sqrt(2*pi);
%% Define your integrand function
p=1;%parameter in the integrand, |x|^p
f= @(x) abs(x).^p;

trueI=sqrt((2^p)/pi)*gamma((p+1)/2); %true integral value, only for plotting the absolute error  

for ni=2:M
n=2^ni+1; % number of points
err(1,ni)=n;
acc=0;
for j=1:n
theta=2*pi*(j)/n+1/(2*n);% equidistant point with a shift
x=phi(theta);
acc=acc+f(x)*rho(x)*phi_(theta)*2*pi/(n);

end
disp(['current estimate: ',num2str(acc,10)])
err(2,ni)=abs(trueI-acc); 
end

figure;
loglog(err(1,:),abs(err(2,:)))
title('Convergence graph') 
xlabel('Number of points') 
ylabel('Absolute error') 

rate = polyfit(log(err(1,4:end)),log(abs(err(2,4:end))),1)


