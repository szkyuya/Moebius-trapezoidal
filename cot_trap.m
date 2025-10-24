%{
This is a sample Matlab code of the Möbius-Transformed Trapezoidal Rule, appeared in the paper

Y. Suzuki, N. Hyvönen, and T. Karvonen. “Möbius-Transformed Trapezoidal Rule”. AMS Mathematics of Computiation (2025, published online). DOI:10.1090/mcom/40

Please cite as the following (bibtex):

@article{SHK2024, doi = {10.1090/mcom/4084}, author = {Suzuki, Yuya and Hyv\"onen, Nuutti and Karvonen, Toni}, title = {M\"obius-Transformed Trapezoidal Rule}, JOURNAL = {AMS Mathematics of Computiation}, YEAR = {2025, published online}, }

Copyright (c) 2025 Yuya Suzuki

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
%}

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




