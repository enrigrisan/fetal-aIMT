function e = aIMTcardiac_GT2(beta,t,y)

% f1 = sin(2*pi*beta(1)*(t+beta(2)));
f1 = cos(2*pi*beta(1)*t - 2*pi*beta(2));
w = ones(size(y));%abs(y)/sum(abs(y));
beta1 = lscov(f1',y',sqrt(w'));

f = beta1 * f1;
e = sum((f-y).^2.*w.^2);

% e = sum(tanh((f-y).^2));
