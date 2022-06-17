function e = aIMTcardiac_GT(beta,t,y)

f = beta(1)*sin(2*pi*beta(2)*(t+beta(3)));

w = abs(y)/sum(abs(y));
e = sum((f-y).^2.*w.^2);

% e = sum(tanh((f-y).^2));
