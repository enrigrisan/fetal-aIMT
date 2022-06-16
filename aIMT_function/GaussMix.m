function y=GaussMix(x,w,m,s,dbf)

% Inizializzo il Vettore dei risultati
y=zeros(size(x));

% Numero di gaussiane
N=length(w);

% Sommo le varie gaussiane 
for ct=1:N,
g=Gauss(x,m(ct),s(ct),dbf);
y=y+w(ct)*g;
end;

% Normalizzo perché sommi a 1
y=y/sum(y); 

