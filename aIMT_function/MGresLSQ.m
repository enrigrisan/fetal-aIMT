function res=MGresLSQ(par,x,data,dbf)

% Calcolo il numero di gaussiane
N=length(par)/3;

% Estraggo i parametri w,m,s separati
w=par(1:N);
m=par(N+1:2*N);
s=par(2*N+1:3*N);

% Calcolo della mistura
MG=GaussMix(x,w,m,s,dbf);

% Calcolo della somma dei residui quadrati
res=(MG-data)';
