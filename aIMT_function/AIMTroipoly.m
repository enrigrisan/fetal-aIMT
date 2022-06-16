function xseg=AIMTroipoly(xim,xseg,xroi,xp,yp,diam,frac,dbf);

dx=diff(xp);
dx=[dx(1),dx];
dy=diff(yp);
dy=[dy(1),dy];
nd=sqrt(dx.^2+dy.^2);
dx=dx./nd;
dy=dy./nd;

[c,r]=find(xseg.*xroi);
y1=min(c);
y2=max(c);
x1=min(r);
x2=max(r);

d=frac*0.5*min(diam);

step1=sqrt((xp(1)-x1)^2+(yp(1)-y1)^2):-1:1;
step2=1:sqrt((xp(end)-x2)^2+(yp(end)-y2)^2);
xp1=[xp(1)-dx(1)*step1,xp,xp(end)+dx(end)*step2];
yp1=[yp(1)-dy(1)*step1,yp,yp(end)+dy(end)*step2];
n=find(xp1>x2 | xp1<x1);
xp1(n)=[];
yp1(n)=[];

prof=improfile(xim,xp1,yp1,length(xp1));
mp=nanmedian(prof)+nanstd(prof);

plab=bwlabel(imclose(prof<=mp,ones(21,1)));
for ct=1:max(plab),
    a(ct)=sum(plab==ct);
end;
[vmax,nmax]=max(a);

n=find(plab==nmax);
xp=xp1(n);
yp=yp1(n);
dx=diff(xp);
dx=[dx(1),dx];
dy=diff(yp);
dy=[dy(1),dy];
nd=sqrt(dx.^2+dy.^2);
dx=dx./nd;
dy=dy./nd;

for ct=1:length(xp),
    rowUp(ct)=xp(ct)-dy(ct)*d;
    colUp(ct)=yp(ct)+dx(ct)*d;
    rowLow(ct)=xp(ct)+dy(ct)*d;
    colLow(ct)=yp(ct)-dx(ct)*d;
end;

xs=[rowUp,fliplr(rowLow)];
ys=[colUp,fliplr(colLow)];
xseg=roipoly(zeros(size(xim)),xs,ys);