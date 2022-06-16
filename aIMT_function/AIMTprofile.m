function [diam,sel]=AIMTprofile(pf,thmax,pstart,mdiam,gradmin, graymin, dbf)

[vmax,nmax]=max(pf);
[vmin,nmin]=min(pf(1:nmax));

dv=vmax-vmin;
mval=mdiam-pstart;
thdgray=0.15;
%th=thmax*dv;

[center, U, obj_fcn] = fcm_my(pf', 2);
%th=mean(center);
th=min(center)+0.5*abs(diff(center));

%fprintf('%f\t%f\t%f\n',thmax,min(center),th);
bl=bwlabel(pf(1:nmax)<th);
nrise=find(bl==1,1,'last');

bl=bwlabel(pf(1:nmax)>th);
nsel=find(bl==1);
vmax=max(pf(nsel));
nmax=mean(nsel(find(pf(nsel)==vmax)));

if(isempty(nrise))
    nrise=mval;
end;

cvd=abs(nrise-mval)/mval;
if(cvd>0.2)
    nrise=mval;
end;

dp=nmax-nrise;
dv=vmax-vmin;

nint=fix(nrise-mdiam*0.4:nrise);
next=fix(nrise:nrise+mdiam*0.4);
dgray=abs(mean(pf(nint))-mean(pf(next)));

if(dgray<thdgray | dv/dp<gradmin | dv<graymin)
    sel=-1;
    diam=mval;
else
    sel=1;
    diam=nrise;
end;