function [lim,peak]=aIMTfindpeak(prof,tol,imax,Nmix,dbf)

lim=[];
peak=[];

N=length(prof);

int=mean(prof(1:fix(N/2)));
ext=mean(prof(end-fix(N/2)+1:end));
if int<tol*ext %punto superiore affidabile
    peaksn=fpeak([1:length(prof)],prof,10,[1,N,0,imax]);
    if (~isempty(peaksn) && peaksn(1,1)==1)
        peaksn=peaksn(2:end,:);
    end
    
    nlim=find(prof<mean(prof));
    peaksn(peaksn(:,1)<nlim(end),:)=[];
    
    if (~isempty(peaksn) && any(peaksn(:,2)<0.05)) % tolgo i minimi sotto una certa soglia
        to_del=find(peaksn(:,2)<0.05);
        for k=1:length(to_del)
            peaksn(to_del(k),:)=[0 0];
        end
        peaksntmp1=nonzeros(peaksn(:,1));
        peaksntmp2=nonzeros(peaksn(:,2));
        clear peaksn
        peaksn=[peaksntmp1';peaksntmp2']';
        clear peaksntmp1 peaksntmp2
    end
    
    peak=peaksn(1,1);
end;