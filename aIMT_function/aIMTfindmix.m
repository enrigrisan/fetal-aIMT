function [lim,peak]=aIMTfindmix(prof,Nmix,dbf)


if(isempty(Nmix))
    Nmix=3;
end;
step=1;

lim=[];
peak=[];
try
    if(sum(prof)>0)
        N=length(prof);
        m=sum(prof'.*[1:N])/sum(prof);
        s=sqrt(sum(prof'.*([1:N]-m).^2)/sum(prof));
        
        mixspace=[-(Nmix-1)/2:-1,1:(Nmix-1)/2];
        m0=[m,m+mixspace*20];
        s0=[s,ones(1,Nmix-1)*s/10];
        w0=[0.4,ones(1,Nmix-1)*0.6/(Nmix-1)];
        
        lb=zeros(1,3*Nmix);
        ub=[ones(1,Nmix),N*ones(1,Nmix),N/2*ones(1,Nmix)];
        par=[w0,m0,s0];
        
        options=optimset('Display','Off','Diagnostics','Off');
        %try
            paropt=lsqnonlin(@MGresLSQ,par,lb,ub,options,[1:N]',prof/sum(prof),0);
        %catch
        %    keyboard;
        %end
        
        mo=paropt(Nmix+1:2*Nmix);
        so=paropt(2*Nmix+1:3*Nmix);
        wo=paropt(1:Nmix);
        
        xp=[1:step:N];
        yg=GaussMix(xp,wo,mo,so,0);
        
        [sorts,sind]=sort(so,'ascend');
        
        if(mo(sind(1))<mo(sind(2)))
            m1=mo(sind(1));
            m2=mo(sind(2));
            s1=so(sind(1));
            s2=so(sind(2));
            w1=wo(sind(1));
            w2=wo(sind(2));
        else
            m1=mo(sind(2));
            m2=mo(sind(1));
            s1=so(sind(2));
            s2=so(sind(1));
            w1=wo(sind(2));
            w2=wo(sind(1));
        end;
        
        g1=w1*Gauss(xp,m1,s1,dbf);
        g2=w2*Gauss(xp,m2,s2,dbf);
        
        if abs(m1-m2)>2*sorts(1) & abs(m1-m2)>3 & m1>N/2 & m2>N/2
            nvalley=find(xp>m1 & xp<m2);
            
            nmin1=find(prof(nvalley)==min(prof(nvalley)));
            nmin2=find(yg(nvalley)==min(yg(nvalley)));
            nmin=nanmin(nmin1,nmin2);
            %nmin=find(g1(nvalley)>g2(nvalley));
            try
                peak=nvalley(nmin(end));
            catch
                keyboard
            end;
            
            nlumen=find(xp<m1 & yg<0.75*yg(fix(m1)));
            %ydist=(yg(nlumen)-yg(peak)).^2;
            %nmin=find(ydist==min(ydist));
            
            ec1=(peak-nmin(end))<sorts(1)*0.5;
            ec2=(yg(round(m1))-yg(peak))<0.05*max(yg);
            ec3=(yg(round(m2))-yg(peak))<0.05*max(yg);
            ec4=isempty(nlumen);
            
            ec1=0;
            ec3=0;
            if(ec1 | ec2 | ec3 | ec4)
                peak=[];
                lim=[];
            else
                lim=xp(nlumen(end));
            end;
            
            if(dbf)
                figure;
                plot(prof)
                hold on
                plot(yg*sum(prof),'r')
                if(~isempty(peak))
                    plot([peak,peak],[0,max(prof)],'m')
                    plot([lim,lim],[0,max(prof)],'m')
                end
            end
        elseif(abs(m1-m2)>2*sorts(1) & abs(m1-m2)>3 & m2>N/2)
            
            try
                if(m1>1 && yg(round(m1))>yg(round(m2)))
                    ms=m1;
                else
                    ms=m2;
                end;
            catch
                keyboard
            end;
            
            [mw,nw]=max(yg);
            ms=nw;
            nlumen=find(xp<ms & yg<0.75*yg(fix(ms)));
            %lim=xp(nlumen(end));
            
            nvalley=find(xp>ms & yg<0.75*yg(fix(ms)));
            if(isempty(nvalley))
                nvalley=xp(end);
            end;
            
            ec4=isempty(nlumen);
            if(ec4)
                peak=[];
                lim=[];
            else
                try
                    lim=xp(nlumen(end));
                    peak=nvalley(1);
                catch
                    keyboard;
                end
            end;
            
            if(dbf)
                figure;
                plot(prof)
                hold on
                plot(yg*sum(prof),'r')
                if(~isempty(peak))
                    plot([peak,peak],[0,max(prof)],'m')
                    plot([lim,lim],[0,max(prof)],'m')
                end
            end
        end;
    end;
end



