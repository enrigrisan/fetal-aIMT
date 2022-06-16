%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Phys Med Biol
%% 2014 Nov 7;59(21):6355-71. doi: 10.1088/0022-3727/59/21/6355. Epub 2014 Oct 8.
%% Estimation of prenatal aorta intima-media thickness from ultrasound examination
%% E Veronese 1, G Tarroni, S Visentin, E Cosmi, M G Linguraru, E Grisan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('aIMT_function\')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AORTA TRACKING AND AIMT ESTIMATION IN EACH FRAME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters setting
step=5;         % Tracking step
smoothxy=0.001; % Smoothing spline regularization
passo=5;        % Spline resampling
nker=2;         % Size of morphological operator
dbf=0;          % Toggle intermediate visualization on (vessel tracking, intima-media locations)


%% Selecting and loading US video
if(~exist('xvideo'))
    [filename1,pathname1]=uigetfile('*.avi','Select avi file');
    xObj=VideoReader([pathname1,filename1]);
    fprintf('Loading frames of video %s \n', filename1)
    for ct=1:xObj.NumberofFrames,
        xtmp=read(xObj, ct);
        xim=double(xtmp(:,:,2))/255;
        xvideo(:,:,ct)=xim;
    end
end;

%% Manual initialization
%% Draw a line along the aorta
xtmp=xvideo(:,:,1);
imagesc(xtmp)
roi = drawline();
pos=roi.Position;
bc=mean(pos(:,1));
br=mean(pos(:,2));
v1=diff(pos)/norm(diff(pos));
vd=[v1(2),v1(1)];


%% Track aorta on each frame
Nframe=size(xvideo,3);
for frame =1:Nframe,

    %% Track aorta
    fprintf('Tracking abdominal aorta in frame %0.4i of %0.4i\n',frame,Nframe);
    xtmp=xvideo(:,:,frame);
    [x1,y1,diam1,dir1]=AIMTtrack(xtmp,bc,br,vd,3,200,50,step,dbf);
    [x2,y2,diam2,dir2]=AIMTtrack(xtmp,bc,br,-vd,3,200,50,step,dbf);
    xtot=[fliplr(x2),x1];
    ytot=[fliplr(y2),y1];
    
    %% Smooth centerlines and boundaries
    fprintf('\t ...Smoothing spline \n');
    l=sqrt(diff(xtot).^2+diff(ytot).^2);
    l=[0,l];
    param=cumsum(l);
    ppx = csaps(param,xtot,smoothxy);
    ppy = csaps(param,ytot,smoothxy);
    xp=fnval(ppx,[ppx.breaks(1):passo:ppx.breaks(length(ppx.breaks))]);
    yp=fnval(ppy,[ppy.breaks(1):passo:ppy.breaks(length(ppy.breaks))]);
    
    dir=atan2(yp(2:end)-yp(1:end-1),xp(2:end)-xp(1:end-1));
    dir(end+1)=dir(end);
    
    %% Refine aIMT estimate
    fprintf('\t ...Refining aIMT step 1 \n');
    mdiam(frame)=nanmedian([diam1,diam2]);
    [diamd,diamu,seld,selu]=AIMTrefine(xtmp,xp,yp,dir,mdiam(frame),nker,dbf);
    
    fprintf('\t ...Refining aIMT step 2 \n');
    nsel=find(seld==1 | selu==1);
    aIMT_vec=aIMT_thickness2(xtmp,frame,xp,yp,diamd+diamu,nsel,median(mdiam)*0.33,dbf);
    
    %% Evaluate average distal and proximal IMT
    fprintf('\t ...Evaluate proximal and distal aIMT step 1 \n');
    for ct=1:4,
        nIMT=find(aIMT_vec(:,ct)>0);
        aIMT(frame,ct)=nanmedian(aIMT_vec(nIMT,ct));
    end;
    
    ndsel=find(seld==1 & selu==1);
    if(isempty(seld))
        mdiamf(frame)=2*nanmedian(diamu(find(selu)));
    elseif(isempty(selu))
        mdiamf(frame)=2*nanmedian(diamd(find(seld)));
    else
    mdiamf(frame)=nanmedian(diamd(find(seld)))+nanmedian(diamu(find(selu)));%nanmedian(diamd(nsel)+diamu(nsel));
    end;
    %keyboard
    
    close all

    %% Save results
    fprintf('\t ...Store results \n');
    aseg(frame).x=xtot;
    aseg(frame).y=ytot;
    aseg(frame).ppx=ppx;
    aseg(frame).ppy=ppy;
    aseg(frame).xp=xp;
    aseg(frame).yp=yp;
    aseg(frame).dir=dir;
    aseg(frame).diamd=diamd;
    aseg(frame).diamu=diamu;
    aseg(frame).seld=seld;
    aseg(frame).selu=selu;
    aseg(frame).aIMT=aIMT_vec;
    
    %bc=mean(xp);
    %br=mean(yp);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimate end-systole and end-diastole values for aIMT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Ultrasound video parameters
fps=25;
f=6;
off=1;
cal=0.14;

%% Maxima and minima detection 
T=size(xvideo,3)/f;
bpm=fps/T*size(xvideo,3);

tmax=off:T:size(xvideo,3);
tmin=(off+T/2):T:size(xvideo,3);
 neigh=2;
for ct=1:length(tmax)
    nt=max(round(tmax(ct))-neigh,1):min(round(tmax(ct))+neigh,size(xvideo,3));
    nnew=find(mdiamf(nt)==max(mdiamf(nt)));
    tmax(ct)=tmax(ct)-1-neigh+nnew(1);
end;

for ct=1:length(tmin)
    nt=max(round(tmin(ct))-neigh,1):min(round(tmin(ct))+neigh,size(xvideo,3));
    nnew=find(mdiamf(nt)==min(mdiamf(nt)));
    tmin(ct)=tmin(ct)-1-neigh+nnew(1);
end;

%% Estimating pulse period
t=0:length(mdiamf)-1;
ft = @(b,t)  b(1)+b(2).*(cos(2*pi*(t-b(3))./b(4)));
fcn = @(b) sum((ft(b,t) - mdiamf).^2);
s0=[median(mdiamf),median(mdiamf(round(tmax))-mdiamf(round(tmin)))/2, tmax(1),median(diff(tmax))];
s = fminsearch(fcn, s0)

%% Estimating times of maxima and times of minima from model
if(s(2)>0),
    tmax=s(3):s(4):length(mdiamf);
    tmin=s(3)-s(4)/2:s(4):length(mdiamf)
else
    tmax=s(3)-s(4)/2:s(4):length(mdiamf);
    tmin=s(3):s(4):length(mdiamf);
end

plot(mdiamf)
hold on
plot(tmax,mdiamf(round(tmax)),'or')
plot(tmin,mdiamf(round(tmin)),'*r')
plot(ft(s,t))

sys_aIMT=[];
sys_t=[];
dia_aIMT=[];
dia_t=[];
for ct=1:length(tmax),
    n=round(tmax(ct)-5):round(tmax(ct)+5);
    [val,pos]=max(mdiamf(n));
    sys_aIMT=[sys_aIMT,val];
    sys_t=[sys_t,tmax(ct)-5+pos-1];
end
for ct=1:length(tmin),
    n=round(tmin(ct)-5):round(tmin(ct)+5);
    [val,pos]=min(mdiamf(n));
    dia_aIMT=[dia_aIMT,val];
    dia_t=[dia_t,tmin(ct)-5+pos-1];
end

%nanmedian(mdiamf(round(tmax)))*cal
%nanmedian(mdiamf(round(tmin)))*cal
%nanmedian(aIMT(round(tmax),4))*cal
%nanmedian(aIMT(round(tmin),4))*cal

