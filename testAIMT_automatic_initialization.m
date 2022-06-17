%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Phys Med Biol
%% 2014 Nov 7;59(21):6355-71. doi: 10.1088/0022-3727/59/21/6355. Epub 2014 Oct 8.
%% Estimation of prenatal aorta intima-media thickness from ultrasound examination
%% E Veronese 1, G Tarroni, S Visentin, E Cosmi, M G Linguraru, E Grisan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('aIMT_function\')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AORTA TRACKING AND AIMT ESTIMATION IN EACH FRAME WITH AUTOMATIC INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Tracking parameters setting
step=5;         % Tracking step
smoothxy=0.001; % Smoothing spline regularization
passo=5;        % Spline resampling
nker=2;         % Size of morphological operator
maxd=100;       % maximum size of diameter
maxfail=5;      % maximum number of consecutive fails before stopping tracking

Th_grad=0.5;
Trel=0.3;

%% Multiscale aorta detection parameters
pars.dim=100;
pars.sigma=5:5:20;
pars.theta=-24:8:24;
pars.resize=0.5;
smooth = 0.02;
struelem=ones(21);

%% Thickness refinement and acceptance parameters
parthick.N= 11;
parthick.NN=15;
parthick.tol=0.2;     % Internal-external intensity tolerance
parthick.tol_aIMT=10; % maximum thickness in pixels

%% INTERMEDIATE RESULTS FLAG
%% If dbf=1 the functions will provide intermediate results, figures and information
dbf=0;

%% OPTIONAL: DEFINITION OF THE ANALYSIS INTERVAL FOR VIDEOS OF EACH PATIENT
%% Manually insert the first frame in which the aorta can be seen and analyzed
%% The first value will be used for the first file, the second for the second and so on
framestart=-1;
framefinish=-1;

%% LOAD VIDEO
[filename1,pathname1]=uigetfile('*.avi','Select avi file');
xObj=VideoReader([pathname1,filename1]);
fprintf('Loading frames of video %s \n', filename1)
for ct=1:xObj.NumberofFrames,
    xtmp=read(xObj, ct);
    xim=double(xtmp(:,:,2))/255;
    seq(:,:,ct)=xim;
end

% Initialize the output structure
aseg=[];
mdiamf=[];
mdiam=[];
aIMT=[];

% Ask the user to delineate a region of interest (ROI) around the aorta
% on an image obtained as the mean of the video
frame=1;
if(framefinish==-1)
    framefinish=size(seq,3);
end
if(framestart==-1)
    framestart=size(seq,3);
end
xtmp=mean(seq(:,:,framestart:framefinish),3);
imagesc(xtmp);
xroi=roipoly();


% Define the cropping coordinates around the valid region of the ROI
[xcrop,ycrop]=find(xroi);
crop=[min(xcrop),max(xcrop),min(ycrop),max(ycrop)];

if((crop(2)-crop(1))<200)
    pars.resize=1;
    maxd=100;
else
    pars.resize=0.5;
    maxd=100;
end


maxd=input('maxd:');

%% Analyze each frame of the video
for frame=framestart(ctfile):framefinish(ctfile)

    % Extract the image xtmp
    xtmp=seq(:,:,frame);
    dims = (-80:10:80);
    rots = 0;
    th_area = 800;
    th_ecc = [0.92,1];

    % Find a starting position for tracking the aorta
    % and evaluate if the aorta is visible within the image under analysis
    [aorta_coor, xaorta]=AIMTinitialize(xtmp, xroi,pars, dbf);
    cc = bwconncomp(xaorta);
    stats = regionprops(cc, 'Area','Eccentricity','Centroid');
    areas = [stats.Area];
    eccs = [stats.Eccentricity];
    valid_areas = areas>=th_area;
    valid_eccs = (eccs>th_ecc(1))&(eccs<th_ecc(2));

    % If the detected aorta has characteristics not compatible with a
    % vessel set the lumen_id_flag to 0
    % If the lumen_id_flag=0 the image will not be analayzed
    if sum(valid_areas.*valid_eccs)>0
        lumen_id_flag = 1;
    else
        lumen_id_flag = 0;
    end

    if lumen_id_flag==1

        % Compute the coordinate [br,bc] of the center of the approximate aorta
        % Compute the direction V
        [r,c]=find(xaorta);
        Vc=cov(r,c);
        [V,D]=eig(Vc);
        if(D(1,1)>D(2,2))
            vd=V(:,1);
        else
            vd=V(:,2);
        end;
        ns=find(abs(c-mean(c))<3);
        bc=mean(c);
        br=mean(r(ns));

        % Track the lumen center and diamater of tha aorta starting from [br,bc] along the
        % positive direction
        [x1,y1,diam1,dir1]=AIMTtrack(xtmp.*xroi,xroi,bc,br,V(:,2),maxfail,200,maxd,step,dbf);

        % Track the lumen of tha aorta starting from [br,bc] along the
        % negative direction
        [x2,y2,diam2,dir2]=AIMTtrack(xtmp.*xroi,xroi,bc,br,-V(:,2),maxfail,200,maxd,step,dbf);

        % Check the position of the estimated center point sample
        n1=find(x1>crop(3) & x1<crop(4));
        n2=find(x2>crop(3) & x2<crop(4));

        % Merge the positive and negative tracking samples
        xtot=[fliplr(x2(n2)),x1(n1)];
        ytot=[fliplr(y2(n2)),y1(n1)];

        % If the number of samples is tto small something has gone
        % wrong
        if(length(xtot)<10)
            keyboard
        end

        % Interpolate the samples along the lumen center with a cubic
        % smoothing spline for a denser and smoother representation
        l=sqrt(diff(xtot).^2+diff(ytot).^2);
        l=[0,l];
        param=cumsum(l);
        ppx = csaps(param,xtot,smoothxy);
        ppy = csaps(param,ytot,smoothxy);

        xp=fnval(ppx,[ppx.breaks(1):passo:ppx.breaks(length(ppx.breaks))]);
        yp=fnval(ppy,[ppy.breaks(1):passo:ppy.breaks(length(ppy.breaks))]);

        dir=atan2(yp(2:end)-yp(1:end-1),xp(2:end)-xp(1:end-1));
        dir(end+1)=dir(end);

        % use a running average as initial estimate of the diameter
        mdiam(frame)=nanmedian([diam1(diam1>0),diam2(diam2>0)]);
        if(frame>1)
            mind=max(max(mdiamf)-mean(mdiamf),0.5*max(mdiamf));
            if( isnan(mdiam(frame)) | mdiam(frame)<mind)
                mdiam(frame)=mdiamf(end);
            end;
        end;

        % Refine the estimation of centerpoints and diamaters.
        % evaluating the points for which a reliable estimate of the
        % vessel widh might be obtained
        [diamd,diamu,seld,selu]=AIMTrefine2(xtmp,xp,yp,dir,mdiam(frame),nker,Th_grad,dbf);


        % Find the selected points and estimate the aIMT only at their
        % position
        nsel=find(seld==1 | selu==1);
        nsel(nsel<0.1*length(xp))=[];
        nsel(nsel>0.9*length(xp))=[];
        aIMT_vec=aIMT_thickness3(xtmp,frame,xp,yp,diamd+diamu,nsel,max(0.3*mean(diamd+diamu),25),Trel,dbf);

        % Evaluate median values of aIMT
        for ct=1:4,
            nIMT=find(aIMT_vec(:,ct)>0);
            aIMT(frame,ct)=nanmedian(aIMT_vec(nIMT,ct));
        end;

        % Evaluate median value of the lumen diameter
        ndsel=find(seld==1 & selu==1);
        if(isempty(seld))
            mdiamf(frame)=2*nanmedian(diamu(find(selu)));
        elseif(isempty(selu))
            mdiamf(frame)=2*nanmedian(diamd(find(seld)));
        else
            mdiamf(frame)=nanmedian(diamd(find(seld)))+nanmedian(diamu(find(selu)));%nanmedian(diamd(nsel)+diamu(nsel));
        end;
        %keyboard

        % Save all results in the output structure
        aseg(frame).bc=bc;
        aseg(frame).br=br;
        aseg(frame).V=V;
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

    else
        aIMT_vec=-1;
    end
    fprintf('Frame %i, aIMT %f %f %f %f\n', frame, nanmedian(aIMT_vec))
end
%meanseg=aseg;
save(fullfile('..\Risultati',[filename(1:end-4),'.mat']),'aseg');%,'-append')

