%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Phys Med Biol
%% 2014 Nov 7;59(21):6355-71. doi: 10.1088/0022-3727/59/21/6355. Epub 2014 Oct 8.
%% Estimation of prenatal aorta intima-media thickness from ultrasound examination
%% E Veronese 1, G Tarroni, S Visentin, E Cosmi, M G Linguraru, E Grisan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('aIMT_function\')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ESTIMATION OF CARDIAC CYCLE AND STIFFNESS FROM TRACKED AORTA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run testAIMT on a video to track aorta and save the corresponding structure "aseg"
if(not(exist('aseg')))
    [filename1,pathname1]=uigetfile('*.mat','Select mat file containing the aseg structure');
    load(fullfile(pathname1, filename1));
end

for ct=1:length(aseg)
    if(or(isempty(aseg(ct).diamd),isempty(aseg(ct).diamu)))
        mdiam(ct)=NaN;
    else
        mdiam(ct)=mean(aseg(ct).diamd(aseg(ct).seld==1))+mean(aseg(ct).diamu(aseg(ct).selu==1));
    end
end

ind_nan=find(not(isnan(mdiam)));
%% Remove linear trend
t=[1:length(mdiam)]-1;
pl=polyfit(t(ind_nan),mdiam(ind_nan),2);
mfit=polyval(pl,t(ind_nan));
mdiam_detrend=mdiam(ind_nan)-mfit;

torig=t;
t=t(ind_nan);

%% Estimate cycle
fdiam=(fft(mdiam_detrend));
[fmax,pmax]=max(abs(fdiam(1:fix(length(t)/2))));

f0=pmax/length(t);
phi0=phase(fdiam(pmax));

beta0 = [f0,phi0];
f = @(beta_x) aIMTcardiac_GT2(beta_x,t,mdiam_detrend);

options = saoptimset('simulannealbnd');
options.TolFun = 1e-12;

[beta_opt,fval,exitFlag,output] = simulannealbnd(f,beta0,[0 0],[0.5 0.5],options);

%% Estimating amplitude
f_opt = cos(2*pi*beta_opt(1)*t - 2*pi*beta_opt(2));
y=mdiam_detrend;
w = abs(y)/sum(abs(y));
beta1_opt = lscov(f_opt',y',sqrt(w'));

figure()
plot(torig,mdiam)
hold on
plot(t,beta1_opt*f_opt+mfit)


%% Stiffnes using mean amplitude

%median_lumen_width=mean(mdiam);
mean_dist_lumen_stim = beta1_opt * f_opt + mfit;

lumen_stiff(1) = 2*abs(beta1_opt)./nanmedian(mdiam);

%% Stiffness estimating maxs and mins
T = 1/beta_opt(1);
if beta1_opt > 0
    t_maxs = beta_opt(2)/beta_opt(1):T:t(end);
    t_mins = (beta_opt(2)/beta_opt(1) + T/2):T:t(end);
else
    t_maxs = (beta_opt(2)/beta_opt(1) + T/2):T:t(end);
    t_mins = beta_opt(2)/beta_opt(1):T:t(end);
end

if length(t_maxs)~=length(t_mins)
    t_maxs = t_maxs(1:min([length(t_maxs) length(t_mins)]));
    t_mins = t_mins(1:min([length(t_maxs) length(t_mins)]));
end

n_vals = 1;

ind_maxs = zeros(1,length(t_maxs));
ind_mins = zeros(1,length(t_mins));
maxs = zeros(1,length(t_maxs));
mins = zeros(1,length(t_mins));

for k = 1:length(t_maxs)
    [~, ind_maxs(k)] = min(abs(t-t_maxs(k)));
    [~, ind_mins(k)] = min(abs(t-t_mins(k)));

    aux_n = n_vals;
    while true
        if (ind_maxs(k)-aux_n>=1) && (ind_maxs(k)+aux_n<=length(mean_dist_lumen_stim))
            maxs(k) = nanmean(mdiam(ind_maxs(k)-aux_n:ind_maxs(k)+aux_n));
            break
        else
            aux_n = aux_n-1;
        end
    end

    aux_n = n_vals;
    while true
        if (ind_mins(k)-aux_n>=1) && (ind_mins(k)+aux_n<=length(mean_dist_lumen_stim))
            mins(k) = nanmean(mdiam(ind_mins(k)-aux_n:ind_mins(k)+aux_n));
            break
        else
            aux_n = aux_n-1;
        end
    end

end

lumen_diffs = abs(maxs-mins);
lumen_means = mean([maxs; mins]);
lumen_stiff(2) = mean(lumen_diffs./lumen_means);

diam_sys=nanmedian(mins);
diam_dia=nanmedian(maxs);
diam_med=nanmedian(mdiam);

fprintf('--------------------------------------------- \n Results for file: %s \n ---------------------------------------------\n',filename1);
fprintf('Fetal Abdominal Aorta Measurements\n')
fprintf('\n\t Diameters: \n')
fprintf('\t\t Median diameter: %.4f\n',diam_med)
fprintf('\t\t Median systolic diameter: %.4f\n',diam_sys)
fprintf('\t\t Median diastolic diameter: %.4f\n',diam_dia)
fprintf('\n\t Stiffnes: \n')
fprintf('\t\t Stiffness using data range: %.4f\n',lumen_stiff(2))
fprintf('\t\t Stiffness using model range: %.4f\n',lumen_stiff(1))
fprintf('--------------------------------------------- \n')


