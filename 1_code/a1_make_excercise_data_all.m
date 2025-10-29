clc                      
clear           
%% Config
ups_dir = "E:\PM_foraging-performance\0_data"
out_Dir = 'E:\PM_foraging-performance\2_pipeline\make_excersise_data_all';
alpha=25; %degrees, cuts off magnetometry when the sensor is too closely aligned with the earth's magnetic field

%% LOAD DATA
load(fullfile(ups_dir,'SW_tab_all.mat'));
tab_all.caldir(:) = {'E:\Shared_raw\3S\cal\'}
tab_all.prhdir(:) = {'E:\Shared_raw\3S\prh\'}
tab_all.audiodir(:) = {''}
save(fullfile(ups_dir,'SW_tab_UH.mat'), 'tab_all')
writetable(tab_all, fullfile(ups_dir,'SW_tab_UH.csv'))
load(fullfile(ups_dir,'SW_tab_UH.mat'));


%in for loop
for i = 13:45
tagid = char(tab_all.ind(i))
calDir = char(tab_all.caldir(i));
prhDir = char(tab_all.prhdir(i));
settagpath('cal', calDir, 'prh', prhDir, 'audio', '')

[Aw, M, p, pitch, roll, head, t, tt, fs, tot, Aframe] = loadPrhData_D3_AB(tagid, prhDir, calDir, 'sw', '5');
    
% check for slips
[CAL,DEPLOY,ufname] = d3loadcal(tagid)
DEPLOY.OTAB
%trim
t_off = tab_all.tagoff(i)
Aw = Aw(1:round(t_off*fs),:);
M= M(1:round(t_off*fs),:);
p = (p(1:round(t_off*fs)));


%%
FR = tab_all.flukefs(i)
cutoff = tab_all.flukefs_lo(i)
f = cutoff/FR
J = 0.03
tmax = (1/FR)*1.5
k=1:length(p);


[MagAcc,pry,Sa,GLm,KKm] = magnet_rot_sa(Aw,M,fs,FR,f,alpha,1,k,J,tmax);
OSA = sum(abs(Sa),2); %OSA per Martin Lopez et al.

sp1 = subplot(5,1,1:4);
scatter(k,p,[],log(OSA),'filled')
sp2= subplot(5,1,5);
ax2= plot(k,pry(:,1));
hold on
plot(KKm(:,1)*fs,pry(round(KKm(:,1)*fs),1),'r*')
linkaxes([sp1 sp2 ax2], 'x'); % links x axes
hold off
saveas(gcf,fullfile(out_Dir,join([tagid,'_exc_5hz.fig'])))
saveas(gcf,fullfile(out_Dir,join([tagid,'_exc_5hz.png'])))

KKm_gmt = 

OSA_1hz = mean(reshape([OSA(:); nan(mod(-numel(OSA),fs),1)],fs,[]), "omitnan")'; % OSA @1hz (block averaging)
p_1hz = p(1:fs:length(p)); %depth at 1hz
Aw_1hz = Aw(1:fs:length(p)); %whale-frame acceleration at 1hz
Mw_1hz = M(1:fs:length(p)); % whale-frame magnetometry at 1hz

save(fullfile(out_Dir,join([tagid,'_exc_5hz.mat'])),"OSA","p")
end


% now for D2s
for i = 1 :12
tagid = char(tab_all.ind(i))
calDir = char(tab_all.caldir(i));
prhDir = char(tab_all.prhdir(i));
settagpath('cal', calDir, 'prh', prhDir, 'audio', '')

[Aw, M, p, pitch, roll, head, t, tt, fs, tot, Aframe]=loadPrhData(tagid, prhDir, calDir, 'sw');

    
    
% check for slips
loadcal(tagid)
OTAB


%%
FR = tab_all.flukefs(i)
cutoff = tab_all.flukefs_lo(i)
f = cutoff/FR
J = 0.03
tmax = (1/FR)*1.5
k=1:length(p);

[MagAcc,pry,Sa,GLm,KKm] = magnet_rot_sa(Aw,M,fs,FR,f,alpha,1,k,J,tmax);
OSA = sum(abs(Sa),2); %OSA per Martin Lopez et al.


%trim
t_off = min(round(tab_all.tagoff(i))*fs,length(p))
Aw = Aw(1:t_off,:);
M= M(1:t_off,:);
pry = pry(1:t_off,:);
p = p(1:t_off);
OSA = OSA(1:t_off);
KKm = KKm(KKm(:,1)*fs < t_off);

k=1:length(p);

sp1 = subplot(5,1,1:4);
scatter(k,p,[],log(OSA),'filled')
sp2= subplot(5,1,5);
ax2= plot(k,pry(:,1));
hold on
plot(KKm(:,1)*fs,pry(round(KKm(:,1)*fs),1),'r*')
linkaxes([sp1 sp2 ax2], 'x'); % links x axes
hold off
saveas(gcf,fullfile(out_Dir,join([tagid,'_exc_5hz.fig'])))
saveas(gcf,fullfile(out_Dir,join([tagid,'_exc_5hz.png'])))



OSA_1hz = mean(reshape([OSA(:); nan(mod(-numel(OSA),fs),1)],fs,[]), "omitnan")'; % OSA @1hz (block averaging)
p_1hz = p(1:fs:length(p)); %depth at 1hz
Aw_1hz = Aw(1:fs:length(p)); %whale-frame acceleration at 1hz
Mw_1hz = M(1:fs:length(p)); % whale-frame magnetometry at 1hz




save(fullfile(out_Dir,join([tagid,'_exc_5hz.mat'])),"OSA","p")
end