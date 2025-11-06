%% process D3 data

species = 'sw';
ups_dir = 'E:\PM_foraging-performance\0_data';
downs_dir = 'E:\PM_foraging-performance\2_pipeline';
load(fullfile(ups_dir,'SW_tab_UH.mat'));
bin_len = 30 %in s, for turning angle
func = @(x) circ_mean(x)

thdeg = 30;
pitch_thresh = 80;
tab = tab_all;
dep_tim_3 = struct; % 4 dimensional data structure [deploy, dive, t, p] for plotting
divetab_3 = cell2table(cell(0,16), 'Variablenames', {'Deploy', 'Dive', 'Startcue', 'Endcue' , 'Max_depth', 'Surf_int', 'Head_var','Kappa', 'Tagoff', 'Turn_ang', 'H_tort', 'res_ind', 'Start_tasc','Dstart_tasc','check', 'Cross_20'});
replace = true;
alpha=25; %degrees, cuts off magnetometry when the sensor is too closely aligned with the earth's magnetic field
mindivedef = 10;
%df = 1;%working decimation factor
plot = true;


for i = 13:45
    tagid = char(tab.ind(i))
    calDir = char(tab.caldir(i));
    prhDir = char(tab.prhdir(i));
    settagpath('cal', calDir, 'prh', prhDir, 'audio', '')
    
    if strcmp(tab.ind(i), 'sw16_126a') 
        'skipped'
        continue;
    end
    
    %clf
    [Aw, M, p, pitch, roll, head, t, tt, fs, tot, Aframe] = loadPrhData_D3_AB(tagid, prhDir, calDir, species, '5');
    
    %Aw = Aw(1:df:length(Aw),:)
    
    if strcmp(tab.ind(i), 'sw19_243a') | strcmp(tab.ind(i), 'sw19_255d')
        load(['\\cfs.st-andrews.ac.uk\shared\smru-bk\MillerLab\Final\3S_2019\DTAG\prh\Depth_updated\', char(tab.ind(i)),'.mat'])
        p = depth;
        'depth corrected'
    end
    
    dep_tim_3(i).deploy = tagid;
    dep_tim_3(i).dives = {};
    
    FR = tab.flukefs(i);
    cutoff = tab.flukefs_lo(i);
    k=1:length(p);
    f = cutoff/FR;
    
    J = 0.03
    tmax = (1/FR)*1.5
    k=1:length(p);
    [MagAcc,pry,Sa,GLm,KKm] = magnet_rot_sa(Aw,M,fs,FR,f,alpha,1,k,J,tmax);
    OSA = sum(abs(Sa),2); %OSA per Martin Lopez et al.
    ODBA = odba(Aw,5,f);

    OSA_1hz = mean(reshape([OSA(:); nan(mod(-numel(OSA),fs),1)],fs,[]), "omitnan")'; % OSA @1hz (block averaging)
    ODBA_1hz = mean(reshape([ODBA(:); nan(mod(-numel(ODBA),fs),1)],fs,[]), "omitnan")'; % OSA @1hz (block averaging)
    
    save(fullfile(downs_dir,join([tagid,'_perf.mat'])),"OSA_1hz","ODBA_1hz","KKm")
    
    [Anlf,~,~,~] = Ahf_Anlf(Aw,fs,FR,f,3,k,[],[]) ;
    [smoothpitch,~]=a2pr(Anlf(k,:));
    %speed = inst_speed(p,smoothpitch,fs,FR,f,k,thdeg);
    
    T = finddives(p,fs,mindivedef,2,0);
    itab = array2table(NaN(size(T,1),17), 'Variablenames',{'Dive','Startcue', 'Endcue' , 'Max_depth','Surf_int', 'Head_var','Kappa', 'Turn_ang', 'Tagoff', 'H_tort', 'res_ind', 'End_desc', 'D_end_desc', 'Start_tasc', 'Dstart_tasc', 'check', 'Cross_20'});
    itab.Deploy = repmat(cellstr(tagid), size(T,1),1);
    itab.Startcue = T(1:size(T,1),1);
    itab.Endcue = T(1:size(T,1),2);
    itab.Max_depth = T(1:size(T,1),3);
    itab.Surf_int = [T(2:end,1)-T(1:end-1,2);NaN];
    % [start_cue end_cue max_depth cue_at_max_depth mean_depth mean_compression]
    hold off
    for dive = 1:size(T,1)
        %dive
        itab.Dive(dive) = dive;
        if tab.tagoff(i)< T(dive, 2) %skip the dive if the tag falls off 
            itab.Tagoff(dive) = 1;
            itab.Max_depth(dive) = NaN;
        else
            itab.Tagoff(dive) = 0;
            kk=round(fs*T(dive,1)):round(fs*T(dive,2));
            kk2=round(fs*T(dive,1)):round(fs*(T(dive,2)-30));
            pdive = p(kk);
            %dive_speed = mean(speed(kk),'omitnan');
            tsd = (kk/fs) - T(dive,1); %time since dive start;
            %  (in samples):
            cross_20 = find(pdive > 20, 1, 'first') + T(dive,1)*fs;
            if length(cross_20) == 1
                itab.Cross_20(dive)=cross_20/fs; %in seconds
            end
            startasc=round(find((smoothpitch(kk&p>mindivedef)*180/pi)<0,1,'last')+T(dive,1)*fs);%search for the last point beroe diving above mindivedef at which the pitch is negative
            if length(tasc) == 1
                itab.Start_tasc(dive)=startasc/fs; %in seconds
                itab.Dstart_tasc(dive)= p(round(startasc));
            end
            enddes=round(find((smoothpitch(kk&p>mindivedef)*180/pi)>0,1,'first')+T(dive,1)*fs);% search for the first point after diving below mindivedef at which pitch is positive
            if length(enddes) == 1
                itab.End_desc(dive)=enddes/fs; %in seconds
                itab.D_end_desc(dive)= p(round(enddes));
            end
            if p(round(tasc)) / itab.Max_depth(dive) < 0.1
                itab.check(dive) = true;
                 p(round(tasc))
                 figure(1);
                 sp1 = subplot(5,1,1:4);
                 plott(pdive, fs);
                 xline(((itab.Start_tasc(dive)*fs)-(T(dive,1)*fs))/fs);
                 xline(((itab.Cross_20(dive)*fs)-(T(dive,1)*fs))/fs);
                 sp2=subplot(5,1,5);
                 plott(pitch(kk),fs);
                 hold on
                 plott(smoothpitch(kk),fs);
                 hold off
                 if replace == true
                 [new_tasc, ~] = ginput(1);
                 new_tasc = new_tasc + T(dive,1);
                 itab.Start_tasc(dive) = new_tasc;
                 itab.Dstart_tasc(dive) = p(round(itab.Start_tasc(dive)*fs));
                 end 
            end
            dep_tim_3(i).dives{dive} = [tsd', pdive]; % time series for plotting dives in R
            if ~strncmp('sw16_', tagid, 4)
                pseud_track = ptrack(pitch(kk),head(kk),p(kk),fs,cutoff);
                tort = tortuosity(pseud_track(:,1:2), fs,floor(length(pseud_track)/fs));
                itab.H_tort(dive) = tort(1,1);
                RI = residence_index(pseud_track, fs, 1800,50);
                itab.res_ind(dive) = mean(RI);
                headdive = (head(kk));
                pitchdive = (pitch(kk));
                itab.Head_var(dive) = circ_var(headdive);
                itab.Kappa(dive) = circ_kappa(headdive);
                if T(dive,2)-T(dive,1)>bin_len*2 % if dive is long enough
                    edges = (floor(min(kk)):bin_len*fs:max(kk));
                    bins = discretize(kk,edges);
                    avg_head = splitapply(func, headdive, bins');
                    max_pitch = splitapply(@max, pitchdive, bins');
                    avg_head(rad2deg(abs(max_pitch))>pitch_thresh) = NaN; % exclude gimbal-locked bins
                    avg_head_prev = [NaN; avg_head];
                    avg_head = [avg_head;NaN];
                    turn_ang = angdiff(avg_head,avg_head_prev);
                    itab.Turn_ang(dive) = (mean(abs(turn_ang(~isnan(turn_ang)))))*180/pi;
                end
                %
            end
            %if i/3 == round(i/3)
            %plott(smoothpitch(kk),fs)
            %plott(headdive, fs)
            %dive
            %end
        end
    end
    divetab_3 = [divetab_3;itab];
    %input('next?')
end

save(fullfile(downs_dir, 'SW_divetab_d3.mat'), 'divetab_3')
save(fullfile(downs_dir, 'SW_depth-time_d3.mat'), 'dep_tim_3')

%%  D2 data
clear all
species = 'sw';
ups_dir = 'D:\Data\Sonar_response\';
downs_dir = ups_dir;
load(fullfile(ups_dir,'SW_tab_all.mat'));
func = @(x) circ_mean(x)

bin_len = 30 %in s, for turning angle
thdeg = 30
pitch_thresh = 80
tab = tab_all;
dep_tim_2 = struct; % 4 dimensional data structure [deploy, dive, t, p] for plotting
divetab2 = cell2table(cell(0,16), 'Variablenames', {'Deploy', 'Dive', 'Startcue', 'Endcue' , 'Max_depth', 'Surf_int', 'Head_var','Kappa', 'Tagoff', 'Turn_ang', 'H_tort', 'res_ind', 'Start_tasc','Dstart_tasc','check', 'Cross_20' });
tab = tab_all;
% tagid = 'sw16_126a'
replace = true;
%hold on
for i = 1:12
    tagid =char(tab.ind(i))
    calDir = char(tab.caldir(i));
    prhDir = char(tab.prhdir(i));
    settagpath('cal', calDir, 'prh', prhDir, 'audio', 'D:\Analysis\SW sonar chapter\make_dive_data')
    %clf
    [Aw, M, p, pitch, roll, head, t, tt, fs, tot, Aframe]=loadPrhData(tagid, prhDir, calDir, species);
    
    dep_tim_2(i).deploy = tagid;
    dep_tim_2(i).dives = {};
    
    FR = tab.flukefs(i);
    cutoff = tab.flukefs_lo(i);
    k=1:length(p);
    f = cutoff/FR;
    
    [Anlf,~,~,~] = Ahf_Anlf(Aw,fs,FR,f,3,k,[],[]) ;
    [smoothpitch,~]=a2pr(Anlf(k,:));
    
    T = finddives(p,fs,10,2,0);
    itab = array2table(NaN(size(T,1),15), 'Variablenames',{'Dive','Startcue', 'Endcue' , 'Max_depth','Surf_int', 'Head_var','Kappa', 'Turn_ang', 'Tagoff', 'H_tort', 'res_ind', 'Start_tasc', 'Dstart_tasc', 'check', 'Cross_20'});    
    itab.Deploy = repmat(cellstr(tagid), size(T,1),1);
    itab.Startcue = T(1:size(T,1),1);
    itab.Endcue = T(1:size(T,1),2);
    itab.Max_depth = T(1:size(T,1),3);
    itab.Surf_int = [T(2:end,1)-T(1:end-1,2);NaN];
    % [start_cue end_cue max_depth cue_at_max_depth mean_depth mean_compression]
    hold off
    for dive = 1:size(T,1)
        itab.Dive(dive) = dive;
        if tab.tagoff(i)< T(dive, 2) %skip the dive if the tag falls off.
            itab.Tagoff(dive) = 1;
            itab.Max_depth(dive) = NaN;
        else
            itab.Tagoff(dive) = 0;
            kk=round(fs*T(dive,1)):round(fs*T(dive,2));
            kk2=round(fs*T(dive,1)):round(fs*(T(dive,2)-30));
            pdive = p(kk);
            %smoothpitch_dive = smoothpitch(kk);
            %dive_speed = mean(speed(kk),'omitnan');
            tsd = (kk/fs) - T(dive,1); %time since dive start;
            %start of terminal ascent (absolute time in seconds):
            cross_20 = find(pdive > 20, 1, 'first') + T(dive,1)*fs;
            if length(cross_20) == 1
                itab.Cross_20(dive)=cross_20/fs; %in seconds
            end
            tasc = find((smoothpitch(kk2)*180/pi)<-5,1,'last')+ T(dive,1)*fs;
            if length(tasc) == 1
                itab.Start_tasc(dive)=tasc/fs;
                itab.Dstart_tasc(dive)= p(round(tasc));
            end
            if p(round(tasc)) / itab.Max_depth(dive) < 0.1;
                itab.check(dive) = true;
                p(round(tasc))
                figure(1);
                sp1 = subplot(5,1,1:4);
                plott(pdive, fs);
                xline(((itab.Start_tasc(dive)*fs)-(T(dive,1)*fs))/fs);
                sp2=subplot(5,1,5);
                plott(pitch(kk),fs);
                hold on
                plott(smoothpitch(kk),fs)
                hold off
                 if replace == true
                 [new_tasc, ~] = ginput(1);
                 new_tasc = new_tasc + T(dive,1);
                 itab.Start_tasc(dive) = new_tasc;
                 itab.Dstart_tasc(dive) = p(round(itab.Start_tasc(dive)*fs));
                 end 
            end
            dep_tim_3(i).dives{dive} = [tsd', pdive]; % time series for plotting dives in R
            if ~strncmp('sw16_', tagid, 4)
                pseud_track = ptrack(pitch(kk),head(kk),p(kk),fs,cutoff);
                tort = tortuosity(pseud_track(:,1:2), fs,floor(length(pseud_track)/fs));
                itab.H_tort(dive) = tort(1,1);
                RI = residence_index(pseud_track, fs, 1800);
                itab.res_ind(dive) = mean(RI);
                headdive = (head(kk));
                pitchdive = (pitch(kk));
                itab.Head_var(dive) = circ_var(headdive);
                itab.Kappa(dive) = circ_kappa(headdive);
                if T(dive,2)-T(dive,1)>bin_len*2 % if dive is long enough
                    edges = (floor(min(kk)):bin_len*fs:max(kk));
                    bins = discretize(kk,edges);
                    avg_head = splitapply(func, headdive, bins');
                    max_pitch = splitapply(@max, pitchdive, bins');
                    avg_head(rad2deg(abs(max_pitch))>pitch_thresh) = NaN; % exclude gimbal-locked bins
                    avg_head_prev = [NaN; avg_head];
                    avg_head = [avg_head;NaN];
                    turn_ang = angdiff(avg_head,avg_head_prev);
                    itab.Turn_ang(dive) = (mean(abs(turn_ang(~isnan(turn_ang)))))*180/pi;
                end
                %
            end
            %plott(pdive, fs)
            %plott(headdive, fs)
            %dive
        end
        %plott(pdive, fs)
    end
    divetab2 = [divetab2;itab];
    %input('next?')
    %clf
end

save(fullfile(downs_dir, 'SW_divetab_d2.mat'), 'divetab2')
save(fullfile(downs_dir, 'SW_depth-time_d2.mat'), 'dep_tim_2')


%%  combine & save
clear all
ups_dir = 'D:\Data\Sonar_response\';
downs_dir = ups_dir;
load(fullfile(downs_dir, 'SW_divetab_d2.mat'));
load(fullfile(downs_dir, 'SW_divetab_d3.mat'));
SW_divetab_all = [divetab2;divetab_3];
save(fullfile(downs_dir, 'SW_divetab_all.mat'), 'SW_divetab_all')
writetable(SW_divetab_all, fullfile(downs_dir,'SW_divetab_all.csv'))



