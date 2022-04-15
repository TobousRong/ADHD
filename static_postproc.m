clc;clear;close all
%% Static analysis: This function calculates the calibrated Hin, Hse and HB from fMRI time series
[~,Subj]=xlsread('E:\ADHD200\ADHD.xlsx','data','a2:a195');
N=100;N_sub=194;
%mypool=parpool('local',24,'IdleTimeout',240);
Long_fmri=[];IN=[];IM=[];rIN=[];rIM=[];
parfor sub=1:N_sub
    path=strcat('E:\ADHD200\output\FC\',Subj(sub),'_FC_combat0.mat');
    MRI=load(char(path));
    FC=abs(MRI.FC);
    [Clus_num,Clus_size,mFC] = Functional_HP(FC,N);
    [Hin,Hse,R_Hin,R_Hse,Hin_inter,Hse_inter] =Seg_Int_component(FC,N,Clus_size,Clus_num);
    IN=[IN;Hin];IM=[IM;Hse];rIN=[rIN,R_Hin'];rIM=[rIM,R_Hse'];
    file=strcat(Subj(sub),'_HF_combat0_abs.mat');
    par_save(file,Hin,Hse,R_Hin,R_Hse,Hin_inter,Hse_inter)
end
load('adhd_stable_FC.mat');
[adhd_Hin,adhd_Hse] = Stable_correct(abs(sFC),IN(1:97),IM(1:97),N);
load('hc_stable_FC.mat');
[hc_Hin,hc_Hse] = Stable_correct(abs(sFC),IN(98:end),IM(98:end),N);
Hin=[adhd_Hin;hc_Hin];Hse=[adhd_Hse;hc_Hse];
rHin=[];rHse=[];
for sub=1:194
    rHin=[rHin,rIN(:,sub)*Hin(sub)/mean(rIN(:,sub))];
    rHse=[rHse,rIM(:,sub)*Hse(sub)/mean(rIM(:,sub))];
end
save('Combat0_abs.mat','Hin','Hse','rHin','rHse')

