function Fcn_Palies_burner_shape

L1Sp = [100,164,228]-4;
L3Sp = [100,150,200,400];

for ss = 1:length(L1Sp)
    L1 = L1Sp(ss);
    for kk = 1:length(L3Sp)
        L3 = L3Sp(kk);
        Fcn_create_txt_file(L1,L3)
    end
end

function Fcn_create_txt_file(L1,L3)
% ------------------------------
% plenum information
xPlenumSp = [- L1, 0];
DPlenumSp = [65,65];
x_sample(1:length(xPlenumSp)-1)     = xPlenumSp(1:length(xPlenumSp)-1)./1e3;
r_sample(1:length(xPlenumSp)-1)     = DPlenumSp(1:length(xPlenumSp)-1)./2./1e3;
SectionIndex(1:length(xPlenumSp)-1) = 0;
TubeIndex(1:length(xPlenumSp)-1)    = 0;
% ------------------------------
% 
DRod = 6;                                           % diameter of the rod
%
% -------------------------------
% convergent nozzle
LConv           = 60;                               % length of convergent nozzle
xConvSp         = linspace(0,LConv,61);            % the entrance of the nozzle is set to '0'
DConv_inlet     = (65.^2 - DRod.^2).^0.5;           % equivalent diameter of the entrance of the nozzle
DConv_outlet    = (30.^2 - DRod.^2).^0.5;           % equivalent diameter of the outlet of the nozzle
Ratio           = DConv_inlet/DConv_outlet;         % a ratio
DConvSp         = DConv_inlet.*(1 - (1-Ratio.^2)*(1-((LConv - xConvSp)./LConv).^2).^2./(1 + (LConv - xConvSp).^2./(3*LConv.^2)).^3).^(-0.5);
%
x_sample        = cat(2,x_sample, xConvSp(1:end-1)./1e3);
r_sample        = cat(2,r_sample, DConvSp(1:end-1)./1e3./2);
SectionIndex    = cat(2,SectionIndex,0.*xConvSp(1:end-1));
TubeIndex       = cat(2,TubeIndex,0.*xConvSp(1:end-1)+2);
%
% --------------------------------
% adapter, which is consider to be diameter-linearly-gradually-increased 
DInject         = (22.^2 - DRod.^2).^0.5;
xAdaptSp        = linspace(LConv, LConv+20,2);
DAdaptSp        = linspace(DInject, DInject, 2);
%
x_sample        = cat(2,x_sample, xAdaptSp(1:end-1)./1e3);
r_sample        = cat(2,r_sample, DAdaptSp(1:end-1)./1e3./2);
SectionIndex    = cat(2,SectionIndex,0.*xAdaptSp(1:end-1));
TubeIndex       = cat(2,TubeIndex,0.*xAdaptSp(1:end-1));
% --------------------------------
% injector
xInjectorSp     = [LConv+20,LConv+20+56];
DInjectorSp     = DInject*[1,1];
%
x_sample        = cat(2,x_sample, xInjectorSp(1:end-1)./1e3);
r_sample        = cat(2,r_sample, DInjectorSp(1:end-1)./1e3./2);
SectionIndex    = cat(2,SectionIndex,0.*xInjectorSp(1:end-1));
TubeIndex       = cat(2,TubeIndex,0.*xInjectorSp(1:end-1));
% --------------------------------
% chamber
xChamberSp      = [LConv+20+56, LConv+20+56 + L3];
DChamberSp      = 70*[1,1];
%
x_sample        = cat(2,x_sample, xChamberSp(1:end)./1e3);
r_sample        = cat(2,r_sample, DChamberSp(1:end)./1e3./2);
SectionIndex    = cat(2,SectionIndex,0.*xChamberSp(1:end));
TubeIndex       = cat(2,TubeIndex,0.*xChamberSp(1:end));
SectionIndex(end-1) = 11;
% --------------------------------
% figure
% plot(x_sample,r_sample)

% figure
% hold on
% plot(xPlenumSp,DPlenumSp./2,'-k')
% plot(xConvSp,DConvSp./2,'-k')
% plot(xAdaptSp,DAdaptSp./2,'-k')
% plot(xInjectorSp,DInjectorSp./2,'-k')
% plot(xChamberSp,DChamberSp./2,'-k')
% plot([xChamberSp(1),xChamberSp(1)],[DInjectorSp(end), DChamberSp(1)]./2,'-k')
% % -----
% plot(xPlenumSp,-DPlenumSp./2,'-k')
% plot(xConvSp,-DConvSp./2,'-k')
% plot(xAdaptSp,-DAdaptSp./2,'-k')
% plot(xInjectorSp,-DInjectorSp./2,'-k')
% plot(xChamberSp,-DChamberSp./2,'-k')
% plot([xChamberSp(1),xChamberSp(1)],-[DInjectorSp(end), DChamberSp(1)]./2,'-k')


% ---------------------------------
data        = cat(1,x_sample,r_sample,SectionIndex,TubeIndex);
% ---------------------------------
currentFolder   = pwd;
filename = ['CD_Palies_Real_Shape_L1_' num2str(L1) '_L3_' num2str(L3) '.txt'];
currentFolder   = fullfile(currentFolder,filename);
fid             = fopen(currentFolder,'wt');
data_title      = {'x[m]','r[m]','SectionIndex','TubeIndex'};
fprintf(fid,'%s\b',data_title{1:end});
fprintf(fid,'\n');
fprintf(fid,'%6.5f\b%6.5f\b%6.0f\b%6.0f\b\n',data);
fclose(fid);
% -----------------------------end-----------------------------------------
