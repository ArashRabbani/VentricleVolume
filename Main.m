function Main
% Please first run Main.py and then Main.m
close all
load('seg.mat'); %  Loading the segmentated geometry calculated by Main.py code. 
Meta.RES=[1.2,1.2,10]; % spatial resolution [dx,dy,dz] (mm/pixel)
Meta.ZLOC=[1.01:10:100]; % location of the short axis slices in the vertical direction (mm)
% segmentation labels [0,1,2,3]=['RV','MY','LV','BG']=['Right Ventricle','Myocardium','Left Ventricle','Background']
LabList={'LV','MY','RV'}; % list of labels
Result1=volcalc(AAG,Meta,'LV'); % calculates the volume of LV using gaussian process
Result2=volcalc(AAG,Meta,'MY'); % calculates the volume of MY using gaussian process
Result3=volcalc(AAG,Meta,'RV'); % calculates the volume of RV using gaussian process
figure; hold on; % Plotting the volumes
plot(Result1.Time,Result1.VolCurve);
plot(Result2.Time,Result2.VolCurve);
plot(Result3.Time,Result3.VolCurve);
xlabel('Dimentionless time'); ylabel('Volume (mL)'); legend(LabList); box on;
MakeAnimation(AAG,Result1,Result2,Result3,LabList); % Making a gif animation of volume changes
% save the results as an excel file: 'Results.xlsx'
data_cells=num2cell([Result1.Time',Result1.VolCurve,Result2.VolCurve,Result3.VolCurve]);
col_header={'Dimentionless time','Left ventricle volume (ml)','Myocardium volume (ml)','Right ventricle volume (ml)'};
data_cells=[ col_header; data_cells];  writecell(data_cells,'Results.xlsx')
end
function make_video_gif
filename='Final.gif';
for I = 1:numel(dir('Export/*.png'))
    A=imread(['Export/F' num2str(I) '.png']);
    [imind,cm] = rgb2ind(A,256);
    if I == 1;
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.05);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.05);
    end
end
end
function MakeAnimation(AAG,Result1,Result2,Result3,LabList)
figure;
load('Raw.mat')
for I=1:size(AAG,4)
    clf; hold on;
    raw=squeeze(M(:,:,4,I));
    [raw,Rat]=prep(raw);
    ax(1)=subplot(2,2,1);
    imagesc(raw); colormap gray; axis equal tight off;
    ax(2)=subplot(2,2,3);
    imagesc(squeeze(AAG(4,:,:,I))); colormap parula; axis equal tight off;
    subplot(2,2,[2,4]);
    ID=round(I/size(AAG,4)*numel(Result3.VolCurve));
    b=bar([Result1.VolCurve(ID),Result2.VolCurve(ID),Result3.VolCurve(ID)]); ylim([0,max(Result1.VolCurve)*1.2])
    set(gca,'xticklabel',LabList); ylabel('Volume (mL)');
    colormap(ax(1),gray);
    colormap(ax(2),parula);
    CM=parula(4);
    b.FaceColor = 'flat'; b.CData=CM(1:3,:);
    drawnow;
    printx(['Export/F' num2str(I) ])
end
make_video_gif;
end
function []=printx(Name)
set(gca, 'units', 'normalized');
print([Name],'-dpng','-r200');
end
function out=volcalc(AAG,Meta,Lab)
LabList={'RV','MY','LV'};
if Lab==LabList{1}; Val=0; end
if Lab==LabList{2}; Val=1; end
if Lab==LabList{3}; Val=2; end
II=1;
show1=0;
show2=0;
show3=0;
RES=Meta.RES;
try
    ZLOC=Meta.ZLOC;
catch
    try
        Meta.ZLOC= [0.01:Meta.RES(3):300];
        Meta.ZLOC=Meta.ZLOC(1:size(AAG,1));
        ZLOC=Meta.ZLOC;
    catch
        error('You should provide either vertical spatial resolution or locations of the slices.')
    end
end
Frames=size(AAG,4);
SLs=size(AAG,1);
for I=1:Frames
    for J=1:SLs
        V(I,J,1)=sum_(AAG(J,:,:,I)==Val);
    end
end
v=cat(1,V(:,:,1),V(:,:,1),V(:,:,1));
v=v.*RES(II,1).*RES(II,2)./100; % cross-sectional area in cm^2
v=imgaussfilt(v,1,'Padding','replicate');
ZL=ZLOC(II,:); ZL(ZL==0)=[]; ZL=ZL-min(ZL); ZL=ZL./10; if ZL(1)>ZL(end); ZL=max(ZL)-ZL; end
if show3
    figure; hold on
    for J=1:SLs
        B=squeeze((AAG(J,:,:,TF(II))==2));
        [x,y]=find(B); x=x.*RES(II,1)./10; y=y.*RES(II,2)./10;
        x=x-mean(x); y=y-mean(y);
        z=-repmat(ZL(J),numel(x),1);
        scatter3(x,y,z,'MarkerEdgeColor','b'); view(3)
    end
    try ;load(['../Data/' D(II).name(1:end-4) '/abaqusInputData.mat']);
        endoNodes = abaqusInputData.endoNodes;
        node = abaqusInputData.node(:,1:3);
        nodes_x = node(endoNodes,1);
        nodes_y = node(endoNodes,2);
        nodes_z = node(endoNodes,3);
        scatter3(nodes_x,nodes_y,nodes_z,'MarkerEdgeColor','k');
    catch; end
    axis equal; drawnow;
end
acc=2;
v(:,end+1)=v(:,end)-mean(abs(v(:,end)-v(:,end-1)))*acc; ZL(end+1)=ZL(end)+(ZL(end)-ZL(end-1));
v(:,end+1)=v(:,end)-mean(abs(v(:,end)-v(:,end-1)))*acc; ZL(end+1)=ZL(end)+(ZL(end)-ZL(end-1));
v(:,end+1)=v(:,end)-mean(abs(v(:,end)-v(:,end-1)))*acc; ZL(end+1)=ZL(end)+(ZL(end)-ZL(end-1));
[x,y]=meshgrid(ZL,1:size(v,1)); x=x(:); y=y(:); z=v(:);
gprMdl = fitrgp([x,y],z);
dx=.1;
dy=.2;
xbin=0:dx:max(ZL);
ybin=1:dy:size(v,1);
[x2,y2]=meshgrid(xbin,ybin); x2=x2(:); y2=y2(:);
zpred = predict(gprMdl,[x2,y2]);
z0=zpred;
if show1
    figure; scatter3(x,y,z,'MarkerEdgeColor','red'); hold on
    p=find(zpred<0);
    zpred(p)=[]; x2(p)=[];y2(p)=[];
    scatter3(x2,y2,zpred,'MarkerEdgeColor','blue','SizeData',2); pbaspect([2 3 1])
    zlim([0,max_(zpred)]); ylabel('Time frame'); xlabel('Slice location from base (cm)'); zlabel('Cross-section Area (cm^2)');
    drawnow;
end
sum_(zpred);
z0=reshape(z0,[numel(ybin),numel(xbin)]);
z0(z0<0)=0;
t=sum(z0,2);
z_bad=v; z_bad(z_bad<0)=0; t_bad=sum(z_bad,2);
m1=min(t);
m2=max(t);
if show2
    timestamps=(1:numel(t)).*dy;
    figure; hold on; plot(timestamps,t.*dx); ylabel('LV volume (ml)'); xlabel('Time frames');
    
end
ESV(II)=m1.*dx;
[Feat,FeatNames]=Features(AAG,RES(II,1:2));
X(II,:)=Feat;
EDV(II)=t(Frames/dy+1).*dx;
out.Vols=[ EDV(II) ESV(II)]; % [end diastolic volume , end systolic volume]
out.Features=X(II,:);
out.FeatNames=FeatNames;
out.VolCurve=t(Frames/dy+1:Frames/dy*2).*dx;
out.Time=[1:numel(out.VolCurve)]./numel(out.VolCurve);
end
function [B,Rat]=prep(A);
ss=128;
A=squeeze(A);
S=size(A);
if S(1)<S(2)
    m=round((S(2)-S(1))/2); x1=1; x2=S(1); y1=m+1; y2=m+S(1);
end
if S(1)>=S(2)
    m=round((S(1)-S(2))/2); x1=m+1; x2=m+S(2); y1=1; y2=S(2);
end
A=A(x1:x2,y1:y2);
S=size(A);
m=round(S(1)*.1);
A=A(m+1:end-m,m+1:end-m);
A=flatout(A,.01);
S=size(A);
Rat=S(1)/ss;
B=imresize(A, [ss,ss]);
end
function [Feat,FeatNames]=Features(AAG,Res)
A=AAG;
Frames=size(A,4);
SLs=size(A,1);
F=zeros(9,Frames,SLs);
for I=1:Frames
    for J=1:SLs
        B=squeeze(A(J,:,:,I));
        F(1,I,J)=sum_(B==0)*Res(1)*Res(2);
        F(2,I,J)=sum_(B==1)*Res(1)*Res(2);
        F(3,I,J)=sum_(B==2)*Res(1)*Res(2);
        F(4,I,J)=sum_(bwdist(B~=0))*Res(1)*Res(2);
        F(5,I,J)=sum_(bwdist(B~=1))*Res(1)*Res(2);
        F(6,I,J)=sum_(bwdist(B~=2))*Res(1)*Res(2);
        F(7,I,J)=max_(bwdist(B~=0))*Res(1);
        F(8,I,J)=max_(bwdist(B~=1))*Res(1);
        F(9,I,J)=max_(bwdist(B~=2))*Res(1);
    end
end
F(find(isnan(F)))=0;
t1=mean(mean(F,2),3);
t2=mean(std(F,[],2),3);
t3=std(mean(F,2),[],3);
t4=std(std(F,[],2),[],3);
Types={'Volumes','Distances','Thicknesses'};
Loc={'RV','MY','LV'};
Method={'Mean of Mean of','Mean of Std of','Std of Mean of','Std of Std of'};
a=1;
for I=1:numel(Method)
    for J=1:numel(Types)
        for K=1:numel(Loc)
            FeatNames{a}=[Method{I} ' ' Loc{K} ' ' Types{J}];
            a=a+1;
        end
    end
end
Feat=[t1; t2 ;t3; t4];
end
function m=max_(A)
m=max(A(:));
end
function m=sum_(A)
m=sum(A(:));
end
function [A]=flatout(A,prc)
low=quantilex(A,prc);
high=quantilex(A,1-prc);
A(A<low)=low;
A(A>high)=high;
end
function [q]=quantilex(x,p)
x=sort(x(:)); id=ceil(p.*max(size(x)));
id(id==0)=1;
q=x(id);
end
