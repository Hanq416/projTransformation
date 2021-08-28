% Circular Fisheye Image projection transformation
% Hankun Li, University of Kansas
% 08/16/2021
% Reference: 
% [1] https://github.com/k-machida/360-degree-image-processing), GitHub.
clear all; %#ok<*CLALL>
close all;
clc;

%% Camera parameters: Canon T2i + sigma 4.5-f2.8, change to your own camera
sx = 22.3; sy = 14.9; f = 4.5; Lfov = 180; % camera parameters
% sensor size, focal length, lens FOV

%% circular fisheye projection transformation code:
% 0: Equisolid-angle
% 1: Equidistant
% 2: Orthographic
% 3: Stereographic
trans_from = 0; % source proj code
trans_to = 2; % target proj code

%% 1. load HDR, pre_processing
[fn,pn]=uigetfile('*.HDR','load 180 FOV hdr image.');
str=[pn,fn]; imF = hdrread(str); cf1 = ui_cf();
exp = getExpValue(str); cf = cf1./exp;
[y,x] = size(imF(:,:,1));
clear pn fn str;

%% pre_processing: crop, vc_correction
imF = fov180_crop(imF,sy,x,y,f,trans_from);
[~, imF] = VC_correction(imF,Lfov); % VC correction

%% display source HDR image [optional]
gamma = 0.5; % chnage gamma to adjust viewing brightness
[hdr_src] = hdrGammaShow(imF,0.5); 
figure(10); imshow(hdr_src);

%% Projection Transformation
imF_NEW = projTrans_hdr(imF,trans_from,trans_to,Lfov);

%% show luminance map, avg luminance
lmap = luminanceMap_gen(imF_NEW,cf);
Luminance_map_show(lmap,'Luminance map');
fprintf('Average luminance: %.2f cd/m2 \n',mean(nonzeros(lmap)));


%% Show difference of transformation, [optional]
close all;
K = imabsdiff(imF,imF_NEW);
figure(3); imshow(K); 
title('Difference between original/transformed image');












%end here

%% functions lib
function imF_NEW = projTrans_hdr(imF,trans_from,trans_to,Lfov)
imE = imfish2equ(imF,trans_from,Lfov);
imF_NEW = imequ2fish(imE,trans_to,Lfov);
[hdr_tar] = hdrGammaShow(imF_NEW,0.5);
f = figure(1); imshow(hdr_tar); uiwait(f);
end

function Luminance_map_show(lmap,name)
cv = std(std(lmap))/mean(mean(lmap));
lmap(lmap<0) = 0;lumimg = (lmap - min(min(lmap)))/(max(max(lmap))-min(min(lmap)));
if  (1.5<cv)&&(cv<10)
    gm = round(1/cv,2);
elseif cv>10
    gm = 0.09;
else
    gm = 1;
end
lumimg = uint8((lumimg.^gm).*256);
rg = max(max(lmap))-min(min(lmap)); crange = jet(256);crange(1,:) = 0;
cb1 = round(rg.*(0.03316.^(1/gm)),7);cb2 = round(rg.*(0.26754.^(1/gm)),2);
cb3 = round(rg.*(0.50191.^(1/gm)),2);cb4 = round(rg.*(0.73629.^(1/gm)),2);
cb5 = round(rg.*(0.97066.^(1/gm)),2);
figure(2);imshow(lumimg,'Colormap',crange);title(name);
hcb = colorbar('Ticks',[8,68,128,188,248],'TickLabels',{cb1,cb2,cb3,cb4,cb5});
title(hcb,'luminance (cd/m2)');
end

function lmap = luminanceMap_gen(hdr,cf)
lmap = (hdr(:,:,1).*0.265 + hdr(:,:,2).*0.670 + hdr(:,:,3).*0.065).*179.*cf;
end

function [vcmask, vc_hdr] = VC_correction(hdr,fov)
[yc,xc] = size(hdr(:,:,1));
[xf,yf] = meshgrid(1:xc,1:yc); 
xf = (xf - round(xc/2))./round(xc/2);
yf = (yf - round(yc/2))./round(yc/2);
phiS = 2*asind(sqrt(yf.^2+xf.^2)*sind(fov/4));
vcmask = canonVC(phiS); vcmask(vcmask<0) = 0;
bound = imbinarize(rgb2gray(hdr),0); 
vcmask = vcmask.*bound;
vc_hdr = hdr.*vcmask;
end

% VC correction, require calibration for every Camera + Lens
function [vcf] = canonVC(angle) % vc correction for Canon t2i + sigma f2.8
vcf = 1./(-4.3909e-09.*angle.^(4) - 3.9024e-07.*angle.^(3) +...
    3.3680e-05.*angle.^(2)-0.0018.*angle + 1.0018);
end

function cropped = fov180_crop(hdr,sy,x,y,f,projCode)
if projCode == 0
    pr180 = round(sind(45)*f/sy*y)*2; %projection: equisolid-angle
elseif projCode == 1
    pr180 = round(f/sy*y*0.25); %projection: equidistant
elseif projCode == 2
    pr180 = round(f/sy*y); %projection: orthogonal
elseif projCode == 3
    pr180 = round(tand(45)*f/sy*y)*2; %projection: equisolid-angle
else
    error("specified wrong input projection type!\n")
end
cropped = circularCrop_hdr(hdr,x,y,pr180);
end

function Icropped = circularCrop_hdr(I,x,y,r)
xc = round(x/2);yc = round(y/2);
c = zeros(y,x); [L(:,1),L(:,2)] = find(c == 0);
L(:,3) = sqrt((L(:,1) - yc).^2 + (L(:,2) - xc).^2);
L(L(:, 3) > r, :) = [];
for i = 1: size(L,1)
   c(y+1-L(i,1),L(i,2)) = 1;
end
msk = imbinarize(c,0);
ir = double(I(:,:,1)).*msk;
ig = double(I(:,:,2)).*msk;
ib = double(I(:,:,3)).*msk;
Icc = cat(3,ir,ig,ib);
[mski(:,1), mski(:,2)] = find(msk == 1);
Icropped = imcrop(Icc,[min(mski(:,2)),min(mski(:,1)),...
    max(mski(:,2))-min(mski(:,2)),max(mski(:,1))-min(mski(:,1))]);
end

function [ans1] = ui_cf()
prompt = {'Global luminance calibration factor?'};
dlgtitle = 'User Input'; dims = [1 50];definput = {'1.0'};
answer = str2double(inputdlg(prompt,dlgtitle,dims,definput));
if isempty(answer)
    ans1 = 1.0;
else
    ans1 = answer(1);
end
end

function exposure = getExpValue(filename)
fid = fopen(filename); ct = 0;
while ct < 16
    line = fgetl(fid);
    if contains(line, 'EXPOSURE')
        line = erase(line, ' ');
        break
    end
end
fclose(fid); exposure = str2double(erase(line, 'EXPOSURE='));
end

function fig1 = hdrGammaShow(imHDR,gamma)
fr = imHDR(:,:,1);fg = imHDR(:,:,2);fb = imHDR(:,:,3);
fr = single(fr).^gamma;fg = single(fg).^gamma;fb = single(fb).^gamma;
fig1 = cat(3,fr,fg,fb);
end

function [xf,yf] = equ2fish(xe,ye,fov,roll,tilt,pan,tkey)
thetaE = xe*180; phiE = ye*90; cosdphiE = cosd(phiE); 
xs = cosdphiE.*cosd(thetaE); ys = cosdphiE.*sind(thetaE); zs = sind(phiE);   
xyzsz = size(xs); xyz = xyzrotate([xs(:),ys(:),zs(:)],[roll tilt pan]);
xs = reshape(xyz(:,1),xyzsz(1),[]); 
ys = reshape(xyz(:,2),xyzsz(1),[]);
zs = reshape(xyz(:,3),xyzsz(1),[]);
thetaF = atan2d(zs,ys);
if tkey == 1
    r = 2*atan2d(sqrt(ys.^2+zs.^2),xs)/fov; % equidistant
elseif tkey == 0
    r = sind(atan2d(sqrt(ys.^2+zs.^2),xs)/2)/sind(fov/4); % equisolid-angle
elseif tkey == 2
    r = sind(atan2d(sqrt(ys.^2+zs.^2),xs))/sind(fov/2); % orthographic proj*
elseif tkey == 3
    r = tand(atan2d(sqrt(ys.^2+zs.^2),xs)/2)/tand(fov/4); % Stereographic*
end
xf = r.*cosd(thetaF); yf = r.*sind(thetaF);
end

function [xe,ye] = fish2equ(xf,yf,roll,tilt,pan,fov,tkey)
thetaS = atan2d(yf,xf);
if tkey == 1
    phiS = sqrt(yf.^2+xf.^2)*fov/2; % equidistant
elseif tkey == 0
    phiS = 2*asind(sqrt(yf.^2+xf.^2)*sind(fov/4)); % equisolidangle proj
elseif tkey == 2
    phiS = asind(sqrt(yf.^2+xf.^2)*sind(fov/2)); % orthographic proj*
elseif tkey == 3
    phiS = 2*atand(sqrt(yf.^2+xf.^2)*tand(fov/4)); % Stereographic proj*
end
sindphiS = sind(phiS);
xs = sindphiS.*cosd(thetaS); ys = sindphiS.*sind(thetaS); zs = cosd(phiS);
xyzsz = size(xs);
xyz = xyzrotate([xs(:),ys(:),zs(:)],[roll tilt pan]);
xs = reshape(xyz(:,1),xyzsz(1),[]); 
ys = reshape(xyz(:,2),xyzsz(1),[]);
zs = reshape(xyz(:,3),xyzsz(1),[]);
thetaE = atan2d(xs,zs); phiE = atan2d(ys,sqrt(xs.^2+zs.^2));
xe = thetaE/180; ye = 2*phiE/180;
end

function imgF = imequ2fish(imgE,tkey,varargin)
p = inputParser;
addRequired(p,'imgE'); addRequired(p,'tkey');
addOptional(p,'fov' ,  180); % defaul value of fov
addOptional(p,'roll',  0); % defaul value of roll
addOptional(p,'tilt',  0); % defaul value of tilt
addOptional(p,'pan' ,  0); % defaul value of pan
parse(p,imgE,tkey,varargin{:});
we = size(imgE,2); he = size(imgE,1); ch = size(imgE,3);
wf = round(we/2); hf = he;
roll = p.Results.roll; tilt = p.Results.tilt; 
pan  = p.Results.pan; fov = p.Results.fov;
[xf,yf] = meshgrid(1:wf,1:hf);
xf = 2*((xf-1)/(wf-1)-0.5); yf = 2*((yf-1)/(hf-1)-0.5); 
idx = sqrt(xf.^2+yf.^2) <= 1; xf = xf(idx); yf = yf(idx);
[xe,ye] = fish2equ(xf,yf,roll,tilt,pan,fov,tkey);
Xe = round((xe+1)/2*(we-1)+1); % rescale to 1~we
Ye = round((ye+1)/2*(he-1)+1); % rescale to 1~he
Xf = round((xf+1)/2*(wf-1)+1); % rescale to 1~wf
Yf = round((yf+1)/2*(hf-1)+1); % rescale to 1~hf
Ie = reshape(imgE,[],ch); If = zeros(hf*wf,ch,'double');
idnf = sub2ind([hf,wf],Yf,Xf);idne = sub2ind([he,we],Ye,Xe);
If(idnf,:) = Ie(idne,:);imgF = reshape(If,hf,wf,3);
end

function imgE = imfish2equ(imgF,tkey,varargin)
p = inputParser;
addRequired(p,'imgF');
addRequired(p,'tkey');
addOptional(p,'fov' ,180); % defaul value of fov
addOptional(p,'roll',  0); % defaul value of roll
addOptional(p,'tilt',  0); % defaul value of tilt
addOptional(p,'pan' ,  0); % defaul value of pan
parse(p,imgF,tkey,varargin{:});
wf = size(imgF,2); hf = size(imgF,1); ch = size(imgF,3);
we = wf*2; he = hf;
fov  = p.Results.fov; roll = p.Results.roll;
tilt = p.Results.tilt; pan  = p.Results.pan;
[xe,ye] = meshgrid(1:we,1:he);
xe = 2*((xe-1)/(we-1)-0.5); ye = 2*((ye-1)/(he-1)-0.5); 
[xf,yf] = equ2fish(xe,ye,fov,roll,tilt,pan,tkey); idx = sqrt(xf.^2+yf.^2) <=1; 
xf = xf(idx); yf = yf(idx); xe = xe(idx); ye = ye(idx);
Xe = round((xe+1)/2*(we-1)+1); Ye = round((ye+1)/2*(he-1)+1); 
Xf = round((xf+1)/2*(wf-1)+1); Yf = round((yf+1)/2*(hf-1)+1); 
Ie = reshape(imgF,[],ch); If = zeros(he*we,ch,'double');
idnf = sub2ind([hf,wf],Yf,Xf); idne = sub2ind([he,we],Ye,Xe);
If(idne,:) = Ie(idnf,:);imgE = reshape(If,he,we,3);
end

function [xyznew] = xyzrotate(xyz,thetaXYZ)
tX =  thetaXYZ(1); tY =  thetaXYZ(2); tZ =  thetaXYZ(3);
Tm = [ cosd(tY)*cosd(tZ),- cosd(tY)*sind(tZ), sind(tY); ...
      cosd(tX)*sind(tZ) + cosd(tZ)*sind(tX)*sind(tY), cosd(tX)*cosd(tZ) - sind(tX)*sind(tY)*sind(tZ), -cosd(tY)*sind(tX); ...
      sind(tX)*sind(tZ) - cosd(tX)*cosd(tZ)*sind(tY), cosd(tZ)*sind(tX) + cosd(tX)*sind(tY)*sind(tZ),  cosd(tX)*cosd(tY)];
xyznew = xyz*Tm;
end
