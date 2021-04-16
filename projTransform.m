% Circular Fisheye Image projection transformation
% currently support equidistant and equisolid
% Hankun Li, University of Kansas
% 04/16/2021

% Reference: 
% [1] 360-degree-image-processing (https://github.com/k-machida/360-degree-image-processing), GitHub. Retrieved August 18, 2020.
% [2] Tuan Ho, Madhukar Budagavi,  "2DUAL-FISHEYE LENS STITCHING FOR 360-DEGREE IMAGING"
% Copyright of Original Functions: Kazuya Machida (2020)

%%
clear all; %#ok<*CLALL>
close all;
clc;
sx = 22.3; sy = 14.9; f = 4.5; Lfov = 180; % camera parameters
% sensor size, focal length, lens FOV
%%
transform1 = 0; % 0: equisolid to equidrectangular 1: equidistant to equidrectangular...
transform2 = 1; % 0: equidrectangular to equisolid 1: equidrectangular to equidistant...
%
%% 1.load image
[fn,pn]=uigetfile('*.jpg','load 180 FOV image.');
str=[pn,fn]; imF = imread(str); clear pn fn str;
[y,x] = size(imF(:,:,1));

%% crop[optional], projection function is required.
pr180 = round(sind(45)*2*f/sy*y/2)*2; %projection: equisolid-angle
imF = circularCrop(imF,x,y,pr180);

%% 2. Projection Transformation
imE = imfish2equ(imF,transform1,Lfov);
imF_NEW = imequ2fish(imE,transform2,Lfov);

%% Show difference of transformation
K = imabsdiff(imF,imF_NEW);
figure(2); imshow(K);



%% functions

function [xf,yf] = equ2fish(xe,ye,fov,roll,tilt,pan,tkey)
thetaE = xe*180; phiE = ye*90; cosdphiE = cosd(phiE); 
xs = cosdphiE.*cosd(thetaE); ys = cosdphiE.*sind(thetaE); zs = sind(phiE);   
xyzsz = size(xs); xyz = xyzrotate([xs(:),ys(:),zs(:)],[roll tilt pan]);
xs = reshape(xyz(:,1),xyzsz(1),[]); 
ys = reshape(xyz(:,2),xyzsz(1),[]);
zs = reshape(xyz(:,3),xyzsz(1),[]);
thetaF = atan2d(zs,ys);
if tkey
    r = 2*atan2d(sqrt(ys.^2+zs.^2),xs)/fov; % equidistant
else
    r = 2*(sind(atan2d(sqrt(ys.^2+zs.^2),xs)/2))/(2*sind(fov/4)); % equisolid-angle
end
xf = r.*cosd(thetaF); yf = r.*sind(thetaF);
end

function [xe,ye] = fish2equ(xf,yf,roll,tilt,pan,fov,tkey)
thetaS = atan2d(yf,xf);
if tkey
    phiS = sqrt(yf.^2+xf.^2)*fov/2; % equidistant
else
    phiS = 2*asind(sqrt(yf.^2+xf.^2)/(2*sind(0.5*fov/2))); % equisolid-angle
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
% Get index of valid fisyeye image area
idx = sqrt(xf.^2+yf.^2) <= 1; xf = xf(idx); yf = yf(idx);
[xe,ye] = fish2equ(xf,yf,roll,tilt,pan,fov,tkey);
Xe = round((xe+1)/2*(we-1)+1); % rescale to 1~we
Ye = round((ye+1)/2*(he-1)+1); % rescale to 1~he
Xf = round((xf+1)/2*(wf-1)+1); % rescale to 1~wf
Yf = round((yf+1)/2*(hf-1)+1); % rescale to 1~hf
Ie = reshape(imgE,[],ch); If = zeros(hf*wf,ch,'uint8');
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
%fisheye image size
wf = size(imgF,2); hf = size(imgF,1); ch = size(imgF,3);
%equirectangular image size
we = wf*2; he = hf;
fov  = p.Results.fov; roll = p.Results.roll;
tilt = p.Results.tilt; pan  = p.Results.pan;
[xe,ye] = meshgrid(1:we,1:he);
xe = 2*((xe-1)/(we-1)-0.5); ye = 2*((ye-1)/(he-1)-0.5); 
[xf,yf] = equ2fish(xe,ye,fov,roll,tilt,pan,tkey); idx = sqrt(xf.^2+yf.^2) <=1; 
xf = xf(idx); yf = yf(idx); xe = xe(idx); ye = ye(idx);
Xe = round((xe+1)/2*(we-1)+1); Ye = round((ye+1)/2*(he-1)+1); 
Xf = round((xf+1)/2*(wf-1)+1); Yf = round((yf+1)/2*(hf-1)+1); 
Ie = reshape(imgF,[],ch); If = zeros(he*we,ch,'uint8');
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

function Icropped = circularCrop(I,x,y,r)
xc = round(x/2);yc = round(y/2);
c = zeros(y,x); [L(:,1),L(:,2)] = find(c==0);
L(:,3) = sqrt((L(:,1) - yc).^2 + (L(:,2) - xc).^2);
L(L(:, 3) > r, :) = [];
for i = 1: size(L,1)
   c(y+1-L(i,1),L(i,2)) = 1;
end
msk = imbinarize(c,0);
ir = uint8(double(I(:,:,1)).*msk);
ig = uint8(double(I(:,:,2)).*msk);
ib = uint8(double(I(:,:,3)).*msk);
Icc = cat(3,ir,ig,ib);
[mski(:,1), mski(:,2)] = find(msk==1);
Icropped = imcrop(Icc,[min(mski(:,2)),min(mski(:,1)),...
    max(mski(:,2))-min(mski(:,2)),max(mski(:,1))-min(mski(:,1))]);
end
