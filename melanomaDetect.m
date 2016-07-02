
% Skin Cancer Risk Quantification using ABCD Rule
% Aditya D. Pai
% Influenced and built upon the work by John Breneman

clear all
% Load image
file='FriendX.jpg';
img=imread(file);

% Choose Outputs
write_image=0;
write_file=0;

%Preprocessing
img=imresize(img,512/size(img,1));
imggray=rgb2gray(img);

%% Analysis

% Binarize 
imgbin=~(im2bw(imggray,graythresh(imggray))); % using Otsu's method

%% Small region removal

imgSRR=imgbin;

posthr=500; 
negthr=1500; 
label1=bwlabel(imgSRR,4);
binarys=1:max(label1(:));
ahist=histc(label1(:),binarys);
binlist=binarys(ahist<posthr);
imgSRR(ismember(label1,binlist))=0;

label1=bwlabel(~imgSRR,4);
binarys=1:max(label1(:));
ahist=histc(label1(:),binarys);
binlist=binarys(ahist<negthr);
imgSRR(ismember(label1,binlist))=1;

%% Identify primary region of interest

imgx0=size(imggray,2)/2;
imgy0=size(imggray,1)/2;
label1=bwlabel(imgSRR,4);
cent=regionprops(label1,'Centroid');
mindistance=inf;
for rgn=1:numel(cent)
    distance=sqrt(((cent(rgn).Centroid(1)-imgx0)^2) + ...
        ((cent(rgn).Centroid(2)-imgy0)^2));
    if distance<mindistance
        mindistance=distance;
        minrgn=rgn;
    end
end

primereg=imgSRR;
primereg(label1~=minrgn)=0;

%% Boundary detection

regionedge=edge(primereg);
[edgex, edgey]=ind2sub(size(regionedge),find(regionedge));

%% Analyze for symmetry

index=find(primereg);
[regionx, regiony]=ind2sub(size(primereg),index);
xrotate=round(2*cent(minrgn).Centroid(2)-regionx);
yrotate=round(2*cent(minrgn).Centroid(1)-regiony);
indexrotate=sub2ind(size(primereg),xrotate,yrotate);

% Calculate symmetry MSE
A1=mean((imgSRR(index)-imgSRR(indexrotate)).^2);

%% Calculate image intensity variance within region
C1=std(double(imggray(primereg)));

%% Analyze border strength


G=fspecial('gaussian',30,2); % Calculate gradient magnitude
[dx, dy]=gradient(G);

deltax=imfilter(double(imggray),dx,'symmetric');
deltay=imfilter(double(imggray),dy,'symmetric');

bordergrad=sqrt(deltax.*deltax + deltay.*deltay);

indexedge=find(regionedge);
B1=mean(bordergrad(indexedge));

%% Diameter Calculation

Img1 = rgb2gray(img);

strucel1 = strel('disk', 20);
Ierode = imerode(Img1, strucel1);
Irecon = imreconstruct(Ierode, Img1);
Idil = imdilate(Irecon, strucel1);
Icomp = imreconstruct(imcomplement(Idil), imcomplement(Irecon));
Icomp = imcomplement(Icomp);
Icompbw = im2bw(Icomp, graythresh(Icomp));
bwcomp=imcomplement(Icompbw);

major=regionprops(bwcomp,'majoraxislength');
minor=regionprops(bwcomp,'MinorAxisLength');
max=major.MajorAxisLength;
min=minor.MinorAxisLength;
if(max>=min)
    D1=min;
end

%% Normalization and Weighted Score Calculation
A=(A1-0.023511)/0.173179;
if(A<0)
    A=0;
elseif(A>1)
    A=1;
end
B=(B1-2.4553)/7.0059;
if(B<0)
    B=0;
elseif(B>1)
    B=1;
end
C=(C1-11.9032)/29.4863;
if(C<0)
    C=0;
elseif(C>1)
    C=1;
end
D=(D1-96.2135)/261.054;
if(D<0)
    D=0;
elseif(D>1)
    D=1;
end
Total=2.4*A+2.4*B+3*C+2*D

%% OUTPUT GENERATION

f=figure(2);
imshow(img)
hold on
plot(edgey,edgex,'b.')

ac=colorbar;
colormap jet;
caxis([0,9.8]);

if (Total>=4.4 & Total<=4.6)
    ac.TicksMode='manual';
    ac.Ticks=[0,Total,9.8];
    ac.TickLabelsMode='manual';
    ac.TickLabels={'Very Low Risk','Your Risk Score','Very High Risk'};

elseif (Total>4.5)
    ac.TicksMode='manual';
    ac.Ticks=[0,4.5,Total,9.8];
    ac.TickLabelsMode='manual';
    ac.TickLabels={'Very Low Risk','Threshold Point','Your Risk Score','Very High Risk'};
    
else
    ac.TicksMode='manual';
    ac.Ticks=[0,Total,4.5,9.8];
    ac.TickLabelsMode='manual';
    ac.TickLabels={'Very Low Risk','Your Risk Score','Threshold Point','Very High Risk'};
end

ac.FontWeight='bold';

hold off
title('Your Risk Report', 'FontSize',16)

dat={['Parameter values:'],...
    ['Asymmetry MSE: ' num2str(A1)],...
    ['Border Strength: ' num2str(B1)],...
    ['Color Deviation: ' num2str(C1)],...
    ['Diameter: ' num2str(D1)]}

msg1={['Total Score: ' num2str(Total) ' out of 9.8'],...
      ['Precautionary threshold point breached. IMMEDIATE further action warranted.']};
msg2={['Total Score: ' num2str(Total) ' out of 9.8'],...
      ['Nevus is likely benign. No further action warranted.']};

%ylabel(dat,'FontSize',14,'FontWeight','bold')

if Total>=4.5
   xlabel(msg1,'FontSize',14,'FontWeight','bold','Color','r')
else
   xlabel(msg2,'FontSize',14,'FontWeight','bold')
end

[~, name, ~]=fileparts(file);

if write_image
    print(f,'-djpeg','-r400',[name '_Report.jpg'])
end

formatSpec = '%s\n';

[n1rows,n1cols] = size(dat);
[n2rows,n2cols] = size(msg1);

if write_file
    fid = fopen( [name '_Results.txt'], 'wt' );
    fprintf(fid, '***************  YOUR RESULTS  ***************\n\nName:\t%s\n\n',name);
    for row = 1:n1rows
            fprintf(fid,formatSpec,dat{row,:});
    end
    fprintf(fid, '\n\n');
    if Total>=4.5
        for row1 = 1:n2rows
            fprintf(fid,formatSpec,msg1{row1,:});
        end
    else
        for row1 = 1:n2rows
            fprintf(fid,formatSpec,msg2{row1,:});
        end
    end
    fclose(fid);
end