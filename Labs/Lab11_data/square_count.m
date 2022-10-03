% Yield Calculator
% Written by Philip Petersen (pfp@caltech.edu)
% Qian Lab, Caltech 
% -----------------------------------------------------------------
% left click on the image to select; right click on the image to open
% the menu in the Matlab command window; middle click to exit
% Yield = selected_pixels / all_pixels_over_threshold

clear all % funcationalize if desired
I=imread('Inkedsquare_tile_004_matlab.jpg'); % input image

LB = 20; % Lower bound on size in pixels; removes random pixel bumps in image (do not set to a value larger than the smallest tile fragment desired in the calculation)
alpha = .55; % image and color mixing level
tile_color = [0 .8 1]; %blue background tiles
select_color = [1 1 0]; %yellow selected tiles

imshow(I)
I2=rgb2gray(I);
% I2 = imadjust(I2,[0 .1]); % in case you are inputting a really bad image
% I2=medfilt2(I2,[10 1]); % removes horizontal AFM streaking noise if neccesary, comment out this line otherwise

level = graythresh(I2); % automatic threshold for well flattened images
%level = .12; % manual threshold level if automatic fails


%--------- Begin code-------------

bw = im2bw(I2,level); % apply threshold
bw = bwareaopen(bw, LB); % remove all background spikes

cc = bwconncomp(bw); % start making a list of tile identities and color them
labeled = labelmatrix(cc);
colorm=repmat(tile_color,[cc.NumObjects 1]);
RGB_label = label2rgb(labeled, colorm, [.7 .7 .7], 'noshuffle');

 imshow(I); hold on;
 h=imagesc(RGB_label);
 set( h, 'AlphaData', alpha ); % create blended image

  s = regionprops(labeled,'PixelList');
  
  temp=zeros(1,size(s,1));
  tempo= zeros(1,size(s,1));
  z=1;
  sel=2;
while(z~=2)
[x,y z] = mginput(1);

if (z==3)
    sel = input('Input: \n0-exit;\n1-cut tool (click vertices to remove everything in the polygon mask; double click to confirm mask);\n2-click to select;\n3-click to cancel selection;\n4-click to delete object;\n5-separate touching tiles (use for touching objects; create a polygon around object to separate)\n6-bonus tools: ');
end

if (sel==0) z=2; end

if (sel==2) %click to add to count; default after cutting
    if (z==1)
        ID = labeled(round(y),round(x));
    for i=1:size(s,1)
        temp(i)=any(ismember(round([x y]),s(i).PixelList,'rows'));
    end
    temp=temp|tempo;
    colorm(find(temp),:)=repmat(select_color,[sum(temp), 1]);
    end
end

if (sel==1) %cut areas out manually
         cut_mask = roipoly; %select region to remove
         IDs = unique(labeled(cut_mask));IDs=IDs(IDs~=0); %objects cut
         temp(IDs)=0; %clear any object cut
         fbw = ismember(labeled,IDs); %ROI on area manipulated
         bw = bw&(~cut_mask); %new masks 
         fbw = fbw&(~cut_mask); %ROI new masked
         cc = bwconncomp(fbw);
         labeled2 = labelmatrix(cc);
         labeled = (double(labeled).*double(~fbw).*double(~cut_mask))+(double(labeled2)+double((labeled2>0)*length(temp)));
         temp=[temp zeros(1,cc.NumObjects)]; %grow list of objects
         colorm=repmat(tile_color,[length(temp) 1]);
         colorm(find(temp),:)=repmat(select_color,[sum(temp), 1]);
         s = regionprops(labeled,'PixelList');
         sel=2; %default to add tiles on click
end

if (sel==3) %click to cancel selection
    for i=1:size(s,1)
        temp(i)=any(ismember(round([x y]),s(i).PixelList,'rows'));
    end
    temp=(~temp)&tempo;
    colorm(find(temp),:)=repmat(select_color,[sum(temp), 1]);
    colorm(find(~temp),:)=repmat(tile_color,[sum(~temp), 1]);
end

if (sel==4) %click to delete object (removes blob from calculations completely)
    for i=1:size(s,1)
        temp(i)=any(ismember(round([x y]),s(i).PixelList,'rows'));
    end
    ID = labeled(round(y),round(x));
    temp=(~temp)&tempo; %clear data from calcs;
    labeled(labeled==ID)=0; 
    bw(labeled==ID)=0;
end

if (sel==5) %separate areas out manually
         cut_mask = roipoly; %select region to make unique
         labeled(cut_mask) = (labeled(cut_mask)>0)*(length(temp)+1);
         temp=[temp 0]; %grow like of objects
         colorm=repmat(tile_color,[length(temp) 1]);
         colorm(find(temp),:)=repmat(select_color,[sum(temp), 1]);
         s = regionprops(labeled,'PixelList');
         sel=2; %default to add tiles on click
end

if (sel==6)
    sel2 = input(['Bonus: \n0-exit menu\n1-delete all objects touching image edges:\n2-delete all object smaller than X% of average selected object (meant to remove peppercorn noise)\n' ...
    '3-help me find similar objects (just highlighted-not part of calculation; using a tool like cut [1] will remove the highlights)\n4-select all objects (when your yield is high; then cancel bad objects with [3]): ']);
    if (sel2==1)
        labeled=imclearborder(labeled);
        bw=imclearborder(bw);
    end
    if (sel2==2)
        sel3 = input('Input deletion fraction (ie, .05) (this uses the average of selected tiles):');
        R = regionprops(labeled,'Area');
        avg_obj_size = mean([R(find(temp)).Area]);
        labeled=double(bwareaopen(labeled,round(avg_obj_size*sel3))).*double(labeled);
        bw=bwareaopen(bw,round(avg_obj_size*sel3));
    end
    if (sel2==3)
        R = regionprops(labeled,'Area');
        R1 = regionprops(labeled,'Eccentricity');
        avg_obj_size = mean([R(find(temp)).Area])
        std_area = std([R(find(temp)).Area])
        avg_obj_ecc = mean([R1(find(temp)).Eccentricity])
        std_ecc = std([R1(find(temp)).Eccentricity])
        t1=[R1.Eccentricity];
        t2=[R.Area];
        l=(t2<(avg_obj_size+4*std_area)) & (t2>(avg_obj_size-4*std_area)) & (t1<(avg_obj_ecc+3*std_ecc)); %similar objects eccentricity and area
        l(temp(1:length(l)))=0; %don't overlay already selected tiles
        colorm(find(l),:)=repmat([1 0 0],[sum(l), 1]);
    end
    if (sel2==4)
        temp(unique(labeled(labeled>0)))=1; %set all to true
        colorm=repmat(select_color,[length(temp) 1]);
    end
    sel=2;
end



 RGB_label = label2rgb(labeled, colorm, [.7 .7 .7], 'noshuffle');
 hold off
 imshow(I); hold on;
 h=imagesc(RGB_label);
 set( h, 'AlphaData', alpha );
 
 R = regionprops(labeled,'Area');

 tempo=temp;

 %------- Some Output Math Calculation --------------------
 Total_Area = prod(size(I2))
 Area_tiles = sum([R.Area]) %area with tile
 Area_selected = sum([R(find(temp(1:length(R)))).Area]) %area of selected
 ratio_squares = Area_selected / Area_tiles % ratio

end