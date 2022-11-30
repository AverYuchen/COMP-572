%Load images_1:Brooklyn Bridge
Img = imread('https://user-images.githubusercontent.com/83625350/204363127-1000fabc-aeec-4426-bfe3-87712962cc2d.jpeg');
Img = imresize(Img, 0.25);
%Find energy of image
energyimage = energy(Img);

vertical_minEnergy = vert_minEnergymap(energyimage);
horizontal_minEnergy = hori_minEnergymap(energyimage);

vertical_Seam = optimal_vertseam(vertical_minEnergy);
horizontal_Seam = optimal_horiseam(horizontal_minEnergy);
vertical_Seamline = seamMark(Img,vertical_Seam,'VERTICAL');
horizontal_Seamline = seamMark(Img,horizontal_Seam,'HORIZONTAL');

%seam_carving with horizontal seam line
final_img_ver = Img;
final_energyimg_ver = energyimage;
for i = 1:120
    [final_img_ver, final_energyimg_ver] = removeHeight(final_img_ver,final_energyimg_ver);
end

%seam_carving with vertical seam line
final_img_hor = Img;
final_energyimg_hor = energyimage;
for i = 1:150
    [final_img_hor, final_energyimg_hor] = removeWidth(final_img_hor,final_energyimg_hor);
end

%show result
figure(1);
subplot(2,2,1);imshow(Img);
title('Original image');
subplot(2,2,2);imshow(final_img_ver);
title('Result image: removed 120 height');
subplot(2,2,3);imshow(energyimage);
title('Original energy image');
subplot(2,2,4);imshow(final_energyimg_ver);
title('Result energy image');

figure(2);
subplot(1,2,1);imshow(horizontal_minEnergy);
title('Horizontal carved energy map');
subplot(1,2,2);imshow(horizontal_Seamline);
title('Horizontal carved seam line');

figure(3);subplot(2,2,1);imshow(Img);
title('Original image');
subplot(2,2,2);imshow(final_img_hor);
title('Result image: removed 150 width');
subplot(2,2,3);imshow(energyimage);
title('Original energy image');
subplot(2,2,4);imshow(final_energyimg_hor);
title('Result energy image');

figure(4);
subplot(1,2,1);imshow(vertical_minEnergy);
title('Vertical carved energy map');
subplot(1,2,2);imshow(vertical_Seamline);
title('Vertical carved seam line');

%Load images_2: Milan
Img = imread('https://user-images.githubusercontent.com/83625350/204571760-2fa3f102-9342-4e2e-adfd-061ef22aca1a.jpeg');
Img = imresize(Img, 0.25);
%Find energy of image
energyimage = energy(Img);

vertical_minEnergy = vert_minEnergymap(energyimage);
horizontal_minEnergy = hori_minEnergymap(energyimage);

vertical_Seam = optimal_vertseam(vertical_minEnergy);
horizontal_Seam = optimal_horiseam(horizontal_minEnergy);
vertical_Seamline = seamMark(Img,vertical_Seam,'VERTICAL');
horizontal_Seamline = seamMark(Img,horizontal_Seam,'HORIZONTAL');

%seam_carving with horizontal seam line
final_img_ver = Img;
final_energyimg_ver = energyimage;
for i = 1:120
    [final_img_ver, final_energyimg_ver] = removeHeight(final_img_ver,final_energyimg_ver);
end

%seam_carving with vertical seam line
final_img_hor = Img;
final_energyimg_hor = energyimage;
for i = 1:150
    [final_img_hor, final_energyimg_hor] = removeWidth(final_img_hor,final_energyimg_hor);
end

%show result
figure(5);
subplot(2,2,1);imshow(Img);
title('Original image');
subplot(2,2,2);imshow(final_img_ver);
title('Result image: removed 120 height');
subplot(2,2,3);imshow(energyimage);
title('Original energy image');
subplot(2,2,4);imshow(final_energyimg_ver);
title('Result energy image');

figure(6);
subplot(1,2,1);imshow(horizontal_minEnergy);
title('Horizontal carved energy map');
subplot(1,2,2);imshow(horizontal_Seamline);
title('Horizontal carved seam line');

figure(7);subplot(2,2,1);imshow(Img);
title('Original image');
subplot(2,2,2);imshow(final_img_hor);
title('Result image: removed 150 width');
subplot(2,2,3);imshow(energyimage);
title('Original energy image');
subplot(2,2,4);imshow(final_energyimg_hor);
title('Result energy image');

figure(8);
subplot(1,2,1);imshow(vertical_minEnergy);
title('Vertical carved energy map');
subplot(1,2,2);imshow(vertical_Seamline);
title('Vertical carved seam line');

%Load images_3: worldcup
Img = imread('https://user-images.githubusercontent.com/83625350/204574306-f5eb14e0-2542-4dee-b004-85d88c143686.jpeg');
Img = imresize(Img, 0.25);
%Find energy of image
energyimage = energy(Img);

vertical_minEnergy = vert_minEnergymap(energyimage);
horizontal_minEnergy = hori_minEnergymap(energyimage);

vertical_Seam = optimal_vertseam(vertical_minEnergy);
horizontal_Seam = optimal_horiseam(horizontal_minEnergy);
vertical_Seamline = seamMark(Img,vertical_Seam,'VERTICAL');
horizontal_Seamline = seamMark(Img,horizontal_Seam,'HORIZONTAL');

%seam_carving with horizontal seam line
final_img_ver = Img;
final_energyimg_ver = energyimage;
for i = 1:40
    [final_img_ver, final_energyimg_ver] = removeHeight(final_img_ver,final_energyimg_ver);
end

%seam_carving with vertical seam line
final_img_hor = Img;
final_energyimg_hor = energyimage;
for i = 1:150
    [final_img_hor, final_energyimg_hor] = removeWidth(final_img_hor,final_energyimg_hor);
end

%show result
figure(9);
subplot(2,2,1);imshow(Img);
title('Original image');
subplot(2,2,2);imshow(final_img_ver);
title('Result image: removed 40 height');
subplot(2,2,3);imshow(energyimage);
title('Original energy image');
subplot(2,2,4);imshow(final_energyimg_ver);
title('Result energy image');

figure(10);
subplot(1,2,1);imshow(horizontal_minEnergy);
title('Horizontal carved energy map');
subplot(1,2,2);imshow(horizontal_Seamline);
title('Horizontal carved seam line');

figure(11);subplot(2,2,1);imshow(Img);
title('Original image');
subplot(2,2,2);imshow(final_img_hor);
title('Result image: removed 150 width');
subplot(2,2,3);imshow(energyimage);
title('Original energy image');
subplot(2,2,4);imshow(final_energyimg_hor);
title('Result energy image');

figure(12);
subplot(1,2,1);imshow(vertical_minEnergy);
title('Vertical carved energy map');
subplot(1,2,2);imshow(vertical_Seamline);
title('Vertical carved seam line');

%Load images_4:Beijing
Img = imread('https://user-images.githubusercontent.com/83625350/204575512-e62c1c51-84f3-4ef6-b98f-428bcc20d66f.jpeg');
Img = imresize(Img, 0.25);
%Find energy of image
energyimage = energy(Img);

vertical_minEnergy = vert_minEnergymap(energyimage);
horizontal_minEnergy = hori_minEnergymap(energyimage);

vertical_Seam = optimal_vertseam(vertical_minEnergy);
horizontal_Seam = optimal_horiseam(horizontal_minEnergy);
vertical_Seamline = seamMark(Img,vertical_Seam,'VERTICAL');
horizontal_Seamline = seamMark(Img,horizontal_Seam,'HORIZONTAL');

%seam_carving with horizontal seam line
final_img_ver = Img;
final_energyimg_ver = energyimage;
for i = 1:40
    [final_img_ver, final_energyimg_ver] = removeHeight(final_img_ver,final_energyimg_ver);
end

%seam_carving with vertical seam line
final_img_hor = Img;
final_energyimg_hor = energyimage;
for i = 1:150
    [final_img_hor, final_energyimg_hor] = removeWidth(final_img_hor,final_energyimg_hor);
end

%show result
figure(13);
subplot(2,2,1);imshow(Img);
title('Original image');
subplot(2,2,2);imshow(final_img_ver);
title('Result image: removed 40 height');
subplot(2,2,3);imshow(energyimage);
title('Original energy image');
subplot(2,2,4);imshow(final_energyimg_ver);
title('Result energy image');

figure(14);
subplot(1,2,1);imshow(horizontal_minEnergy);
title('Horizontal carved energy map');
subplot(1,2,2);imshow(horizontal_Seamline);
title('Horizontal carved seam line');

figure(15);subplot(2,2,1);imshow(Img);
title('Original image');
subplot(2,2,2);imshow(final_img_hor);
title('Result image: removed 150 width');
subplot(2,2,3);imshow(energyimage);
title('Original energy image');
subplot(2,2,4);imshow(final_energyimg_hor);
title('Result energy image');

figure(16);
subplot(1,2,1);imshow(vertical_minEnergy);
title('Vertical carved energy map');
subplot(1,2,2);imshow(vertical_Seamline);
title('Vertical carved seam line');

%Load images_5: self image taking in art musuem

Img = imread('https://user-images.githubusercontent.com/83625350/204576450-e4242d07-d355-45ee-8bf4-67bee3a04bc0.jpg');
Img = imresize(Img, 0.25);
%Find energy of image
energyimage = energy(Img);

vertical_minEnergy = vert_minEnergymap(energyimage);
horizontal_minEnergy = hori_minEnergymap(energyimage);

vertical_Seam = optimal_vertseam(vertical_minEnergy);
horizontal_Seam = optimal_horiseam(horizontal_minEnergy);
vertical_Seamline = seamMark(Img,vertical_Seam,'VERTICAL');
horizontal_Seamline = seamMark(Img,horizontal_Seam,'HORIZONTAL');

%seam_carving with horizontal seam line
final_img_ver = Img;
final_energyimg_ver = energyimage;
for i = 1:40
    [final_img_ver, final_energyimg_ver] = removeHeight(final_img_ver,final_energyimg_ver);
end

%seam_carving with vertical seam line(fail case: image distortion)
final_img_hor = Img;
final_energyimg_hor = energyimage;
for i = 1:150
    [final_img_hor, final_energyimg_hor] = removeWidth(final_img_hor,final_energyimg_hor);
end

%seam_carving with vertical seam line
final_img_hor_2 = Img;
final_energyimg_hor_2 = energyimage;
for i = 1:15
    [final_img_hor_2, final_energyimg_hor_2] = removeWidth(final_img_hor_2,final_energyimg_hor_2);
end

%show result
figure(17);
subplot(2,2,1);imshow(Img);
title('Original image');
subplot(2,2,2);imshow(final_img_ver);
title('Result image: removed 40 height');
subplot(2,2,3);imshow(energyimage);
title('Original energy image');
subplot(2,2,4);imshow(final_energyimg_ver);
title('Result energy image');

figure(18);
subplot(1,2,1);imshow(horizontal_minEnergy);
title('Horizontal carved energy map');
subplot(1,2,2);imshow(horizontal_Seamline);
title('Horizontal carved seam line');

figure(19);subplot(2,2,1);imshow(Img);
title('Original image');
subplot(2,2,2);imshow(final_img_hor);
title('Result image: removed 150 width');
subplot(2,2,3);imshow(energyimage);
title('Original energy image');
subplot(2,2,4);imshow(final_energyimg_hor);
title('Result energy image');

figure(20);
subplot(1,2,1);imshow(vertical_minEnergy);
title('Vertical carved energy map');
subplot(1,2,2);imshow(vertical_Seamline);
title('Vertical carved seam line');

%show reduce width carving 
figure(21);subplot(2,2,1);imshow(Img);
title('Original image');
subplot(2,2,2);imshow(final_img_hor_2);
title('Result image: removed 15 width');
subplot(2,2,3);imshow(energyimage);
title('Original energy image');
subplot(2,2,4);imshow(final_energyimg_hor_2);
title('Result energy image');


%Load images_6: self image taking in MOMA

Img = imread('https://user-images.githubusercontent.com/83625350/204892296-29842de8-6f4e-4ee6-b98b-4bcd5872b743.jpg');
Img = imresize(Img, 0.25);
%Find energy of image
energyimage = energy(Img);

vertical_minEnergy = vert_minEnergymap(energyimage);
horizontal_minEnergy = hori_minEnergymap(energyimage);

vertical_Seam = optimal_vertseam(vertical_minEnergy);
horizontal_Seam = optimal_horiseam(horizontal_minEnergy);
vertical_Seamline = seamMark(Img,vertical_Seam,'VERTICAL');
horizontal_Seamline = seamMark(Img,horizontal_Seam,'HORIZONTAL');

%seam_carving with horizontal seam line
final_img_ver = Img;
final_energyimg_ver = energyimage;
for i = 1:70
    [final_img_ver, final_energyimg_ver] = removeHeight(final_img_ver,final_energyimg_ver);
end

%seam_carving with vertical seam line(fail case: image distortion)
final_img_hor = Img;
final_energyimg_hor = energyimage;
for i = 1:50
    [final_img_hor, final_energyimg_hor] = removeWidth(final_img_hor,final_energyimg_hor);
end

%show result
figure(22);
subplot(2,2,1);imshow(Img);
title('Original image');
subplot(2,2,2);imshow(final_img_ver);
title('Result image: removed 70 height');
subplot(2,2,3);imshow(energyimage);
title('Original energy image');
subplot(2,2,4);imshow(final_energyimg_ver);
title('Result energy image');

figure(23);
subplot(1,2,1);imshow(horizontal_minEnergy);
title('Horizontal carved energy map');
subplot(1,2,2);imshow(horizontal_Seamline);
title('Horizontal carved seam line');

figure(24);subplot(2,2,1);imshow(Img);
title('Original image');
subplot(2,2,2);imshow(final_img_hor);
title('Result image: removed 50 width');
subplot(2,2,3);imshow(energyimage);
title('Original energy image');
subplot(2,2,4);imshow(final_energyimg_hor);
title('Result energy image');

figure(25);
subplot(1,2,1);imshow(vertical_minEnergy);
title('Vertical carved energy map');
subplot(1,2,2);imshow(vertical_Seamline);
title('Vertical carved seam line');

%%%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Energy calculation(dual gradient energy function) %
%Energy of pixel (x, y) = x_gradient + y_gradient  %
%x_gradient = Rx (x, y)2 + Gx (x, y)2 + Bx (x, y)2 %
%y_gradient = Ry (x, y)2 + Gy (x, y)2 + By (x, y)2 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function energyImg = energy(img)
    result = im2double(im2gray(img));
    filterX = [-1,1];
    filterY = [-1;1];
    gradientX = abs(imfilter(result, filterX));
    gradientY = abs(imfilter(result,filterY));
    energyImg = sqrt(gradientX.^2+gradientY.^2);
end

function vert_minEnergy = vert_minEnergymap(energyImg)
    [rows, cols] = size(energyImg);
    %create a zeros matrix with dimension of energyImg matrix
    Map = zeros(rows,cols);
    %set first col = energyImg col
    Map(1,:) = energyImg(1,:);
    for i = 2:rows
        for j = 1: cols
            if j==1
                Map(i,j) = energyImg(i,j) + min([Map(i-1,j),Map(i-1,j+1)]);
            elseif j ==cols
                Map(i,j) = energyImg(i,cols) + min([Map(i-1,cols),Map(i-1,cols-1)]);
            else 
                Map(i,j) = energyImg(i,j) + min([Map(i-1,j-1),Map(i-1,j),Map(i-1,j+1)]);
            end
        end
    end
    vert_minEnergy = Map;
end

function hori_minEnergy = hori_minEnergymap(energyImg)
    [rows, cols] = size(energyImg);
    %create a zeros matrix with dimension of energyImg matrix
    Map = zeros(rows,cols);
    %set first row = energyImg row
    Map(:,1) = energyImg(:,1);
    for i = 2:cols
        for j = 1:rows
            if j==1
                Map(j,i) = energyImg(j,i) + min([Map(j,i-1),Map(j+1,i-1)]);
            elseif j ==rows
                Map(j,i) = energyImg(j,i) + min([Map(j-1,i-1),Map(j,i-1)]);
            else 
                Map(j,i) = energyImg(j,i) + min([Map(j-1,i-1),Map(j,i-1),Map(j+1,i-1)]);
            end
        end
    end
    hori_minEnergy = Map;
end

%optimal_vertseam helps us to find the lowest vertical energy
function vert_seam = optimal_vertseam(vert_minEnergy)
    [rows,cols] = size(vert_minEnergy);
    vert_seam = zeros(rows,1);
    [min_value,min_idx] = min(vert_minEnergy(rows,:));
    vert_seam(rows,1) = min_idx;
    
    for i = rows-1:-1:1
        %case1: at first col
        if min_idx == 1
            above = vert_minEnergy(i, min_idx);
            diagonal = vert_minEnergy(i, min_idx+1);
            [~,direction] = min([above, diagonal]);
            direction = direction+1;
        %case2: at last col
        elseif min_idx == cols
            above = vert_minEnergy(i, min_idx);
            diagonal = vert_minEnergy(i, min_idx-1);
            [~,direction] = min([above, diagonal]);
        %common case: inside boundaries
        else
            left = vert_minEnergy(i,vert_seam(i+1)-1);
            right = vert_minEnergy(i,vert_seam(i+1)+1);
            above = vert_minEnergy(i,vert_seam(i+1));
            [~,direction] = min([left, right, above]);
        end
        
        if direction == 3
            min_idx = vert_seam(i+1)+1;
        elseif direction == 2
            min_idx = vert_seam(i+1);
        else
            min_idx = vert_seam(i+1)-1;
        end
        vert_seam(i,1) = min_idx;
    end
end

function hori_seam = optimal_horiseam(hori_minEnergy)
    [rows,cols] = size(hori_minEnergy);
    hori_seam = zeros(1,cols);
    [~,min_idx] = min(hori_minEnergy(:,cols));
    hori_seam(1,cols) = min_idx;

    for i = cols-1:-1:1
        %case1: at first col
        if min_idx == 1
            right = hori_minEnergy(min_idx,i);
            diagonal = hori_minEnergy(min_idx+1,i);
            [val,direction] = min([right, diagonal]);
            direction = direction+1;
        %case2: at last col
        elseif min_idx == rows
            right = hori_minEnergy(min_idx-1,i);
            diagonal = hori_minEnergy(min_idx,i);
            [val,direction] = min([right, diagonal]);
        %common case: inside boundaries
        else
            above = hori_minEnergy(hori_seam(i+1)-1,i);
            bottom = hori_minEnergy(hori_seam(i+1)+1,i);
            right = hori_minEnergy(hori_seam(i+1),i);
            [val,direction] = min([above, bottom, right]);
        end
        if direction == 3
            min_idx = hori_seam(i+1)+1;
        elseif direction == 2
            min_idx = hori_seam(i+1);
        else
            min_idx = hori_seam(i+1)-1;
        end
        hori_seam(1,i) = min_idx;
    end
end

function seamline = seamMark (seamimg, seam, type)
    
    if strcmp(type, 'VERTICAL')
        for i =1:size(seamimg,1)
            seamimg(i,seam(i),1)=0;
            seamimg(i,seam(i),2)=0;
            seamimg(i,seam(i),3)=255;
        end
    else
        for i = 1:size(seamimg,2)
            seamimg(seam(i),i,1)=0;
            seamimg(seam(i),i,2)=0;
            seamimg(seam(i),i,3)=255;
        end
    end
    seamline = seamimg;
end

function [result_img, result_energyimg] = removeHeight(Img, energyImg)
min_Energymap = hori_minEnergymap(energyImg);
horizontal_seam = optimal_horiseam(min_Energymap);
rows = size(energyImg,1);
cols = size(energyImg,2);
for i = 1:cols
    Img(horizontal_seam(i):rows-1,i,:) = Img(horizontal_seam(i)+1:rows,i,:);
end
result_img = Img(1:rows-1,:,:);
result_energyimg = energy(result_img);
end

function [result_img, result_energyimg] = removeWidth(Img, energyImg)
min_Energymap = vert_minEnergymap(energyImg);
vertical_seam = optimal_vertseam(min_Energymap);
rows = size(energyImg,1);
cols = size(energyImg,2);
for i = 1:rows
    Img(i,vertical_seam(i,1):cols-1,:) = Img(i,vertical_seam(i,1)+1:cols,:);
end
result_img = Img(:,1:cols-1,:);
result_energyimg = energy(result_img);
end



