

function mondrianNoise(h,w,maxsize,filename)

im = ones(h,w,3);

for ii = 1:10000;
    
    locY = randi(h);
    H = randi(maxsize);
    Y = max((locY-H),1) : min((locY+H),h);
    
    locX = randi(w);
    W = randi(maxsize);
    X = max((locX-W),1) : min((locX+W),w);
    
    im(Y,X,1) = rand;
    im(Y,X,2) = rand;
    im(Y,X,3) = rand;
    
end

image(im);
imwrite(im, filename, 'bmp')