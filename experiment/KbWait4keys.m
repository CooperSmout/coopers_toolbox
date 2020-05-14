function keyOut = KbWait4keys(responseKeys)

found=0;

while ~found

    [~, keyCode] = KbWait([],2);
    
    key = find(keyCode,1);
    
    if find(responseKeys==key)
        
        found=1;
        keyOut = find(responseKeys==key);
    end
    
end
        