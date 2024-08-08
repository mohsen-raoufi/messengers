function isInside = isInSenseRange(relPos,senseRange)

isInside = false;

if(abs(relPos(1))<senseRange)
    
    if(abs(relPos(2))<senseRange)
        
        if(norm(relPos)<senseRange)
            isInside = true;
        end
        
    end
    
end

end