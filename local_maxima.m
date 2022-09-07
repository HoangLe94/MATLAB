function output = local_maxima(P)
    %Not useful function
    output = [];

    for i = 1:length(P(:,1))
        for j = 1:length(P(1,:))
            switch i
                case 1
                    switch j
                        case 1
                            if P(i,j) > P(i,j+1) & P(i,j) > P(i+1,j) & P(i,j) > P(i+1,j+1)
                                output = [output; [P(i,j) i j]];
                            end
                        case length(P(1,:))
                            if P(i,j) > P(i,j-1) & P(i,j) > P(i+1,j) & P(i,j) > P(i+1,j-1)
                                output = [output; [P(i,j) i j]];
                            end
                        otherwise
                            if P(i,j) > P(i,j+1) & P(i,j) > P(i+1,j) & P(i,j) > P(i+1,j+1) & P(i,j) > P(i,j-1) & P(i,j) > P(i+1,j-1) 
                                output = [output; [P(i,j) i j]];
                            end
                    end
                    
                case length(P(:,1))
                    switch j
                        case 1
                        case length(P(1,:))
                        otherwise
                    end
                    
                otherwise
                    switch j
                        case 1
                        case length(P(1,:))
                        otherwise
                    end
            end
        end
    end
end

