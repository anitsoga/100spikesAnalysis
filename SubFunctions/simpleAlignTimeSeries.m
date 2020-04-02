function [Ar, dxs, dys] = simpleAlignTimeSeries (A)

iterate=1;
iterLimit=10;

meanIMG = mean(A,3);
sz = size(A);

disp('First Pass...')
parfor i = 1:sz(3);
%     if mod(i,sz(3)/10)==1
%         fprintf('. ')
%     end
    [dx,dy] = fftalign(A(:,:,i),meanIMG);
    A(:,:,i) = circshift(A(:,:,i),[dx, dy]);
    dxi(i) = dx; %This pass offsets
    dyi(i) = dy;
    
    dxs(i) = dx; %total offsets
    dys(i) = dy;
    %     meanIMG = mean(A,3);
end

if iterate
    iterCounter = 0;
    while ~ all([dxi dyi]==0) & iterCounter <iterLimit;
        iterCounter=iterCounter+1;
        disp(['Itterating Alignment Pass: ' num2str(iterCounter) '. last remaining correction : ' num2str(mean(dxi)) ' ' num2str(mean(dyi)) ] );
        meanIMG = mean(A,3);
        
        parfor i = 1:sz(3);
%             if mod(i,sz(3)/10)==1
%                 fprintf('. ')
%             end
            [dx,dy] = fftalign(A(:,:,i),meanIMG);
            A(:,:,i) = circshift(A(:,:,i),[dx, dy]);
            dxi(i) = dx; %This pass offsets
            dyi(i) = dy;
            
            dxs(i) = dxs(i)+dx; %total offsets
            dys(i) = dys(i)+dy;
            %     meanIMG = mean(A,3);
        end
        
    end
end

disp(['Done. final mean corection: ' num2str(mean(dxs)) ' ' num2str(mean(dys))]);

Ar = A;
% fprintf('\n')