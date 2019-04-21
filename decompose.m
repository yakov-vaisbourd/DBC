function [Trans, TransSizes,g_num] = decompose(m,n,TV,num)
% -----------------------------------------------------------------------
% Copyright (2016): A.Beck, L.Tetruashvili, Yakov Vaisbourd and A.
% Shemtov
% 
% dbc_tv is distributed under the terms of 
% the GNU General Public License 2.0.
% 
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ----------------------------------------------------------------------
%%
modd = rem(m,2);
nodd = rem(n,2);

% Set the linear transformations.

if n==1
    if strcmpi(TV, 'l1')
        if num<2
            error 'For l1 TV num should be greater or equal to 2';
        end
        g_num = num;
        TransSizes = 2*floor((m-(-(num-2):1)')/num);        
        max_size = max(TransSizes);
        Trans = zeros(max_size,g_num);

        for i=1:g_num
            Trans(1:2:TransSizes(i),i) = i:num:((TransSizes(i)/2)*num-2+i);
            Trans(2:2:TransSizes(i),i) = (i:num:((TransSizes(i)/2)*num-2+i))+1;
        end
    else
        if num<3
            error 'For iso TV num should be greater or equal to 3';
        end
        g_num = num;
        TransSizes = 3*floor((m-(-(num-3):2)')/num);
        max_size = max(TransSizes);
        Trans = zeros(max_size,g_num);
               
        for i=1:g_num
            Trans(1:3:TransSizes(i),i) = i:num:((TransSizes(i)/3)*num-3+i);
            Trans(2:3:TransSizes(i),i) = (i:num:((TransSizes(i)/3)*num-3+i))+1;
            Trans(3:3:TransSizes(i),i) = (i:num:((TransSizes(i)/3)*num-3+i))+2;
        end
    end   
else %n>1
    if strcmpi(TV, 'l1')
        g_num = 4;        
        TransSizes = zeros(g_num,1);       
        TransSizes(1:2) = 2*n*floor((m-(0:1)')/2);
        TransSizes(3:4) =2*m*floor((n-(0:1)')/2);
        
        max_size = max(TransSizes);
        Trans = zeros(max_size,g_num);
        

        Trans(1:TransSizes(1),1)= uint32(repmat((1:(m-modd))',n,1)+...
                                      kron((0:(n-1))',ones(m-modd,1))*m);
        Trans(1:TransSizes(2),2)= uint32(repmat((2:(m-1+modd))',n,1)+...
                                      kron((0:(n-1))',ones(m-2+modd,1))*m);
        Trans(1:TransSizes(3),3)= uint32(repmat((1:m:(m*(n-nodd)))',m,1)+...
                                      kron((0:(m-1))',ones(n-nodd,1)));
        Trans(1:TransSizes(4),4)= uint32(repmat(((m+1):m:(m*(n-1+nodd)))',m,1)+...
                                      kron((0:(m-1))',ones(n-2+nodd,1)));                  
    else % TV='iso'
        g_num = 3;
        TransSizes = zeros(g_num,2);
        NinRow = zeros(3);
        % Number of r-shaped triples in a column for each level of 
        % diagonal (that might be cyclic) shift.
        NinCol = [floor((m+1)/3);floor(m/3);floor((m-1)/3)];
        for i=1:3 % column shift = equivalent to change in 
                  % transformation index                                
            % Number of r-shaped triples in a rwo for each level of
            %  diagonal (that might be cyclic) shift.   % i=1 i=2 i=3     
            NinRow(:,i) = [floor((n+1-rem(i-1,3))/3);...  % n+1 n   n-1   m+1
                           floor((n+1-rem(i,3))/3);...    % n   n-1 n+1   m
                           floor((n+1-rem(i+1,3))/3)];    % n-1 n+1 n     m-1

            TransSizes(i,1) = 3*NinCol'*NinRow(:,i); % Compute the length of 
                                                 % the transformed 1D 
                                                 % vector.

            TransSizes(i,2) = 2*(NinCol(rem(3+rem(n,3)-i,3)+1)+NinRow(rem(2+rem(m,3),3)+1,i));
        end    

        max_size = max(max(TransSizes));
        Trans = zeros(max_size,g_num,2);

        % The iso TV is a sum of 3-dimensional vector norms.
        % Each such vector is defined by the r-shaped triple 
        % (x_{i,j+1},x_{i,j},x_{i+1,1}) for each pair of indices i<m,j<n. 
        % The following code produce a one to one transformation of the 
        % entities from 2D to 1D and vise versa according to the type of 
        % decomposition. Throughout of the documentation of this code the 
        % first entity of the triple x_{i,j+1} will be referred as the 
        % column shift entity (cs-entity), the middle one as the base 
        % entity (b-entity) and the last one will be referred as the row 
        % shift entity (rs-entity).

        for i=1:3 %column shift = equivalent to change in 
            % transformation index.

            % Mark the indices that will contain the values of the 
            % rs-entities in the transformed 1D vector.
            indx = zeros(max_size,1);indx(1:3:TransSizes(i))=1;
            indx=logical(indx);
            % Set the boundaries of the domain of each level of 
            % diagonal shift.
            indx_sep = cumsum([1;3*NinCol.*NinRow(:,i)]);


            % rs-entity start index in the left most column for each 
            % level of diagonal shift.
            csSiF = (1:3)';
            csSiF = (rem(i+csSiF+1,3)+1)*m+csSiF; 

            % rs-entity end index in the left most column for each 
            % level of diagonal shift.        
            csEiF = csSiF + 3*(NinCol-1);       

            % b&cs-couple entities start index in the left most column 
            % for each level of diagonal shift.
            brsSiF = (1:3)';
            brsSiF = rem(1+brsSiF+i,3)*m+brsSiF;

            % b&cs-couple entities end index in the left most column 
            % for each level of diagonal shift.        
            bcsEi = brsSiF + 3*(NinCol-1);



            for j=1:3 %row shift
                curr_indx = (indx_sep(j):(indx_sep(j+1)-1))';
                curr_indx = indx(curr_indx).*curr_indx;
                curr_indx = curr_indx(curr_indx>0);
                Trans(curr_indx,i,1) = ...
                     repmat((csSiF(j):3:csEiF(j))',NinRow(j,i),1) + ...
                        kron((0:3:(3*NinRow(j,i)-1))',ones(NinCol(j),1))*m;
                Trans(curr_indx+1,i,1) = ...
                     repmat((brsSiF(j):3:bcsEi(j))',NinRow(j,i),1) + ...
                        kron((0:3:(3*NinRow(j,i)-1))',ones(NinCol(j),1))*m;        
                Trans(curr_indx+2,i,1) = ...
                     repmat(1+(brsSiF(j):3:bcsEi(j))',NinRow(j,i),1) + ...
                        kron((0:3:(3*NinRow(j,i)-1))',ones(NinCol(j),1))*m;
            end

            colIndx = rem(3+rem(n,3)-i,3)+1;
            brsSiL = m*(n-1)+colIndx;
            brsEiL = brsSiL+3*(NinCol(colIndx)-1);
            rowIndx = rem(2+rem(m,3),3)+1;
            bcsSiL = (rem(rowIndx+i-2,3)+1)*m;
            bcsEiL = bcsSiL+3*m*(NinRow(rowIndx,i)-1);


            Trans(1:2:TransSizes(i,2),i,2) = ...
                    [(brsSiL:3:brsEiL)';(bcsSiL:3*m:bcsEiL)'];
            Trans(2:2:TransSizes(i,2),i,2) = ...
                    [(brsSiL+1:3:brsEiL+1)';(bcsSiL+m:3*m:bcsEiL+m)'];

        end    

    end

end


if num~=0 &&  (n>1)
    error('Not supported.');
end

end


