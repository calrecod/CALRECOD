function elemConc=casoConcretoRecRot(a,fdpc,b,h,RotCornerSec,rotCP,gamma)
%------------------------------------------------------------------------
% Syntax:
% elemConc=casoConcretoRecRot(a,fdpc,b,h,RotCornerSec,rotCP,gamma)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: to compute the contribution of resistance of the concrete
% compression zone of a rectangular rotated concrete cross-section,
% regarding axial bending forces.
% 
% INPUT:  a:                is the reduced depth of the neutral axis of the
%                           cross-section in question
%
%         fdpc:             is the factored value of the concrete
%                           compressive concrete f'c as 0.85 * f'c, 
%                           according to the ACI 318 code
%
%         b,h:              are the cross-section dimensions (width and
%                           height, respectively)
%
%         RotCornerSec:     are the cross-section corner coordinates.
%                           Size: 4 x 2, in format:
%
%                           [x-top-right, y-top-right;
%                            x-top-left, y-top-left;
%                            x-bottom-left, y-bottom-left;
%                            x-bottom-right, y-bottom-right]
%
%         rotCP:            are the Plastic Center depth values of the
%                           rotated cross-section
%
%         gamma:            is the angle of rotation of the cross-section
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-02-05
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro, Mexico
%------------------------------------------------------------------------
areaCo=0;
% Slope for boundary 1
m1=(RotCornerSec(1,2)-RotCornerSec(2,2))/...
    (RotCornerSec(1,1)-RotCornerSec(2,1));
% Slope for boundary 2
m2=(RotCornerSec(4,2)-RotCornerSec(1,2))/...
       (RotCornerSec(4,1)-RotCornerSec(1,1));
% Slope for boundary 3
m3=(RotCornerSec(3,2)-RotCornerSec(2,2))/...
       (RotCornerSec(3,1)-RotCornerSec(2,1));
% Slope for boundary 4
m4=(RotCornerSec(4,2)-RotCornerSec(3,2))/...
       (RotCornerSec(4,1)-RotCornerSec(3,1));

% Integration with trapezoidal rule to determine the
% concrete area in compression
n=100;
dh=a/n;
mc=0;
if gamma>=0 && gamma<90 % Case 1
    if RotCornerSec(1,2)>=RotCornerSec(3,2) % case 1.1
        sh1=RotCornerSec(2,2)-RotCornerSec(1,2);
        sh2=RotCornerSec(2,2)-RotCornerSec(3,2);
        y=RotCornerSec(2,2)-a;
        if a<sh1
            x1=RotCornerSec(2,1)+1/m1*(y-RotCornerSec(2,2));
            x3=RotCornerSec(2,1)+1/m3*(y-RotCornerSec(2,2));
            areaCo=(x1-x3)*dh/2;

            mc=areaCo*a;
        elseif a>=sh1 && a<sh2
        
            x2=RotCornerSec(1,1)+1/m2*(y-RotCornerSec(1,2));
            x3=RotCornerSec(2,1)+1/m3*(y-RotCornerSec(2,2));
            areaCo=(x2-x3)*dh/2;

            mc=areaCo*a;
        elseif a>=sh2
        
            x4=RotCornerSec(3,1)+1/m4*(y-RotCornerSec(3,2));
            x2=RotCornerSec(1,1)+1/m2*(y-RotCornerSec(1,2));
            areaCo=(x2-x4)*dh/2;

            mc=areaCo*a;
        end
    elseif RotCornerSec(1,2)<RotCornerSec(3,2) % case 1.2
        sh1=RotCornerSec(2,2)-RotCornerSec(3,2);
        sh2=RotCornerSec(2,2)-RotCornerSec(1,2);
        y=RotCornerSec(2,2)-a;
        if a<sh1

            x1=RotCornerSec(2,1)+1/m1*(y-RotCornerSec(2,2));
            x3=RotCornerSec(2,1)+1/m3*(y-RotCornerSec(2,2));
            areaCo=(x1-x3)*dh/2;

            mc=areaCo*a;
        elseif a>=sh1 && a<sh2
        
            x1=RotCornerSec(2,1)+1/m1*(y-RotCornerSec(2,2));
            x4=RotCornerSec(3,1)+1/m4*(y-RotCornerSec(3,2));
            areaCo=(x1-x4)*dh/2;

            mc=areaCo*a;
        elseif a>=sh2
        
            x4=RotCornerSec(3,1)+1/m4*(y-RotCornerSec(3,2));
            x2=RotCornerSec(1,1)+1/m2*(y-RotCornerSec(1,2));
            areaCo=(x2-x4)*dh/2;

            mc=areaCo*a;
        end
    end
    % Main integration loop
    for i=1:n-1
        xa=dh*i;
        ya=RotCornerSec(2,2)-xa;
        if RotCornerSec(1,2)>=RotCornerSec(3,2) % Case 1.1
            
            if xa<sh1

                x1=RotCornerSec(2,1)+1/m1*(ya-RotCornerSec(2,2));
                x3=RotCornerSec(2,1)+1/m3*(ya-RotCornerSec(2,2));
                areaCo=areaCo+(x1-x3)*dh;

                mc=mc+(x1-x3)*dh*xa;
            elseif xa>=sh1 && xa<sh2

                x2=RotCornerSec(1,1)+1/m2*(ya-RotCornerSec(1,2));
                x3=RotCornerSec(2,1)+1/m3*(ya-RotCornerSec(2,2));
                areaCo=areaCo+(x2-x3)*dh;

                mc=mc+(x2-x3)*dh*xa;
            elseif xa>=sh2

                x4=RotCornerSec(3,1)+1/m4*(ya-RotCornerSec(3,2));
                x2=RotCornerSec(1,1)+1/m2*(ya-RotCornerSec(1,2));
                areaCo=areaCo+(x2-x4)*dh;

                mc=mc+(x2-x4)*dh*xa;
            end
        elseif RotCornerSec(1,2)<RotCornerSec(3,2) % Case 1.2

            if xa<sh1

                x1=RotCornerSec(2,1)+1/m1*(ya-RotCornerSec(2,2));
                x3=RotCornerSec(2,1)+1/m3*(ya-RotCornerSec(2,2));
                areaCo=areaCo+(x1-x3)*dh;

                mc=mc+(x1-x3)*dh*xa;

            elseif xa>=sh1 && xa<sh2

                x4=RotCornerSec(3,1)+1/m4*(ya-RotCornerSec(3,2));
                x1=RotCornerSec(2,1)+1/m1*(ya-RotCornerSec(2,2));
                areaCo=areaCo+(x1-x4)*dh;

                mc=mc+(x1-x4)*dh*xa;
            elseif xa>=sh2

                x4=RotCornerSec(3,1)+1/m4*(ya-RotCornerSec(3,2));
                x2=RotCornerSec(1,1)+1/m2*(ya-RotCornerSec(1,2));
                areaCo=areaCo+(x2-x4)*dh;

                mc=mc+(x2-x4)*dh*xa;
            end
        end
    end
elseif gamma>=90 && gamma<180 % Case 2
    if RotCornerSec(2,2)>=RotCornerSec(4,2) % case 2.1
        sh1=RotCornerSec(3,2)-RotCornerSec(2,2);
        sh2=RotCornerSec(3,2)-RotCornerSec(4,2);
        y=RotCornerSec(3,2)-a;
        if a<sh1

            x4=RotCornerSec(3,1)+1/m4*(y-RotCornerSec(3,2));
            x3=RotCornerSec(2,1)+1/m3*(y-RotCornerSec(2,2));
            areaCo=(x3-x4)*dh/2;

            mc=areaCo*a;
        elseif a>=sh1 && a<sh2
        
            x4=RotCornerSec(3,1)+1/m4*(y-RotCornerSec(3,2));
            x1=RotCornerSec(2,1)+1/m1*(y-RotCornerSec(2,2));
            areaCo=(x1-x4)*dh/2;

            mc=areaCo*a;
        elseif a>=sh2
        
            x1=RotCornerSec(2,1)+1/m1*(y-RotCornerSec(2,2));
            x2=RotCornerSec(1,1)+1/m2*(y-RotCornerSec(1,2));
            areaCo=(x1-x2)*dh/2;

            mc=areaCo*a;
        end
    elseif RotCornerSec(2,2)<RotCornerSec(4,2) % case 2.2
        sh1=RotCornerSec(3,2)-RotCornerSec(4,2);
        sh2=RotCornerSec(3,2)-RotCornerSec(2,2);
        y=RotCornerSec(3,2)-a;
        if a<sh1

            x4=RotCornerSec(3,1)+1/m4*(y-RotCornerSec(3,2));
            x2=RotCornerSec(1,1)+1/m2*(y-RotCornerSec(1,2));
            areaCo=(x2-x4)*dh/2;

            mc=areaCo*a;
        elseif a>=sh1 && a<sh2
        
            x2=RotCornerSec(1,1)+1/m2*(y-RotCornerSec(1,2));
            x3=RotCornerSec(2,1)+1/m3*(y-RotCornerSec(2,2));
            areaCo=(x2-x3)*dh/2;

            mc=areaCo*a;
        elseif a>=sh2
        
            x1=RotCornerSec(2,1)+1/m1*(y-RotCornerSec(2,2));
            x3=RotCornerSec(2,1)+1/m3*(y-RotCornerSec(2,2));
            areaCo=(x1-x3)*dh/2;

            mc=areaCo*a;
        end
    end
    % Main integration loop
    for i=1:n-1
        xa=dh*i;
        ya=RotCornerSec(3,2)-xa;
        if RotCornerSec(2,2)>=RotCornerSec(4,2) % Case 2.1
            
            if xa<sh1

                x4=RotCornerSec(3,1)+1/m4*(ya-RotCornerSec(3,2));
                x3=RotCornerSec(2,1)+1/m3*(ya-RotCornerSec(2,2));
                areaCo=areaCo+(x3-x4)*dh;

                mc=mc+(x3-x4)*dh*xa;
            elseif xa>=sh1 && xa<sh2

                x4=RotCornerSec(3,1)+1/m4*(ya-RotCornerSec(3,2));
                x1=RotCornerSec(2,1)+1/m1*(ya-RotCornerSec(2,2));
                areaCo=areaCo+(x1-x4)*dh;

                mc=mc+(x1-x4)*dh*xa;
            elseif xa>=sh2

                x1=RotCornerSec(2,1)+1/m1*(ya-RotCornerSec(2,2));
                x2=RotCornerSec(1,1)+1/m2*(ya-RotCornerSec(1,2));
                areaCo=areaCo+(x1-x2)*dh;

                mc=mc+(x1-x2)*dh*xa;
            end
        elseif RotCornerSec(2,2)<RotCornerSec(4,2) % Case 2.2

            if xa<sh1

                x4=RotCornerSec(3,1)+1/m4*(ya-RotCornerSec(3,2));
                x3=RotCornerSec(2,1)+1/m3*(ya-RotCornerSec(2,2));
                areaCo=areaCo+(x3-x4)*dh;

                mc=mc+(x3-x4)*dh*xa;

            elseif xa>=sh1 && xa<sh2

                x2=RotCornerSec(1,1)+1/m2*(ya-RotCornerSec(1,2));
                x3=RotCornerSec(2,1)+1/m3*(ya-RotCornerSec(2,2));
                areaCo=areaCo+(x3-x2)*dh;

                mc=mc+(x3-x2)*dh*xa;
            elseif xa>=sh2

                x1=RotCornerSec(2,1)+1/m1*(ya-RotCornerSec(2,2));
                x2=RotCornerSec(1,1)+1/m2*(ya-RotCornerSec(1,2));
                areaCo=areaCo+(x1-x2)*dh;

                mc=mc+(x1-x2)*dh*xa;
            end
        end
    end
    
elseif gamma>=180 && gamma<270 % Case 3
    if RotCornerSec(3,2)>=RotCornerSec(1,2) % Case 3.1
        sh1=RotCornerSec(4,2)-RotCornerSec(3,2);
        sh2=RotCornerSec(4,2)-RotCornerSec(1,2);
        y=RotCornerSec(4,2)-a;
        if a<sh1

            x4=RotCornerSec(3,1)+1/m4*(y-RotCornerSec(3,2));
            x2=RotCornerSec(1,1)+1/m2*(y-RotCornerSec(1,2));
            areaCo=(x4-x2)*dh/2;

            mc=areaCo*a;
        elseif a>=sh1 && a<sh2

            x2=RotCornerSec(1,1)+1/m2*(y-RotCornerSec(1,2));
            x3=RotCornerSec(2,1)+1/m3*(y-RotCornerSec(2,2));
            areaCo=(x3-x2)*dh/2;

            mc=areaCo*a;
        elseif a>=sh2

            x1=RotCornerSec(2,1)+1/m1*(y-RotCornerSec(2,2));
            x3=RotCornerSec(2,1)+1/m3*(y-RotCornerSec(2,2));
            areaCo=(x3-x1)*dh/2;

            mc=areaCo*a;
        end
    elseif RotCornerSec(3,2)<RotCornerSec(1,2) % Case 3.2
        sh1=RotCornerSec(4,2)-RotCornerSec(1,2);
        sh2=RotCornerSec(4,2)-RotCornerSec(3,2);
        y=RotCornerSec(4,2)-a;
        if a<sh1

            x4=RotCornerSec(3,1)+1/m4*(y-RotCornerSec(3,2));
            x2=RotCornerSec(1,1)+1/m2*(y-RotCornerSec(1,2));
            areaCo=(x4-x2)*dh/2;

            mc=areaCo*a;
        elseif a>=sh1 && a<sh2

            x4=RotCornerSec(3,1)+1/m4*(y-RotCornerSec(3,2));
            x1=RotCornerSec(2,1)+1/m1*(y-RotCornerSec(2,2));
            areaCo=(x4-x1)*dh/2;

            mc=areaCo*a;
        elseif a>=sh2

            x1=RotCornerSec(2,1)+1/m1*(y-RotCornerSec(2,2));
            x3=RotCornerSec(2,1)+1/m3*(y-RotCornerSec(2,2));
            areaCo=(x3-x1)*dh/2;

            mc=areaCo*a;
        end
    end
    % Main integration loop
    for i=1:n-1
        xa=dh*i;
        ya=RotCornerSec(4,2)-xa;
        if RotCornerSec(3,2)>=RotCornerSec(1,2) % Case 3.1
            
            if xa<sh1

                x4=RotCornerSec(3,1)+1/m4*(ya-RotCornerSec(3,2));
                x2=RotCornerSec(1,1)+1/m2*(ya-RotCornerSec(1,2));
                areaCo=areaCo+(x4-x2)*dh;

                mc=mc+(x4-x2)*dh*xa;
            elseif xa>=sh1 && xa<sh2

                x2=RotCornerSec(1,1)+1/m2*(ya-RotCornerSec(1,2));
                x3=RotCornerSec(2,1)+1/m3*(ya-RotCornerSec(2,2));
                areaCo=areaCo+(x3-x2)*dh;

                mc=mc+(x3-x2)*dh*xa;
            elseif xa>=sh2

                x1=RotCornerSec(2,1)+1/m1*(ya-RotCornerSec(2,2));
                x3=RotCornerSec(2,1)+1/m3*(ya-RotCornerSec(2,2));
                areaCo=areaCo+(x3-x1)*dh;

                mc=mc+(x3-x1)*dh*xa;
            end
        elseif RotCornerSec(2,2)<RotCornerSec(4,2) % Case 4.2

            if xa<sh1

                x4=RotCornerSec(3,1)+1/m4*(ya-RotCornerSec(3,2));
                x2=RotCornerSec(1,1)+1/m2*(ya-RotCornerSec(1,2));
                areaCo=areaCo+(x4-x2)*dh;

                mc=mc+(x4-x2)*dh*xa;

            elseif xa>=sh1 && xa<sh2

                x4=RotCornerSec(3,1)+1/m4*(ya-RotCornerSec(3,2));
                x1=RotCornerSec(2,1)+1/m1*(ya-RotCornerSec(2,2));
                areaCo=areaCo+(x4-x1)*dh;

                mc=mc+(x4-x1)*dh*xa;
            elseif xa>=sh2

                x1=RotCornerSec(2,1)+1/m1*(ya-RotCornerSec(2,2));
                x3=RotCornerSec(2,1)+1/m3*(ya-RotCornerSec(2,2));
                areaCo=areaCo+(x3-x1)*dh;

                mc=mc+(x3-x1)*dh*xa;
            end
        end
    end
    
elseif gamma>=270 && gamma<=360 % Case 4

    if RotCornerSec(2,2)>=RotCornerSec(4,2) % Case 4.1
        sh1=RotCornerSec(1,2)-RotCornerSec(2,2);
        sh2=RotCornerSec(1,2)-RotCornerSec(4,2);
        y=RotCornerSec(1,2)-a;
        if a<sh1

            x1=RotCornerSec(2,1)+1/m1*(y-RotCornerSec(2,2));
            x2=RotCornerSec(1,1)+1/m2*(y-RotCornerSec(1,2));
            areaCo=(x2-x1)*dh/2;

            mc=areaCo*a;
        elseif a>=sh1 && a<sh2

            x2=RotCornerSec(1,1)+1/m2*(y-RotCornerSec(1,2));
            x3=RotCornerSec(2,1)+1/m3*(y-RotCornerSec(2,2));
            areaCo=(x2-x3)*dh/2;

            mc=areaCo*a;
        elseif a>=sh2

            x4=RotCornerSec(3,1)+1/m4*(y-RotCornerSec(3,2));
            x3=RotCornerSec(2,1)+1/m3*(y-RotCornerSec(2,2));
            areaCo=(x4-x3)*dh/2;

            mc=areaCo*a;
        end
    elseif RotCornerSec(2,2)<RotCornerSec(4,2) % Case 4.2
        sh1=RotCornerSec(1,2)-RotCornerSec(4,2);
        sh2=RotCornerSec(1,2)-RotCornerSec(2,2);
        y=RotCornerSec(1,2)-a;
        if a<sh1

            x1=RotCornerSec(2,1)+1/m1*(y-RotCornerSec(2,2));
            x2=RotCornerSec(1,1)+1/m2*(y-RotCornerSec(1,2));
            areaCo=(x2-x1)*dh/2;

            mc=areaCo*a;
        elseif a>=sh1 && a<sh2

            x4=RotCornerSec(3,1)+1/m4*(y-RotCornerSec(3,2));
            x1=RotCornerSec(2,1)+1/m1*(y-RotCornerSec(2,2));
            areaCo=(x4-x1)*dh/2;

            mc=areaCo*a;
        elseif a>=sh2

            x4=RotCornerSec(3,1)+1/m4*(y-RotCornerSec(3,2));
            x3=RotCornerSec(2,1)+1/m3*(y-RotCornerSec(2,2));
            areaCo=(x4-x3)*dh/2;

            mc=areaCo*a;
        end
    end
    % Main integration loop
    for i=1:n-1
        xa=dh*i;
        if RotCornerSec(2,2)>=RotCornerSec(4,2) % Case 4.1
            ya=RotCornerSec(1,2)-xa;
            if xa<sh1

                x1=RotCornerSec(2,1)+1/m1*(ya-RotCornerSec(2,2));
                x2=RotCornerSec(1,1)+1/m2*(ya-RotCornerSec(1,2));
                areaCo=areaCo+(x2-x1)*dh;

                mc=mc+(x2-x1)*dh*xa;
            elseif xa>=sh1 && xa<sh2

                x2=RotCornerSec(1,1)+1/m2*(ya-RotCornerSec(1,2));
                x3=RotCornerSec(2,1)+1/m3*(ya-RotCornerSec(2,2));
                areaCo=areaCo+(x2-x3)*dh;

                mc=mc+(x2-x3)*dh*xa;
            elseif xa>=sh2

                x4=RotCornerSec(3,1)+1/m4*(ya-RotCornerSec(3,2));
                x3=RotCornerSec(2,1)+1/m3*(ya-RotCornerSec(2,2));
                areaCo=areaCo+(x4-x3)*dh;

                mc=mc+(x4-x3)*dh*xa;
            end
        elseif RotCornerSec(2,2)<RotCornerSec(4,2) % Case 4.2
            ya=RotCornerSec(1,2)-xa;
            if xa<sh1

                x1=RotCornerSec(2,1)+1/m1*(ya-RotCornerSec(2,2));
                x2=RotCornerSec(1,1)+1/m2*(ya-RotCornerSec(1,2));
                areaCo=areaCo+(x2-x1)*dh;

                mc=mc+(x2-x1)*dh*xa;

            elseif xa>=sh1 && xa<sh2

                x4=RotCornerSec(3,1)+1/m4*(ya-RotCornerSec(3,2));
                x1=RotCornerSec(2,1)+1/m1*(ya-RotCornerSec(2,2));
                areaCo=areaCo+(x4-x1)*dh;

                mc=mc+(x4-x1)*dh*xa;
            elseif xa>=sh2

                x4=RotCornerSec(3,1)+1/m4*(ya-RotCornerSec(3,2));
                x3=RotCornerSec(2,1)+1/m3*(ya-RotCornerSec(2,2));
                areaCo=areaCo+(x4-x3)*dh;

                mc=mc+(x4-x3)*dh*xa;
            end
        end
    end
end
if areaCo==0
    areaCo=1e-6;
end
yc=mc/areaCo; % depth of the concrete compressión zone's centroid

ha=max(RotCornerSec(:,2))-min(RotCornerSec(:,2));
if (a>ha)
    areaCo=b*h;
    
    yc=max(RotCornerSec(:,2));
end

frc=-areaCo*fdpc;
mrc=-frc*(rotCP(1)-yc);

elemConc=[frc mrc];
