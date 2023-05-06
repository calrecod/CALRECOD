function plotRCWallDesign(H,m1,m2,toe,heel,hf,b,D,alfa,beta)

%-------------------------------------------------------------------------
% Syntax:
% plotRCWallDesign(H,m1,m2,toe,heel,hf,b,D,alfa,beta)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%-------------------------------------------------------------------------
% PURPOSE: To plot the cross-section of a retaining RC wall
%
% INPUT:  toe:                  is the cross-sectional length of the wall's
%                               foot
%
%         heel:                 is the cross-sectional length dimension of
%                               the wall's heel
%
%         hf:                   is the cross-sectional width of the wall's 
%                               heel and foot
%
%         b:                    is the upper cross-section wall's trunk
%                               width
%
%         H:                    is the wall's stem height dimension
%
%         D:                    is the soil's bakc fill depth
%
%         m1, m2:               are the wall's dowels' from and back slopes
%
%         beta:                 is the front fill soil's upper-grade angle
%
%         alfa:                 is the back fill soil's upper-grade angle
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-04-29
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

%% Wall's cross-section
secWall=[toe+1/m2*H, hf+H;
         toe+1/m2*H+b, hf+H;
         toe+1/m2*H+b+1/m1*H, hf;
         toe+1/m2*H+b+1/m1*H+heel, hf;
         toe+1/m2*H+b+1/m1*H+heel, 0;
         0, 0;
         0, hf;
         toe, hf;
         toe+1/m2*H, hf+H];
     
x=secWall(:,1);
y=secWall(:,2);
figure(4)
plot(x,y,'k -','linewidth',2.0)
legend('Wall´s borders')
hold on
xlabel('x´')
ylabel('y´')
title('Retaining RC wall´ cross-section')
axis([-(toe)*1.1 (toe+1/m2*H+b+1/m1*H+heel)*1.1 -(hf)*1.1 (H+hf)*1.15])
hold on

%% Back and front fill
if D>0
    backfill=[0,D+hf+tan(deg2rad(alfa))*(toe+1/m2*H);
          toe+1/m2*H,D+hf];
      
    plot(backfill(:,1),backfill(:,2),'r -','linewidth',1.5,'DisplayName',...
        'Back Fill')
end
frontfill=[toe+1/m2*H+b,H+hf;
          toe+1/m2*H+b+1/m1*H+heel,H+hf+tan(deg2rad(beta))*(1/m1*H+heel)];
      
plot(frontfill(:,1),frontfill(:,2),'b -','linewidth',1.5,'DisplayName',...
    'Front Fill')
