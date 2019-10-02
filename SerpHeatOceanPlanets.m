function Planets = SerpHeatOceanPlanets(Planets)

m2km = 1e-3;
W2mW = 1e3;
yr2Gyr = 1e-9;
yr2s = 3.15576*1e7;
Age_s = yr2s*4.5e9;

HeadSize = 45;
CaptSize = 28;
LW = 1.25;
DERVS = 1;

% [MOR_F_W_m2,MOR_FH2_molecules_cm2_s]=MOR_SerpProd;
% MOR_FH2_moles_yr = MOR_FH2_molecules_cm2_s.*3.15569e7./6.022e23.*4*pi*(1e2*6678.1e3).^2;    


lPlanet = length(Planets);
if lPlanet
for iPlanet = 1:length(Planets)
%         if ~isempty(Planets(iPlanet).z_bd_m)
%              for jk = 1:length(Planets(iPlanet).d_ocean_km)
%                 Planets(iPlanet).surf_area_seafloor(jk) = 4*pi*(Planets(iPlanet).R_m-Planets(iPlanet).d_ocean_km(jk)*1e3)^2; 
%                 [Planets(iPlanet).F_serp_bd_W_m2,Planets(iPlanet).FH2_bd_molecules_cm2_s,~] =...
%                     getHeatFlux(Planets(iPlanet).R_m,Planets(iPlanet).d_ocean_km(jk)*1e3,Planets(iPlanet).z_bd_m(jk,:,:),DERVS);
%                 Planets(iPlanet).FH2_bd_moles_yr(jk,:) = ...
%                     Planets(iPlanet).FH2_bd_molecules_cm2_s(jk,:).*3.15569e7./6.022e23.*4*pi*(1e2*Planets(iPlanet).R_m).^2;
%                 % the next part should really be consolidated as a function with the
%                 % repeated section below
%                if isfield(Planets(iPlanet),'t_recharge_Mya')
%                    if ~isempty(Planets(iPlanet).t_recharge_Mya)
%                       t_ya = Planets(iPlanet).t_yr;
%                        t_inds = find(t_ya>Planets(iPlanet).t_recharge_Mya*1e6);
%                        add_inds = [t_inds(1)-2 t_inds(1)-1];
%                        F_serp_bd_recharge_W_m2 = sum(Planets(iPlanet).F_serp_bd_W_m2(jk,t_inds))*(sum(t_ya(t_inds)))/(t_ya(t_inds(1))-t_ya(t_inds(1)-1))/2;
%                        Planets(iPlanet).F_serp_bd_recharge_W_m2(jk,:) = Planets(iPlanet).F_serp_bd_W_m2(jk,:);
%                        Planets(iPlanet).F_serp_bd_recharge_W_m2(jk,add_inds)=Planets(iPlanet).F_serp_bd_recharge_W_m2(jk,add_inds) + F_serp_bd_recharge_W_m2;
%                        FH2_bd_recharge_molecules_cm2_s = sum(Planets(iPlanet).FH2_bd_molecules_cm2_s(jk,t_inds))*(sum(t_ya(t_inds)))/(t_ya(t_inds(1))-t_ya(t_inds(1)-1))/2;
%            
%                       Planets(iPlanet).FH2_bd_recharge_molecules_cm2_s(jk,:) = Planets(iPlanet).FH2_bd_molecules_cm2_s(jk,:);
%                        Planets(iPlanet).FH2_bd_recharge_molecules_cm2_s(jk,add_inds) = Planets(iPlanet).FH2_bd_molecules_cm2_s(jk,add_inds)+FH2_bd_recharge_molecules_cm2_s;                       
%                        
%                        Planets(iPlanet).FH2_bd_recharge_moles_yr(jk,:) = Planets(iPlanet).FH2_bd_moles_yr(jk,:);
%                        Planets(iPlanet).FH2_bd_recharge_moles_yr(jk,add_inds) = Planets(iPlanet).FH2_bd_moles_yr(jk,add_inds)+FH2_bd_recharge_molecules_cm2_s.*3.15569e7./6.022e23.*4*pi*(1e2*Planets(iPlanet).R_m).^2;                  
%                    end
%                end
%              end            
%         end
		if find(Planets(iPlanet).z_cracking_1mm_m>0)
            for jk = 1:length(Planets(iPlanet).d_ocean_km)
                Planets(iPlanet).surf_area_seafloor(jk) = 4*pi*(Planets(iPlanet).R_m-Planets(iPlanet).d_ocean_km(jk)*1e3)^2; 
                [F_serp_W_m2,FH2_molecules_cm2_s,~] = getHeatFlux(Planets(iPlanet).R_m,Planets(iPlanet).d_ocean_km(jk)*1e3,Planets(iPlanet).z_cracking_1mm_m(jk,:),DERVS);
                %FS = -diff(F_serp_W_m2)./diff(Planets(iPlanet).t_yr)*4.5e9;H2 = -diff(Planets(iPlanet).FH2_molecules_cm2_s)./diff(Planets(iPlanet).t_yr)*4.5e9;
                %Planets(iPlanet).F_serp_W_m2 = [FS FS(end)];Planets(iPlanet).FH2_molecules_cm2_s = [H2 H2(end)];

                if strcmp(Planets(iPlanet).name,'Enceladus') || strcmp(Planets(iPlanet).name,'Iapetus')
                    Planets(iPlanet).F_serp_W_m2(jk,:) = F_serp_W_m2;
                    Planets(iPlanet).mean_F_serp_W_m2(jk) = mean(F_serp_W_m2);
                    Planets(iPlanet).FH2_molecules_cm2_s(jk,:) = FH2_molecules_cm2_s;
                    Planets(iPlanet).mean_FH2_molecules_cm2_s(jk) = mean(FH2_molecules_cm2_s);                    
                else
                    if DERVS
                        disp(iPlanet)
                        if iPlanet==10
                            x = 1;
                        end
                            pp = spline(Planets(iPlanet).t_yr*yr2s,F_serp_W_m2*Age_s); % this is a kluge to convert heat production averaged over the age of the solar system back to energy of reaction for the depth of cracking. In the calculations at the bottom of the page I divide by the age of the solar system, so this just gets back to the right units.
                        Planets(iPlanet).F_serp_W_m2(jk,:) = -pp1derv(pp,Planets(iPlanet).t_yr*yr2s); % negative is for reverse direction in time
                            pp = spline(Planets(iPlanet).t_yr*yr2s,FH2_molecules_cm2_s*Age_s); % another kluge, see comment two lines up.
                        Planets(iPlanet).FH2_molecules_cm2_s(jk,:) = -pp1derv(pp,Planets(iPlanet).t_yr*yr2s);
                    else
                        Planets(iPlanet).F_serp_W_m2(iPlanet,:) = F_serp_W_m2;
                         Planets(iPlanet).FH2_molecules_cm2_s(jk,:) = FH2_molecules_cm2_s;
                    end
                    Planets(iPlanet).mean_F_serp_W_m2(jk) = mean(F_serp_W_m2);
                    Planets(iPlanet).mean_FH2_molecules_cm2_s(jk) = mean(FH2_molecules_cm2_s);
                    Planets(iPlanet).FH2_moles_yr(jk,:) = Planets(iPlanet).FH2_molecules_cm2_s(jk,:).*3.15569e7./6.022e23.*4*pi*(1e2*Planets(iPlanet).R_m).^2; % previously used value (in 2007 paper: 3.1536e7)
                end
                Planets(iPlanet).F_rad_W_m2(jk,:) = Planets(iPlanet).H_W/Planets(iPlanet).surf_area_seafloor(jk);
                
               if isfield(Planets(iPlanet),'t_recharge_Mya')
                   if ~isempty(Planets(iPlanet).t_recharge_Mya)
                       t_ya = Planets(iPlanet).t_yr;
                       t_inds = find(t_ya>Planets(iPlanet).t_recharge_Mya*1e6);
                       add_inds = [t_inds(1)-2 t_inds(1)-1]; % the time of recharge is just a few steps, so is model dependent. The current setting (for the 2016 paper) supposes that Europa reserpentinizes over 100 Myr.
                       F_serp_recharge_W_m2 = sum(Planets(iPlanet).F_serp_W_m2(jk,t_inds))*(sum(t_ya(t_inds)))/(t_ya(t_inds(1))-t_ya(t_inds(1)-1))/2;
                       Planets(iPlanet).F_serp_recharge_W_m2(jk,:) = Planets(iPlanet).F_serp_W_m2(jk,:);
                       Planets(iPlanet).F_serp_recharge_W_m2(jk,add_inds)=Planets(iPlanet).F_serp_recharge_W_m2(jk,add_inds) + F_serp_recharge_W_m2;
                       FH2_recharge_molecules_cm2_s = sum(Planets(iPlanet).FH2_molecules_cm2_s(jk,t_inds))*(sum(t_ya(t_inds)))/(t_ya(t_inds(1))-t_ya(t_inds(1)-1))/2;
           
                      Planets(iPlanet).FH2_recharge_molecules_cm2_s(jk,:) = Planets(iPlanet).FH2_molecules_cm2_s(jk,:);
                       Planets(iPlanet).FH2_recharge_molecules_cm2_s(jk,add_inds) = Planets(iPlanet).FH2_molecules_cm2_s(jk,add_inds)+FH2_recharge_molecules_cm2_s;                       
                       
                       Planets(iPlanet).FH2_recharge_moles_yr(jk,:) = Planets(iPlanet).FH2_moles_yr(jk,:);
                       Planets(iPlanet).FH2_recharge_moles_yr(jk,add_inds) = Planets(iPlanet).FH2_moles_yr(jk,add_inds)+FH2_recharge_molecules_cm2_s.*3.15569e7./6.022e23.*4*pi*(1e2*Planets(iPlanet).R_m).^2;                  
                   end
               end                       
            end
		end
	end


figure(155);clf; set(gcf,'Name','Dist Serp Heat Flux with Time, z_cracking_1mm');
    set(gca,'YAxisLocation','right','YScale','log','YLim',[0 600],'XDir','reverse','FontSize',CaptSize);hold on;box on;
    for iPlanet = 1:lPlanet
		if find(Planets(iPlanet).z_cracking_1mm_m>0)
            line(Planets(iPlanet).t_yr*yr2Gyr,Planets(iPlanet).F_serp_W_m2*W2mW,'Color',Planets(iPlanet).plot_color,'LineWidth',LW); 
            text(Planets(iPlanet).t_yr(end)*yr2Gyr,Planets(iPlanet).F_serp_W_m2(end)*W2mW,Planets(iPlanet).name,'FontSize',CaptSize);
            h=line(Planets(iPlanet).t_yr*yr2Gyr,Planets(iPlanet).F_rad_W_m2*W2mW);set(h,'Color',Planets(iPlanet).plot_color,'LineWidth',LW+2);
            text(Planets(iPlanet).t_yr(end)*yr2Gyr,Planets(iPlanet).F_rad_W_m2(end)*W2mW,[Planets(iPlanet).name ': Radiogenic'],'FontSize',CaptSize);
            if isfield(Planets(iPlanet),'t_recharge_Mya')
                       if ~isempty(Planets(iPlanet).t_recharge_Mya)
                                   h=line(Planets(iPlanet).t_yr*yr2Gyr,Planets(iPlanet).F_serp_recharge_W_m2*W2mW);set(h,'Color',Planets(iPlanet).plot_color,'LineWidth',LW+2,'LineStyle','--');                               
                       end
            end
		end
    end
%     MOR_F_mW_m2 =   MOR_F_W_m2*W2mW;
%     h = line(Planets(1).t_yr*yr2Gyr,MOR_F_mW_m2*ones(1,length(Planets(1).t_yr)));set(h,'Color','g','LineWidth',LW+2);
    xlabel('Time (Gya)','FontSize',CaptSize);
    ylabel('Heat (mW m^{-2})','FontSize',CaptSize);h

% figure(1555);clf; set(gcf,'Name','Dist Serp Heat Flux with Time, brittle-ductile from Fournier 1999');
%     set(gca,'YAxisLocation','right','YScale','log','YLim',[0 600],'XDir','reverse','FontSize',CaptSize);hold on;box on;
%     for iPlanet = 1:lPlanet
%             if ~isempty(Planets(iPlanet).F_serp_bd_W_m2)
%                 for jk = 1:length(Planets(iPlanet).F_serp_bd_W_m2(:,1))
%                     line(Planets(iPlanet).t_yr*yr2Gyr,squeeze(Planets(iPlanet).F_serp_bd_W_m2(jk,:))*W2mW,'Color',Planets(iPlanet).plot_color,'LineWidth',LW); 
%                     text(Planets(iPlanet).t_yr(end)*yr2Gyr,squeeze(Planets(iPlanet).F_serp_bd_W_m2(jk,end))*W2mW,[Planets(iPlanet).name 'bd'],'FontSize',CaptSize);
%                     if ~isempty(Planets(iPlanet).F_rad_W_m2)
%                         h=line(Planets(iPlanet).t_yr*yr2Gyr,Planets(iPlanet).F_rad_W_m2*W2mW);set(h,'Color',Planets(iPlanet).plot_color,'LineWidth',LW+2);
%                         text(Planets(iPlanet).t_yr(end)*yr2Gyr,Planets(iPlanet).F_rad_W_m2(end)*W2mW,[Planets(iPlanet).name ': Radiogenic'],'FontSize',CaptSize);
%                     end
%                     if isfield(Planets(iPlanet),'t_recharge_Mya') && ~isempty(Planets(iPlanet).F_serp_bd_recharge_W_m2)
%                                if ~isempty(Planets(iPlanet).t_recharge_Mya)
%                                            h=line(Planets(iPlanet).t_yr*yr2Gyr,squeeze(Planets(iPlanet).F_serp_bd_recharge_W_m2(jk,:))*W2mW);set(h,'Color',Planets(iPlanet).plot_color,'LineWidth',LW+2,'LineStyle','--');                               
%                                end
%                     end
%             end
%         end
%     end
% %     MOR_F_mW_m2 =   MOR_F_W_m2*W2mW;
% %     h = line(Planets(1).t_yr*yr2Gyr,MOR_F_mW_m2*ones(1,length(Planets(1).t_yr)));set(h,'Color','g','LineWidth',LW+2);
%     xlabel('Time (Gya)','FontSize',CaptSize);
%     ylabel('Heat (mW m^{-2})','FontSize',CaptSize);



%     ax1 = gca;
% 	ax2 = axes('Position',get(ax1,'Position'),...
% 			   'YAxisLocation','right','YColor','r',...
% 			   'XAxisLocation','top','XColor','r',...
% 			   'Color','none');
% 		   xc = 0.5;
figure(156);clf;
    for iPlanet = 1:lPlanet		
		if find(Planets(iPlanet).z_cracking_1mm_m>0)
        h=line(Planets(iPlanet).t_yr*yr2Gyr, Planets(iPlanet).FH2_molecules_cm2_s);set(h,'Color',Planets(iPlanet).plot_color,'LineWidth',LW);
        text(Planets(iPlanet).t_yr(end)*yr2Gyr,Planets(iPlanet).FH2_molecules_cm2_s(end),Planets(iPlanet).name,'FontSize',CaptSize);
                if isfield(Planets(iPlanet),'t_recharge_Mya')
                   if ~isempty(Planets(iPlanet).t_recharge_Mya)
                        h=line(Planets(iPlanet).t_yr*yr2Gyr, Planets(iPlanet).FH2_recharge_molecules_cm2_s);set(h,'Color',Planets(iPlanet).plot_color,'LineWidth',LW,'LineStyle','--');
                   end
        end
		end    
    end
%         h = line(Planets(1).t_yr*yr2Gyr,MOR_FH2_molecules_cm2_s*ones(1,length(Planets(1).t_yr)));set(h,'Color','g','LineWidth',LW+2);

	% xlimits = get(ax1,'XLim');
	% ylimits = get(ax1,'YLim');
	% xinc = (xlimits(2)-xlimits(1))/5;
	% yinc = (ylimits(2)-ylimits(1))/5;
% 	xlim = [0 170];
% 	set(ax1,'YLim',[1e-3 1e2/3],'XLim',xlim,'YScale','log','FontSize',CaptSize);
 	set(gca,'YAxisLocation','right','YScale','log','XDir','reverse','FontSize',CaptSize);
    
	ylabel('H_2 Production (molecules cm^{-2} s^{-1})','FontSize',CaptSize,...
		'Color','k');
      xlabel('Time (Gya)','FontSize',CaptSize);
	axis tight;
	box on;
% 	gtext('H_2','FontSize',CaptSize)   
    
    figure(157);clf;
    for iPlanet = 1:lPlanet		
		if ~isempty( Planets(iPlanet).FH2_moles_yr)% plot present day values
            h=line(Planets(iPlanet).t_yr*yr2Gyr, Planets(iPlanet).FH2_moles_yr);set(h,'Color',Planets(iPlanet).plot_color,'LineWidth',LW);
            text(Planets(iPlanet).t_yr(end)*yr2Gyr,Planets(iPlanet).FH2_moles_yr(end),Planets(iPlanet).name,'FontSize',CaptSize);
%             if ~isempty(Planets(iPlanet).FH2_bd_moles_yr)
%                 for jk = 1:length(Planets(iPlanet).FH2_bd_moles_yr(:,1))
%                     h=line(Planets(iPlanet).t_yr*yr2Gyr, Planets(iPlanet).FH2_bd_moles_yr(jk,:));set(h,'Color',Planets(iPlanet).plot_color,'LineWidth',LW+1);
%                     text(Planets(iPlanet).t_yr(end)*yr2Gyr,Planets(iPlanet).FH2_bd_moles_yr(jk,end),[Planets(iPlanet).name ' bd'],'FontSize',CaptSize);
%                 end
%             end
            if isfield(Planets(iPlanet),'t_recharge_Mya') 
               if ~isempty(Planets(iPlanet).t_recharge_Mya)
                    h=line(Planets(iPlanet).t_yr*yr2Gyr, Planets(iPlanet).FH2_recharge_moles_yr);set(h,'Color',Planets(iPlanet).plot_color,'LineWidth',LW,'LineStyle','--');
%                     if ~isempty(Planets(iPlanet).FH2_bd_moles_yr)
%                         for jk = 1:length(Planets(iPlanet).FH2_bd_moles_yr(:,1))
%                                 h=line(Planets(iPlanet).t_yr*yr2Gyr, Planets(iPlanet).FH2_bd_recharge_moles_yr(jk,:));set(h,'Color',Planets(iPlanet).plot_color,'LineWidth',LW+1,'LineStyle','--');
%                         end
%                     end
               end
           end
        end
		end    
%         h = line(Planets(1).t_yr*yr2Gyr,MOR_FH2_moles_yr*ones(1,length(Planets(1).t_yr)));set(h,'Color','g','LineWidth',LW+2);

	% xlimits = get(ax1,'XLim');
	% ylimits = get(ax1,'YLim');
	% xinc = (xlimits(2)-xlimits(1))/5;
	% yinc = (ylimits(2)-ylimits(1))/5;
% 	xlim = [0 170];
% 	set(ax1,'YLim',[1e-3 1e2/3],'XLim',xlim,'YScale','log','FontSize',CaptSize);
 	set(gca,'YAxisLocation','right','YScale','log','XDir','reverse','FontSize',CaptSize);
    
	ylabel('H_2 Production (moles yr^{-1})','FontSize',CaptSize,...
		'Color','k');
      xlabel('Time (Gya)','FontSize',CaptSize);
	axis tight;
	box on;
%     	R_Earth = 6371e3;
% 	
% 	R_Europa = 1565e3;
% 	A_Europa = 4*pi*R_Europa^2;
% 	d_Europa = 100*1e3; %m
% 	Dh_Europa = 0:1e2:20e3; % m
% 	Dh_Europa_km = Dh_Europa*m2km;
% 
% 	R_Mars = 3397e3;
% 	A_Mars = 4*pi*R_Mars^2;
% 	d_Mars = 0;
% 	Dh_Mars = 0:100:12e3;%m
% 	Dh_Mars_km = Dh_Mars*m2km;
% 
% 	R_Enceladus = 249e3;
% 	A_Enceladus = 4*pi*R_Enceladus^2;
% 	d_Enceladus = 80*1e3;
% 	Dh_Enceladus = 0:100:169e3;
% 	Dh_Enceladus_km = Dh_Enceladus*m2km;
% 	[J_Earth_1mm,JH2_Earth_1mm] = getHeatFlux(R_Earth,Earth.d_ocean*m2km,Earth.z_cracking_1mm_m*m2km);
% 	[J_Europa_1mm,JH2_Europa_1mm] = getHeatFlux(R_Europa,Europa.d_ocean*m2km,Europa.z_cracking_1mm_m*m2km);
% 	[J_Mars_1mm,JH2_Mars_1mm] = getHeatFlux(R_Mars,Mars.d_ocean*m2km,Mars.z_cracking_1mm_m*m2km);
% 	[J_Enceladus_1mm,JH2_Enceladus_1mm] = getHeatFlux(R_Enceladus,80,169);
% 
% 	figure(155); clf; set(gcf,'Name','Dist Serp Heat Flux with Time, z_cracking_1mm');
% 	
% 	line(Earth.t_yr,[-diff(J_Earth_1mm)./diff(Earth.t_yr) 0]*4.5e9*W2mW,'Color','g','LineWidth',LW); 
%     hold on
% 	line(Europa.t_yr,[-diff(J_Europa_1mm)./diff(Europa.t_yr) 0]*4.5e9*W2mW,'Color','b','LineWidth',LW); 
% 	h = hline(J_Enceladus_1mm*W2mW);set(h,'Color','c','LineWidth',LW);
% 	line(Mars.t_yr,[-diff(J_Mars_1mm)./diff(Mars.t_yr) 0]*4.5e9*W2mW,'Color','r','LineWidth',LW+1);
% % 		line(Mars.t_yr,[-diff(J_Mars_1mm)./diff(Mars.t_yr) J_Mars_1mm(end)]*4.5e9*W2mW,'Color','r','LineWidth',LW+1);
% 	Jradiothermal_Earth = Earth.H_W/4/pi/(R_Earth-Earth.d_ocean).^2;
% 		h=line(Earth.t,Jradiothermal_Earth*W2mW);set(h,'Color','g','LineWidth',LW+5)
% 	Jradiothermal_Europa = Europa.H_W/4/pi/(R_Europa-Europa.d_ocean).^2;
% 		h=line(Europa.t,Jradiothermal_Europa*W2mW);set(h,'Color','b','LineWidth',LW+5)
% 	Jradiothermal_Mars = Mars.H_W/4/pi/(R_Mars-Mars.d_ocean).^2;
% 		h=line(Mars.t,Jradiothermal_Mars*W2mW);set(h,'Color','r','LineWidth',LW+5);
% 	%Jtidal_Europa = 20;h=hline(Jtidal_Europa);set(h,'Color','b','LineWidth',LW);
% 	Jradiothermal_Enceladus = Enceladus.H_W/4/pi/(R_Enceladus-Enceladus.d_ocean).^2;
% 		h=line(Enceladus.t,Jradiothermal_Enceladus*W2mW);set(h,'Color','c','LineWidth',LW+5);
% 
% 	text(Europa.t(end)/1e9,max(J_Europa_1mm)*1e3,{'Enceladus'},'FontSize',CaptSize);
% 	text(Europa.t(end)/1e9,max(J_Europa_1mm)*1e3/2.02/3,{'Europa: Radiogenic Heat'},'FontSize',CaptSize);
% 	text(Europa.t(end)/1e9,20/3,{'Europa: Radiogenic + Tidal Heat'},'FontSize',CaptSize);
% 	text(Europa.t(end)/1e9,20/3,{'Enceladus: Radiogenic'},'FontSize',CaptSize);
% 	ax1 = gca;
% 	xlim = [0 4.5e9];
% 		xlabel('Time Before Present (ya)','FontSize',HeadSize);
% 	ylabel('Heat Flux (mW m^{-2})','FontSize',HeadSize); 
% 	set(ax1,'XLim',xlim,'YScale','log','XScale','log','FontSize',CaptSize);
% 	grid on
% 	
% %%
% 	figure(156);clf
% 	ax2 = axes('YAxisLocation','right','YColor','r',...
% 				'XDir','reverse',...
% 			   'XAxisLocation','bottom','XColor','r' );
% 		   xc = 0.5;
% 		  % line(Earth.t,[-diff(JH2_Earth_1mm)./diff(Earth.t) 0]*4.5e9,'Color','g','LineWidth',LW);
% 		  % line(Europa.t,[-diff(JH2_Europa_1mm)./diff(Europa.t) 0]*4.5e9,'Color',[1 0 xc],'LineWidth',LW);
% 		   line(Mars.t,[-diff(JH2_Mars_1mm)./diff(Mars.t) 0]*4.5e9,'Color','r','LineWidth',LW+1); hold on
% 		   hline((JH2_Enceladus_1mm));set(h,'Color',[1 xc xc],'LineWidth',LW);
% 	% xlimits = get(ax1,'XLim');
% 	% ylimits = get(ax1,'YLim');
% 	% xinc = (xlimits(2)-xlimits(1))/5;
% 	% yinc = (ylimits(2)-ylimits(1))/5;
% 	set(ax2,'XLim',xlim,'YScale','log','XScale','log','FontSize',CaptSize);
% 	ylabel('H_2 Production (molecules cm^{-2} s^{-1})','FontSize',HeadSize,...
% 		'Color','k','Rotation',-90)
% 		xlabel('Time Before Present (ya)','FontSize',HeadSize);
% 	gtext('H_2','FontSize',CaptSize)  
% 	grid on
else
	R_Earth = 6371e3;
	
	R_Europa = 1565e3;
	A_Europa = 4*pi*R_Europa^2;
	d_Europa = 100*1e3; %m
	Dh_Europa = 0:1e2:20e3; % m
	Dh_Europa_km = Dh_Europa*m2km;

	R_Mars = 3397e3;
	A_Mars = 4*pi*R_Mars^2;
	d_Mars = 0;
	Dh_Mars = 0:100:12e3;%m
	Dh_Mars_km = Dh_Mars*m2km;

	R_Enceladus = 249e3;
	A_Enceladus = 4*pi*R_Enceladus^2;
	d_Enceladus = 80*1e3;
	Dh_Enceladus = 0:100:169e3;
	Dh_Enceladus_km = Dh_Enceladus*m2km;
	[J_Europa,JH2_Europa] = getHeatFlux(R_Europa,d_Europa,Dh_Europa,DERVS);

	[J_Mars,JH2_Mars] = getHeatFlux(R_Mars,d_Mars,Dh_Mars,DERVS);


	[J_Enceladus,JH2_Enceladus] = getHeatFlux(R_Enceladus,d_Enceladus,Dh_Enceladus,DERVS);
	
	figure(154);clf; set(gcf,'Name','Dist Serp Heat Flux');clf
	line(Dh_Europa_km,J_Europa*W2mW,'Color','b','LineWidth',LW); hold on
	line(Dh_Enceladus_km,J_Enceladus*W2mW,'Color','c','LineWidth',LW);
	line(Dh_Mars_km,J_Mars*W2mW,'Color','r','LineWidth',LW+1);
	Jradiothermal_Europa = 9;h=hline(Jradiothermal_Europa);set(h,'Color','b','LineWidth',LW)
	Jtidal_Europa = 20;h=hline(Jtidal_Europa);set(h,'Color','b','LineWidth',LW);
	Jradiothermal_Enceladus = 0.005;h=hline(Jradiothermal_Enceladus);set(h,'Color','c','LineWidth',LW);
	xlabel('Depth of Hydration (km)','FontSize',HeadSize);
	ylabel('Heat Flux (mW m^{-2})','FontSize',HeadSize); 
	text(Dh_Enceladus_km(end)/3,max(J_Europa)*1e3,{'Enceladus'},'FontSize',CaptSize);
	text(Dh_Enceladus_km(end)/4,max(J_Europa)*1e3/2.02/3,{'Europa: Radiogenic Heat'},'FontSize',CaptSize);
	text(Dh_Enceladus_km(end)/4,20/3,{'Europa: Radiogenic + Tidal Heat'},'FontSize',CaptSize);
	text(Dh_Enceladus_km(end)/4,20/3,{'Enceladus: Radiogenic'},'FontSize',CaptSize);
	ax1 = gca;
	ax2 = axes('Position',get(ax1,'Position'),...
			   'YAxisLocation','right','YColor','r',...
			   'XAxisLocation','top','XColor','r',...
			   'Color','none');
		   xc = 0.5;
		   line(Dh_Europa_km,JH2_Europa,'Color',[1 0 xc],'LineWidth',LW);
		   line(Dh_Mars_km,JH2_Mars,'Color','r','LineWidth',LW+1);
		   line(Dh_Enceladus_km,JH2_Enceladus,'Color',[1 xc xc],'LineWidth',LW);
	% xlimits = get(ax1,'XLim');
	% ylimits = get(ax1,'YLim');
	% xinc = (xlimits(2)-xlimits(1))/5;
	% yinc = (ylimits(2)-ylimits(1))/5;
	xlim = [0 170];
	set(ax1,'YLim',[1e-3 1e2/3],'XLim',xlim,'YScale','log','FontSize',CaptSize);
	set(ax2,'YLim',[1e7 1e12/3],'XLim',xlim,'YScale','log','FontSize',CaptSize);
	ylabel('H_2 Production (molecules cm^{-2} s^{-1})','FontSize',HeadSize,...
		'Color','k','Rotation',-90)
	gtext('H_2','FontSize',CaptSize)   
end
% %%
% figure(153); set(gcf,'Name','Serp Heat Flux');clf 
% %plot(D_hydr,P_olpyr*1e3,D_hydr,P_ol*1e3,D_hydr,0.1*ones(length(D_hydr)),'g
% %',D_hydr,2*ones(length(D_hydr)),'g','LineWidth',LW)
% plot(Dh_Europa_km,J_Europa*A_Europa,'b','LineWidth',LW); hold on
% plot(Dh_Enceladus_km,J_Enceladus*A_Enceladus,'c','LineWidth',LW);
% plot(Dh_Mars_km,J_Mars*A_Mars,'r','LineWidth',LW);
% title('Heat Flux from Serpentinization','FontSize',HeadSize,'FontWeight','bold');
% xlabel('Depth of Hydration (km)','FontSize',HeadSize,'FontWeight','bold');
% ylabel('Heat Flux (W)','FontSize',HeadSize,'FontWeight','bold');
% %set(gca,'YScale','log',
% set(gca,'FontSize',CaptSize);
% axis tight
% %gtext('Olivine + Pyroxene','FontSize',CaptSize); 


% figure(145);set(gcf,'Name','H2 Production');
% plot(Dh_Europa_km,JH2_Europa,'b',Dh_Mars_km,JH2_Mars,'r',Dh_Enceladus_km,JH2_Enceladus,'c','LineWidth',LW);
% set(gca,'YScale','log','FontSize',CaptSize);
% ylabel('H_2 Production (molecules cm^{-2} s^{-1})','FontSize',CaptSize,'FontWeight','bold')
% xlabel('Depth of Serpentinization (km)','Fontsize',CaptSize,'FontWeight','bold');

%======================================================================
function [F_W_m2,FH2_molecules_cm2_s,mH2O_kg] = getHeatFlux(R_m,d_m,z_m,DERVS)

% weights of elements
W_Mg = 24.305; % g/mol
W_Si = 28.086;
W_O = 15.999;
W_H = 1.0079;

% weights of molecules
%W_pyr = W_Mg + W_Si +3*W_O; % pyroxene: MgSiO3
%W_ol = 2*W_Mg + W_Si + 4*W_O; % olivine: Mg2SiO4
%W_h2o = 2*W_H + W_O;
W_serp = 3*W_Mg + 2*W_Si + 9*W_O + 4*W_H;
W_bruc = W_Mg + 2*W_O + 2*W_H;

% volumes of reaction
%V_pyr = 31.47; % cm3/mol
V_ol = 43.79;
%V_h2o = 1; 
V_serp = 108.5;
V_bruc = 24.63; 

% heats of reaction (converted from kcal to J)
% H_hydr_olpyr = 6.7e4 %16*1000*4.18400; % heat of serp reaction hydrating olivine and pyroxene, in J per mol
H_hydr_ol_J_mol = 41e3; %19.5*1000*4.18400/2;

N_Avogadro = 6.0221415e23;
Age_s = 4.5*pi*3.15576e7*1e9; %age of planet, in seconds
R_seafloor_m = R_m-d_m; % m
A_seafloor_m2 = 4*pi*R_seafloor_m.^2; % m2
Rh_m = R_seafloor_m-z_m; % m
Vh_m3 = ones(1,length(Rh_m));
if DERVS
Vh_m3 = 4/3*pi*(R_seafloor_m.^3-Rh_m.^3); % m3
else
    Vh_m3(end) = 4/3*pi*(R_seafloor_m.^3-Rh_m(end).^3); % m3
    for iPlanet = 1:length(Rh_m)-1
        Vh_m3(iPlanet) = 4/3*pi*(Rh_m(iPlanet+1).^3-Rh_m(iPlanet).^3); % m3
    end
end

cm3m3=1e6;
% n_serp_ol = Vh_m3./(V_serp+V_bruc)*(W_serp/(W_serp+W_bruc))*cm3m3;

Frac_per = 1.0; %Fraction of peridotite in crust
Frac_oliv= 0.7; %Fraction of olivine in peridotite
Frac_fay = 0.1; %Frac_fay(Fe/Mg)=Fraction of fayalite in olivine
Frac_serp= 1.0; %Fraction of fayalite that is serpentinizable (Schulte et al. 2006 use Frac_serp=0.5)
n_hydr_ol = Frac_per* Frac_oliv*Frac_serp* Vh_m3 / (V_ol/cm3m3); % serpentinizable olivine [moles] in crust
n_hydr_fa = Frac_fay * n_hydr_ol; % Serpentinizable fayalite [moles] in crust
mH2O_kg = (2*W_H+W_O)*(2*Frac_fay+3*(1-Frac_fay))*Frac_per*Frac_oliv*Vh_m3/V_ol*cm3m3/1e3;

m2tocm2 = 1e4;
molecules_cm2_per_mol_m2 = N_Avogadro/m2tocm2;

%nH2 = 2.3*M_hydr_ol; % mmol from Allen 2003, Table 2, pp. 1535, 100% olivine
nH2 = (2/3)*n_hydr_fa; % Moles of H2 production from 3 fa + 2 h2o = 2 mag + 3 sio2 + 2 h2
FH2_molecules_cm2_s = nH2./A_seafloor_m2/Age_s*molecules_cm2_per_mol_m2;
% From Zahnle et al 2006
% Table 1: Volcanic fluxes. All quantities given in units of molecules cm?2s?
% Model               Title   SO2         H2S         H2          CO
% V1 - ?modern low?           1 ? 10^9    1 ? 10^8    2 ? 10^9    2 ? 10^8
% V2 - ?modern high?          3.5 ? 10^9  3.5 ? 10^8  1 ? 10^10   1 ? 10^9
% V3 - ?Archean high?         1 ? 10^10   1 ? 10^9    3 ? 10^10   3 ? 10^9

% H_olpyr = n_serp_olpyr*H_hydr_olpyr;
H_ol_J = n_hydr_ol*H_hydr_ol_J_mol; % changed to reflect more accurate serpentinization of both fa and fo

% P_olpyr = H_olpyr/A_seafloor/Age;
F_W_m2 = H_ol_J./A_seafloor_m2/Age_s;
