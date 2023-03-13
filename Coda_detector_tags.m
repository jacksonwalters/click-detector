function [TOA_tag,TOA_other,Coda_Type_on_Whale,Coda_Type_off_Whale]=Coda_detector_tags(F_ds,Y_filtered,Plot_flag,W_seg,SNR_window,SNR_thresh)

% SNR_window=SNR_window_echo;
% SNR_thresh=SNR_thresh_echo;
% consistency_T=consistency_T_echo;
% T_sec=10;
% T=T_sec*F_ds;
% NOI=floor(File_duration/T_sec);

Buffer_ind=1;
   TOA_tag={};
   TOA_other={};
Coda_Type_on_Whale={};
Coda_Type_off_Whale={};

%  for Buffer_ind=1:NOI
     J1=[]; J2=[]; Locs1=[]; Locs2=[];
     ratios=[]; ratios2=[];
     ind_on_Whale=[]; ind_off_Whale=[];

%        Y_filtered=bandpass(Y_decimated(int32((Buffer_ind-1)*T+1):int32((Buffer_ind-1)*T+T)),[F_low, F_high],F_ds);     % Aply band pass filter and extract buffer                             
       
        c=0;                   % counter for the transient detection stage
        Y_zoom=Y_filtered;     % Aply band pass filter               
        t_zoom=[0:1/F_ds:(1/F_ds)*(length(Y_zoom)-1)];                                              % Time calls of the analyzed signal [sec]
        [ey,~]=energyop(Y_zoom,0);              % Apply TKEO (ey is the output signal with enhanced SNR)
        ey_norm=ey/max(ey);                     % Normalize the enhaced signal ey          
        ED_thresh=0.5*raylfit(ey_norm);           % Define threshold of minimum allowed noise (in the enhanced signal)
        ED_thresh=1e-4;
        time=[0:1/F_ds:(1/F_ds)*(length(ey_norm)-1)];   
        [pks,locs] =findpeaks(ey_norm,F_ds,'MinPeakHeight',ED_thresh,'MinPeakDistance',25e-3);  % Apply instantaneous energy detector (find peaks)          
        locs_samples=locs*F_ds;                 % store indices of time readings of the detected impulses

        %% Eliminate transients bellow the minimum allowed SNR
        
        Locs=[]; Pks=[];
        for i=1:length(locs)
            if locs_samples(i)>SNR_window && (locs_samples(i)+SNR_window)<length(ey_norm)
                tmp=ey_norm(int32(locs_samples(i)-SNR_window):int32(locs_samples(i)+SNR_window));
                SNR(i)=pks(i)/median(tmp); 
                if SNR(i)>SNR_thresh
                    c=c+1;
                    Locs(c)=locs(i);
                    Pks(c)=pks(i);
                end
            end
        end


        if length(Pks)>100
           [Pks2,I] = maxk(Pks,100);
           Locs2=Locs(I);
           Pks=Pks2;
           Locs=Locs2;
        end

        
        
 %% validation:
  
        
 
    seg_ds=round(W_seg*F_ds); 
    locs_samples=Locs*F_ds; 
    time=[0:1/F_ds:(1/F_ds)*(length(ey_norm)-1)];
        
    BANK=[];  M_bank=[];
        slice=0.0005*F_ds;
%         slice=0.002*F_ds;
    for ij=1:length(locs_samples)
%           BANK(ij,:)=Y_zoom(int32(locs_samples(ij)-slice):int32(locs_samples(ij)+slice)).^2;
          BANK(ij,:)=ey_norm(int32(locs_samples(ij)-slice):int32(locs_samples(ij)+slice)).^2;
          M_bank(ij)=sum(BANK(ij,:));
    end
    
    [PKS,Inds]=sort(M_bank);

    idx = kmedoids(PKS',2);

    check=unique(idx);
    

    
%    if length(check)==2 
       
    [PKS,Inds]=sort(M_bank);    
%     other_whale=Inds(find(idx==1));
%     on_whale=Inds(find(idx==2));
%     
%           
%  Locs1=Locs(on_whale);    Locs2=Locs(other_whale);
%  Pks_tmp1=Pks(on_whale);  Pks_tmp2=Pks(other_whale); 
 
  
   
   x=[1:3];
   
   for i=1:length(PKS)-2
       seq=PKS(i:2+i);       
       P = polyfit(x,seq,1);
       yfit = polyval(P,x);
       m(i)=yfit(2)-yfit(1);
   end
   

   
   
D=abs(diff(m));
segment_ind=0;

for i=1:length(D)
    if mean(m(1:i))<0.3 && D(i)>0.5
        segment_ind=i;
        break;
    end    
end

%      figure; 
% %    subplot(4,1,1);
%    plot([1:segment_ind+2],PKS(1:segment_ind+2),'o'); hold on;
%    plot([segment_ind+3:length(PKS)],PKS(segment_ind+3:end),'o');
% %    subplot(4,1,2); plot([1:length(idx)],idx,'o')
% %    subplot(4,1,3); plot([1:length(m)],m,'o')
% %    subplot(4,1,4); plot([1:length(m)-1],abs(diff(m)),'o')
%    

if segment_ind>0

    other_whale=Inds(segment_ind+3:length(PKS));
    on_whale=Inds(1:segment_ind+2);
  Locs1=Locs(on_whale);    Locs2=Locs(other_whale);
  Pks_tmp1=Pks(on_whale);  Pks_tmp2=Pks(other_whale); 
    
          
 Locs1=Locs(on_whale);    Locs2=Locs(other_whale);
 Pks_tmp1=Pks(on_whale);  Pks_tmp2=Pks(other_whale); 
 
else
 
    Locs1=Locs;    Locs2=Locs;
    Pks_tmp1=Pks;  Pks_tmp2=Pks;
    
end
    
%     figure; 
%     subplot(3,1,1); plot(Y_zoom); xlabel('time [sec]'); ylabel('Amplitude'); title('Raw signal'); 
%     subplot(3,1,2); plot(time,ey_norm); hold on;
%    plot(Locs1,Pks_tmp1,'ro'); plot(Locs2,Pks_tmp2,'go'); hold off;
 
 roi=F_ds*25e-3;
 t_roi=F_ds*Locs2;
 for i=1:length(t_roi)
 Zoom_roi=Y_zoom(t_roi(i)-roi:t_roi(i)+roi);
 t_zoom=[0:1/F_ds:(1/F_ds)*(length(Zoom_roi)-1)];
 
 [p_zoom,l_zoom] =findpeaks(Zoom_roi,F_ds,'MinPeakHeight',ED_thresh,'MinPeakDistance',2e-3);  % Apply instantaneous energy detector (find peaks)          
% figure; plot(1e3*t_zoom,Zoom_roi); hold on; plot(1e3*l_zoom,p_zoom,'x');
[Ip,It]=sort(p_zoom,'descend');
choose_p1(:,i)=Ip(1:6);
choose_t1(:,i)=l_zoom(It(1:6));
% J1(i)=JF(choose_t1(:,i));
J1(i)=var(choose_t1(:,i));
ratios(:,i)=choose_p1(1:5,i)./choose_p1(2:6,i);

 end
 
 t_roi=F_ds*Locs1;
 for i=1:length(t_roi)
 Zoom_roi=Y_zoom(t_roi(i)-roi:t_roi(i)+roi);
 t_zoom=[0:1/F_ds:(1/F_ds)*(length(Zoom_roi)-1)];
 
 [p_zoom,l_zoom] =findpeaks(Zoom_roi,F_ds,'MinPeakHeight',ED_thresh,'MinPeakDistance',2e-3);  % Apply instantaneous energy detector (find peaks)          
% figure; plot(1e3*t_zoom,Zoom_roi); hold on; plot(1e3*l_zoom,p_zoom,'x');
[Ip,It]=sort(p_zoom,'descend');
choose_p2(:,i)=Ip(1:6);
choose_t2(:,i)=l_zoom(It(1:6));
% J2(i)=JF(choose_t2(:,i));
J2(i)=var(choose_t2(:,i));
ratios2(:,i)=choose_p2(1:5,i)./choose_p2(2:6,i);

 end

% [coeff,score,latent,tsquared,explained] = pca(ratios);

% figure; 
% subplot(3,1,1); plot(choose_p1,'rx-'); hold on; plot(choose_p2,'gx-');
% subplot(3,1,2); plot(J1,'rx-');  hold on; plot(J2,'gx-');
% subplot(3,1,3); plot(ratios,'rx-');  hold on; plot(ratios2,'gx-');

% figure; 
%  plot(J1'*ones(1,size(ratios,1)),ratios','ro');  hold on; plot(J2'*ones(1,size(ratios2,1)),ratios2','go'); grid on;
 
J=[J1 J2];
R1=ratios'; R2=ratios2';
% AV1=mean(R1,1); AV2=mean(R2,1);
load AV1; load AV2;
MSE1=sum((R1-ones(size(R1,1),1)*AV1).^2,2)';
MSE2=sum((R2-ones(size(R2,1),1)*AV2).^2,2)';
%  figure; 
%  plot(J',[MSE1 MSE2],'ko','Linewidth',2);  grid on; hold on;
 
 
 Voting_on_whale(Buffer_ind)=sum(J1<5e-5 & MSE1<2.5)/length(J1);
 Voting_off_whale(Buffer_ind)=sum(J2<5e-5 & MSE2<2.5)/length(J2);
      
 
 
% 
%   subplot(3,1,3)
%   plot(J1',MSE1,'ro','Linewidth',2);  grid on; hold on; 
%   plot(J2',MSE2,'go','Linewidth',2);  grid on; hold off;
  
  
%  end
 
%

 if (Voting_on_whale(Buffer_ind)>0.6 || Voting_off_whale(Buffer_ind)>0.6)

    [on_whale_seq,I1]=sort(Locs1);  
    [other_whale_seq,I2]=sort(Locs2);           
    Pks1=Pks_tmp1(I1);
    Pks2=Pks_tmp2(I2);


     D_on_Whale=diff(on_whale_seq);
     D_off_Whale=diff(other_whale_seq);

     lag_thresh=0.55;
     Segment_on_Whale=find(D_on_Whale>lag_thresh);
     Segment_off_Whale=find(D_off_Whale>lag_thresh);

     Segment_on_Whale=[0 Segment_on_Whale];
     Segment_off_Whale=[0 Segment_off_Whale];

      for i=1:length(Segment_on_Whale)
         if i<length(Segment_on_Whale)
             ind_on_Whale(i,:)=[Segment_on_Whale(i)+1 Segment_on_Whale(i+1)-1];
             Coda_Type_on_Whale(i)={Coda_word(D_on_Whale(Segment_on_Whale(i)+1:Segment_on_Whale(i+1)-1)')};
         elseif Segment_on_Whale(i)+1<length(D_on_Whale)
             ind_on_Whale(i,:)=[Segment_on_Whale(i)+1 length(D_on_Whale)];
             Coda_Type_on_Whale(i)={Coda_word(D_on_Whale(Segment_on_Whale(i)+1:end)')};
         end
     end

     for i=1:length(Segment_off_Whale)
         if i<length(Segment_off_Whale)
             ind_off_Whale(i,:)=[Segment_off_Whale(i)+1 Segment_off_Whale(i+1)-1];
             Coda_Type_off_Whale(i)={Coda_word(D_off_Whale(Segment_off_Whale(i)+1:Segment_off_Whale(i+1)-1)')};
         elseif Segment_off_Whale(i)+1<length(D_off_Whale)
             ind_off_Whale(i,:)=[Segment_off_Whale(i)+1 length(D_off_Whale)];
             Coda_Type_off_Whale(i)={Coda_word(D_off_Whale(Segment_off_Whale(i)+1:end)')};
         end
     end


     if Plot_flag==1
        figure;set(gcf, 'Position', get(0,'Screensize'));
        t_zoom=[0:1/F_ds:(1/F_ds)*(length(Y_zoom)-1)];  
        subplot(4,1,1); plot(t_zoom,Y_zoom); xlabel('time [sec]'); ylabel('Amplitude'); title('Raw signal');
        subplot(4,1,2); plot(time,ey_norm); 
        if segment_ind>0
        hold on; plot(Locs(other_whale),Pks(other_whale),'rx','Linewidth',2)
        hold on; plot(Locs(on_whale),Pks(on_whale),'gx','Linewidth',2)           
        xlabel('time [sec]'); ylabel('TKEO'); ylim([0 1]); 
        legend('','whale 2','whale 1');
        else
        hold on; plot(Locs,Pks,'gx','Linewidth',2)           
        xlabel('time [sec]'); ylabel('TKEO'); ylim([0 1]); 
        legend('','whale 1');
        end
        subplot(4,1,3); plot(time,ey_norm);
        legendInfo{1} = [''];
        for j=1:size(ind_on_Whale,1)
            hold on; plot(on_whale_seq(ind_on_Whale(j,1):ind_on_Whale(j,2)+1),Pks1(ind_on_Whale(j,1):ind_on_Whale(j,2)+1),'o','Linewidth',2)     
            TOA_tag(j)={on_whale_seq(ind_on_Whale(j,1):ind_on_Whale(j,2)+1)};
            xlabel('time [sec]'); ylabel('TKEO'); ylim([0 1]); title('whale 1','Fontsize',14);
            legendInfo{j+1} = ['Coda-' num2str(j) ': ' cell2mat(Coda_Type_on_Whale(j))];
        end
            legend(legendInfo)

            if segment_ind>0
                 subplot(4,1,4); plot(time,ey_norm);
                 legendInfo2{1} = [''];
                for j=1:size(ind_off_Whale,1)
                    hold on; plot(other_whale_seq(ind_off_Whale(j,1):ind_off_Whale(j,2)+1),Pks2(ind_off_Whale(j,1):ind_off_Whale(j,2)+1),'o','Linewidth',2)       
                    TOA_other(j)={other_whale_seq(ind_off_Whale(j,1):ind_off_Whale(j,2)+1)};
                    xlabel('time [sec]'); ylabel('TKEO'); ylim([0 1]); title('whale 2','Fontsize',14); 
                    legendInfo2{j+1} = ['Coda-' num2str(j) ': ' cell2mat(Coda_Type_off_Whale(j))];
                end
                legend(legendInfo2)
            end
     end            

 else
    if Plot_flag==1
        figure;set(gcf, 'Position', get(0,'Screensize'));
        t_zoom=[0:1/F_ds:(1/F_ds)*(length(Y_zoom)-1)];  
        subplot(2,1,1); plot(t_zoom,Y_zoom); xlabel('time [sec]'); ylabel('Amplitude'); title('Raw signal');
        subplot(2,1,2); plot(time,ey_norm);       
     end

 end
 
%  end
 
% %    tab={'type','TOA'};
% for i=1:size(on_whale_seq,1)
%    tab={cell2mat(Coda_Type_on_Whale(i)),num2str(round(on_whale_seq(ind_on_Whale(i,1):ind_on_Whale(i,2)+1),2))};
%    writecell(tab,['On-whale.xls'],'WriteMode','append');
% end
%    writecell(tab,['Other-whale.xls'],'WriteMode','append');
% 
% 
%    
%    TOA_tag(i)={on_whale_seq(ind_on_Whale(i,1):ind_on_Whale(i,2)+1)};
%    TOA_other(j)=other_whale_seq(ind_off_Whale(i,1):ind_off_Whale(i,2)+1);
   
end 
%     
        