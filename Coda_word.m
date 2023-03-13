
function Coda_Type=Coda_word(Detected_pattern)

% Detected_pattern=D_on_Whale(Segment_on_Whale(i)+1:end)';

% data = readcell('Codas_Reference2.xlsx');
%     
% random_letter = 'NOISE';
% word_list = data(:,end);
% has_r = not( cellfun( @isempty, regexp( word_list, random_letter ) ) );
% word_list( has_r ) = [];
% Data=data(~has_r,:);
%   
%     D=Data(:,1);
%     Analysis=Data(:,2:11);
% 
%     ICI=zeros(size(Data,1),9);
%     for i=1:size(Analysis,1) 
%         NOC(i)=cell2mat(D(i,1));
%         for j=1:size(Analysis,2)      
%             ICI(i,j)=cell2mat(Analysis(i,j));       
%         end    
%     end

load ICI; load Analysis; load Data; load NOC;

    tokeep = Detected_pattern ~= 0;   %logical array indicating which elements to keep
    workingP = Detected_pattern(tokeep);
    D_NOC = length(workingP);
    A_inds=find(NOC==D_NOC+1);

    Correspondance=[]; count=0;
    for i=A_inds 
        count=count+1;
        ICI_cand(count,:)=ICI(i,1:D_NOC);
        word(count)=Data(i,end);
        Correspondance(count)=mean((ICI_cand(count,:)-Detected_pattern').^2);
    end

    Thresh=0.03;
    xxx=sort(Correspondance);
    
    check=min(Correspondance);
    if check<Thresh
       id=find(Correspondance==check);
    end
    
    CT=Data(id,end);
    
    Coda_Type=[];
    if ~isempty(CT)
        Coda_Type=cell2mat(CT);
    else
        Coda_Type=['Unseen'];
    end 
    
end