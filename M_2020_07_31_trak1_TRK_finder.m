clear
%You need the IlluminaHiSeq data of gene expression and the Curated survival data to run this
%code
%Input your cancer file names here in this vertical manner, also labeled in
%the same format
cancerList = [
%     "Data_Acute_Myeloid_Leukemia_HiSeqV2" %Problematic, no primary tumor site
    "Data_Adrenocortical_HiSeqV2"
    "Data_Bile_Duct_HiSeqV2"
    "Data_Bladder_HiSeqV2"
    "Data_Breast_HiSeqV2"
    "Data_Cervical_HiSeqV2"
    
    "Data_Colon_HiSeqV2"
    "Data_Endometrioid_HiSeqV2"
    "Data_Esophageal_HiSeqV2"
    "Data_Glioblastoma_HiSeqV2"
    "Data_Head_and_Neck_HiSeqV2"
    "Data_Kidney_Chromophobe_HiSeqV2"
    "Data_Kidney_Clear_Cell_HiSeqV2"
    "Data_Kidney_Papillary_Cell_HiSeqV2"
    "Data_Large_B-cell_Lymphoma_HiSeqV2"
    "Data_Liver_HiSeqV2"
    
    "Data_Lower_Grade_Glioma_HiSeqV2"
    "Data_Lung_Adenocarcinoma_HiSeqV2"
   
    "Data_Lung_Squamous_Cell_HiSeqV2"
    "Data_Melanoma_HiSeqV2"
    "Data_Mesothelioma_HiSeqV2"
    "Data_Ocular_Melanomas_HiSeqV2"
    "Data_Ovarian_HiSeqV2"
    "Data_Pancreatic_HiSeqV2"
    "Data_Pheochromocytoma_Paraganglioma_HiSeqV2"
    "Data_Prostate_HiSeqV2"
    "Data_Rectal_HiSeqV2"
    "Data_Sarcoma_HiSeqV2"
    "Data_Stomach_HiSeqV2"
    "Data_Testicular_HiSeqV2"
    "Data_Thymoma_HiSeqV2"
    "Data_Thyroid_HiSeqV2"
    "Data_Uterine_HiSeqV2"
    ];
% cancerList=[
%     "Data_Pheochromocytoma___Paraganglioma_HiSeqV2"
%     ];
geneOfInterest = [
   "COL1A1"
   "LMNB1"
   "COL1A2"
   "MKI67"
    ]; 

%You can change how you want your Trk and trak1 files to be ultimately
%named here if you want to
nameMod1 = erase(cancerList, "Data");
nameMod1 = erase(nameMod1, "HiSeqV2");


for y = 1:size(geneOfInterest,1)

    trkNames = strcat("Trk", nameMod1, geneOfInterest(y));
    trak1Names = strcat("trak1", nameMod1, geneOfInterest(y));
    
    for z = 1:size(cancerList,1)

        DATA=tdfread(cancerList(z),'\t');
        %fprintf("data read");
        AvPro1=struct2table(DATA);
        c=fieldnames(DATA);
        c=c(2:length(c),1);
        c=char(c);
        siz=size(AvPro1);
        k=1;



        for i=1:siz(2)-1
            %fprintf("Primary Tumor hits %d \n", i);
                if (c(i,23:24)== '01')
                 %main(:,k)=power(2,AvPro1{:,i+1}); 
                 main(:,k)=AvPro1{:,i+1};   
                 k=k+1;
                end

        end

        %Storing the number of patients with primary tumors for each cancer
        %type
        numOfPatientsWithPrimTum(z,:) = [cancerList(z), size(main,2)];

        k=1;
        for i=1:siz(1)
            %fprintf("Storing gene names with more than 75% %d \n\n", i);
           if(sum(main(i,:) == 0)<0.75*siz(2))
               main1(k,:)=main(i,:);
               b(k,:)=AvPro1{i,1}; %storing gene names of patients who have more than 75% of reads
               k=k+1;
           end
        end

        siz=size(main1);
        siz1=size(b);


        for i=1:siz(1)
            %fprintf("Finding gene id to scale with %d\n", i);
            if(string(strtrim(b(i,:)))== geneOfInterest(y))
                protein=i;
                break;
            end
        end

        Trk=zeros(siz(1),3);

        for i=1:siz(1)
            %fprintf("Running robust linear regression %d %d %d\n",y,z,i);
            mdf=fitlm(main1(protein,:),main1(i,:),'linear','RobustOpts','on');
            Trk(i,1)=mdf.Coefficients.Estimate(2); 
            Trk(i,2)=mdf.Coefficients.Estimate(1);
            Trk(i,3)=mdf.Rsquared.Ordinary;
        end

        % genes that scale strongly with the gene of interest
        j=1;
        trak1=cell(5,5);
        for i=1:siz(1)
            %fprintf("genes that scale strongly with gene of interest %d\n", i);
            if(Trk(i,3)>0.5 & Trk(i,1)>0.5)
                trak1{j,1}=Trk(i,1);
                trak1{j,2}=Trk(i,2);
                trak1{j,3}=Trk(i,3);
                trak1{j,4}=i;
                trak1{j,5}=strtrim(b(i,:));
                j=j+1;
            end
        end

        writematrix(Trk, trkNames(z));
        writecell(trak1, trak1Names(z));

    %     scatter(Trk(:,3),Trk(:,1),30,'b','filled');

        clear main
        clear main1
        clear Trk
        clear trak1
        clear b
        clear c

    end
end


