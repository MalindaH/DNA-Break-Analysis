% for plotting aligned sequence peaks
clear all;
x=1;
y=0;
firstChr=1;
lastChr=22;

outputName='outputtest0826';
pq='q';
chrnum=10;

c=int2str(chrnum);

    filenameStem=strcat(outputName,'chr',c,'_',pq,'val.txt');
    STEM=load (filenameStem);
    hfig=figure;
    set(hfig,'Position',[0,0,2000,550]);
    %stem(STEM(:,1),STEM(:,2),'.');
    stem(STEM(:,1),-10*log10(STEM(:,2)),'.')
    ylim([0 50.00000])
    set(gca,'XTick',0)
    set(gca,'XTickLabel','')
    %ylabel('-10.log_1_0({\itq}-value)','fontsize',12);
    ylabel('-10log_{10}(pval)','fontsize',12)
    title(['\fontsize{13}chromosome ',c,' ',pq,'-value plot'])
