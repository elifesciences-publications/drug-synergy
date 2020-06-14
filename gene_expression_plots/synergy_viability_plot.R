plot_synergy_viability_v2<-function(viability_data,platform,
                                    pval_diff,
                                    data,method,coordin,gene_set,nrandom,ngenes)
{
  
  conds=c("TM","TW","MW")
  time=c(3,6,9,12,24)
  
  require(gplots)
  require(ggplot2)
  if(coordin=="t_score")
  {
    coord=3  
  }
  
  if(coordin=="LFC")
  {
    coord=1 
  }
  ndim=10
  
  E_bliss=vector()
  cond=vector()
  per_syn=vector()
  n_syn=vector()
  n_non_syn=vector()
  combo=vector()
  combo_angle=vector()
  combo_angle_n=vector()
  cor_de=vector()
  sd_cor_de=vector()
  sd_combo_angle=vector()
  sd_combo_angle_n=vector()
  Treatment=vector()
  l=1
  
  
  for(i in conds)
  {
    cond1=substr(i,1,1)
    cond2=substr(i,2,2)
    cond3=i
    # datax=data#data[,c(paste0(cond1,"_",time),paste0(cond2,"_",time),paste0(cond3,"_",time))]
    # fittt <- prcomp(t(datax),
    #                 center = TRUE,
    #                 scale. = TRUE) 
    
    par(mfrow=c(2,3))
    
    for(j in time)
    {
      
      genes_de_ae=vector()
      
      d1=get(paste0("limma_eb_",cond1,"_",j))  
      d2=get(paste0("limma_eb_",cond2,"_",j)) 
      d3=get(paste0("limma_eb_",cond3,"_",j))  
      
      
      condss1=paste0(cond1,"_",j)
      condss2=paste0(cond2,"_",j)
      condss3=paste0(cond3,"_",j)
      
      
      diff1=rownames(d1)[which(d1[,5]<pval_diff)]
      diff2=rownames(d2)[which(d2[,5]<pval_diff)]
      diff3=rownames(d3)[which(d3[,5]<pval_diff)]
      
      
      ax=list(diff1,diff2,diff3)
      names(ax)=c(condss1,condss2,condss3)
      venn(ax)
      title(paste0("DE Genes,t=",j))
      
      genes_de_ae=union(diff1,diff2)
      #genes_de_ae=sample(rownames(d1),50)
      #genes_de_ae=union(genes_de_ae,rownames(d3)[which(d3[,5]<pval_diff)] )
      
      
      
      
      
      
      genes_de_ae=intersect(genes_de_ae,gene_set)
      print(length(genes_de_ae))
      
      dummy_cor=vector()
      dummy_combo_angle=vector()
      dummy_combo_angle_n=vector()
      
      for(jj in 1:nrandom)
      {
        genes_r=sample(genes_de_ae,ngenes)
        d1_r=d1[genes_r,coord]
        d2_r=d2[genes_r,coord]
        d3_r=d3[genes_r,coord]
        dummy_cor[jj]=cor(d1_r,d2_r,method=method)
        dummy_combo_angle[jj]=acos(sum(d1_r*d2_r)/(sqrt(sum(d1_r^2))*sqrt(sum(d2_r^2))))
        dummy_combo_angle[jj]=dummy_combo_angle[jj]*(180/pi)
        xls=lsfit(cbind(d1_r,d2_r),d3_r,intercept=FALSE)$residuals
        xls=d3_r-xls
        dummy_combo_angle_n[jj]=acos(sum(d3_r*xls)/(sqrt(sum(d3_r^2))*sqrt(sum(xls^2))))
        dummy_combo_angle_n[jj]=dummy_combo_angle_n[jj]*(180/pi)
      } 
      
      cor_de[l]=mean(dummy_cor)
      sd_cor_de[l]=sd(dummy_cor)
      
      #d1_an=d1[genes_de_ae,coord]
      #d2_an=d2[genes_de_ae,coord]
      #d3_an=d3[genes_de_ae,coord]
      
      combo_angle[l]=mean(dummy_combo_angle)
      sd_combo_angle[l]=sd(dummy_combo_angle)
      
      #d12=as.vector(fittt$x[condss3,1:ndim])
      
      combo_angle_n[l]=mean(dummy_combo_angle_n)
      sd_combo_angle_n[l]=sd(dummy_combo_angle_n)
      
      condd=paste0(i,"_",j)
      combo[l]=i
      cond[l]=condd
      E_bliss[l]=viability_data[intersect(which(drug_viability[,"Assay"]==platform),
                                          which(drug_viability[,"DRUGS"]==condd)),"Excess_over_bliss"]
      
      Treatment[l]=i
      n_syn[l]=length(intersect(gene_set,setdiff(diff3,union(diff1,diff2))))
      n_non_syn[l]=length(intersect(gene_set,diff3))#-length(setdiff(diff3,union(diff1,diff2)))
      per_syn[l]=n_syn[l]/max(n_syn[l],n_non_syn[l])
      
      l=l+1
    }
    
    par(mfrow=c(1,1))
    
  }
  
  abc=cond
  # abc=gsub("T","A",abc)
  # abc=gsub("M","B",abc)
  # abc=gsub("W","C",abc)
  
  print(combo_angle)
  n_non_syn=log(n_non_syn)
  data=data.frame(condition=abc,n_synergy=log(n_syn),percentage_syn=per_syn,
                  ebliss=E_bliss,combo=combo,corr=cor_de,non_syn=n_non_syn,
                  angle=combo_angle,angle_n=combo_angle_n,
                  sd_cor=sd_cor_de,sd_angle=sd_combo_angle,
                  Treatment=Treatment,
                  sd_angle_n=sd_combo_angle_n)
  
  cor_nsyn=cor(log(n_syn),E_bliss,method=method)
  cor_nper=cor(per_syn,E_bliss,method=method)
  cor_non_syn=cor(n_non_syn,E_bliss,method=method)
  cor_an=cor(combo_angle,E_bliss,method=method)
  cor_an_n=cor(combo_angle_n,E_bliss,method=method)
  cor_dee=cor(cor_de,E_bliss,method=method)
  
  print("eren")
  sizen=3
  height=4
  
  
  View(data)
  
  p=ggplot(data, aes(n_synergy, ebliss, label =data$condition ))+ 
    #geom_text(aes(colour = factor(combo)),size=sizen)+
    geom_point(aes(colour = Treatment),size=sizen)+
    ggtitle(paste0("Number of Synergistic Genes vs EOB (Platform=",platform,", Correlation=",cor_nsyn," )"))+
    theme(axis.text=element_text(size=12,face="bold"))+
    theme(title=element_text(size=10,face="bold")) +
    xlab("Number of Synergistic Genes (log)")+
    ylab("EOB (Excess Over Bliss)")+
    scale_color_manual(values=colors_mcf7[c(5,6,7),1]) 
  plot(p)
  
  
  p=ggplot(data, aes(non_syn, ebliss, label =data$condition ))+ geom_text(aes(colour = factor(combo)),size=sizen)+
    ggtitle(paste0("Number of Non-Synergistic Genes vs EOB Platform=",platform," Correlation=",cor_non_syn))
  #plot(p)
  
  p=ggplot(data, aes(percentage_syn, ebliss, label =data$condition ))+ 
    #geom_text(aes(colour = factor(combo)),size=sizen)+
    geom_point(aes(colour = Treatment),size=sizen)+
    ggtitle(paste0("Percentage of Synergistic Genes vs EOB (Platform=",platform,", Correlation=",cor_nper," )"))+
    theme(axis.text=element_text(size=12,face="bold"))+
    theme(title=element_text(size=10,face="bold"))+
    ylab("EOB (Excess Over Bliss)")+
    xlab("Percentage of Synerygistic Genes")+
    scale_color_manual(values=colors_mcf7[c(5,6,7),1])
  
  plot(p)
  
  p=ggplot(data, aes(angle, ebliss, label =data$condition ))+ 
    #geom_text(aes(colour = factor(combo)),size=sizen)+
    geom_point(aes(colour = Treatment),size=sizen)+
    ggtitle(paste0("Angle Between Monotherapies vs EOB (Platform=",platform,", Correlation=",cor_an," )"))+
    theme(axis.text=element_text(size=12,face="bold"))+
    theme(title=element_text(size=10,face="bold"))+
    geom_errorbarh(aes(xmax = angle + sd_angle, xmin = angle - sd_angle, height = height))+
    ylab("EOB (Excess Over Bliss)")+
    xlab("Angle (Degrees)")+
    scale_color_manual(values=colors_mcf7[c(5,6,7),1])
  plot(p)
  
  p=ggplot(data, aes(corr, ebliss, label =data$condition ))+ 
    #geom_text(aes(colour = factor(combo)),size=sizen)+
    geom_point(aes(colour = Treatment),size=sizen)+
    ggtitle(paste0("Correlation Between Monotherapies vs EOB (Platform=",platform,", Correlation=",cor_dee," )"))+
    theme(axis.text=element_text(size=12,face="bold"))+
    theme(title=element_text(size=10,face="bold")) +
    geom_errorbarh(aes(xmax = corr + sd_cor, xmin = corr - sd_cor, height = height))+
    ylab("EOB (Excess Over Bliss)")+
    xlab("Correlation")+
    scale_color_manual(values=colors_mcf7[c(5,6,7),1])
  plot(p)
  
  p=ggplot(data, aes(angle_n, ebliss, label =data$condition ))+ 
    #geom_text(aes(colour = factor(combo)),size=sizen)+
    geom_point(aes(colour = Treatment),size=sizen)+
    ggtitle(paste0("Angle Between Combo to Monotheraphy Plane vs EOB (Platform=",platform,", Correlation=",cor_an_n," )"))+
    theme(axis.text=element_text(size=12,face="bold"))+
    theme(title=element_text(size=10,face="bold")) +
    geom_errorbarh(aes(xmax = angle_n + sd_angle_n, xmin = angle_n - sd_angle_n, height = height))+
    ylab("EOB (Excess Over Bliss)")+
    xlab("Angle (Degrees)")+
    scale_color_manual(values=colors_mcf7[c(5,6,7),1])
  plot(p)
  
  return(data)
  
}  