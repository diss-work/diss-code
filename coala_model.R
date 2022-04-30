library(coala)
library(ggplot2)



model1=coal_model(sample_size=c(0,1607,0,0,0,1615,0,0,0) ,loci_number=1000, loci_length=16000)+ 
  
  #Sampling 1607 and 1615 haploid from populations 2 and 6
  
  feat_mutation(rate=1.25e-8*4*10000*16000) + 
  feat_recombination(rate=1.0e-8*4*10000*16000)+ #Rate=locus_length * per_generation_per_base_pair_mutation_rate * 4 * Ne
  
  feat_migration(0.25*4*10000,pop_from = 3, pop_to = 2 , time = 0) +  #Create 'inland' and 'coastal' 
  feat_migration(0.25*4*10000,pop_from = 4, pop_to = 2 , time = 0) +  #By replacing with 25% from each village
  feat_migration(0.25*4*10000,pop_from = 5, pop_to = 2 , time = 0) +  
  feat_migration(0.25*4*10000,pop_from = 7, pop_to = 6 , time = 0) +
  feat_migration(0.25*4*10000,pop_from = 8, pop_to = 6 , time = 0) +
  feat_migration(0.25*4*10000,pop_from = 9, pop_to = 6 , time = 0) +
  
  feat_size_change(((1452/8)/10000)*0.75, p="all", time=1/40000)+  #Set village pop sizes 
  
  feat_size_change((800/10000)*0.75, p=1, time=1/40000)+                #Pop. 1: 'Melanesian' ancestry
  feat_migration(rate=0.015*4*10000,symmetric=TRUE, time=1/40000)+      #Set migration rate
 
  #*0.75 for all pop sizes as simualting G6PD on X chrom                                                                
  
  feat_migration(0,pop_from=2,pop_to=1,time=1/40000)+                 #Immediately remove migration from others to Pop 1
  feat_migration(0,pop_from=3,pop_to=1,time=1/40000)+                 
  feat_migration(0,pop_from=4,pop_to=1,time=1/40000)+
  feat_migration(0,pop_from=5,pop_to=1,time=1/40000)+
  feat_migration(0,pop_from=6,pop_to=1,time=1/40000)+
  feat_migration(0,pop_from=7,pop_to=1,time=1/40000)+
  feat_migration(0,pop_from=8,pop_to=1,time=1/40000)+
  feat_migration(0,pop_from=9,pop_to=1,time=1/40000)+
  
  feat_migration(0,pop_from=1,pop_to=2,time=1/40000)+                 #Immediately remove migration from Pop 1 to others
  feat_migration(0,pop_from=1,pop_to=3,time=1/40000)+                 
  feat_migration(0,pop_from=1,pop_to=4,time=1/40000)+
  feat_migration(0,pop_from=1,pop_to=5,time=1/40000)+
  feat_migration(0,pop_from=1,pop_to=6,time=1/40000)+
  feat_migration(0,pop_from=1,pop_to=7,time=1/40000)+
  feat_migration(0,pop_from=1,pop_to=8,time=1/40000)+
  feat_migration(0,pop_from=1,pop_to=9,time=1/40000)+
  
  feat_pop_merge(time=161/40000,pop_source=3, pop_target = 2)+      #merge all villages into pop 2 at 160 gens
  feat_pop_merge(time=161/40000,pop_source=4, pop_target = 2)+      
  feat_pop_merge(time=161/40000,pop_source=5, pop_target = 2)+
  feat_pop_merge(time=161/40000,pop_source=6, pop_target = 2)+
  feat_pop_merge(time=161/40000,pop_source=7, pop_target = 2)+
  feat_pop_merge(time=161/40000,pop_source=8, pop_target = 2)+
  feat_pop_merge(time=161/40000,pop_source=9, pop_target = 2)+
  
  feat_size_change(((1425/10000))*0.75, population=2, time=161/40000)+ #Change pop 2 size at 160 gens 
  
  feat_migration(rate=0.26*4*10000, pop_from=1,pop_to =2,time = 162/40000)+  #pop 2 receives 26% from 'Melanesian' pop 1
  feat_migration(rate=0.0,pop_from=1,pop_to=2, time = 163/40000)+            
  #^^ pulse of admixture from pop 1 to 2
  
  feat_size_change((2050/10000)*0.75, population=2, time=164/40000)+     
  
  feat_pop_merge(time=2002/40000, pop_source=2,pop_target=1)+                #Merge pop 2 into 1 at 2000 gens 
  feat_size_change((10500/10000)*0.75, 1, time = 2002/40000)+                #ancestral pop. size at 2002 gens
  sumstat_jsfs("jsfs_2_6", populations=c(2,6), per_locus=FALSE) 

model1


stats=simulate(model1)

nrow(stats$jsfs_2_6)#1608 - inland
stats$jsfs_2_6[47,] #47th row corresponds to variants at freq. 46 in the inland pop

stats$jsfs_2_6[42:52,] #41st to 51st rows of JSFS

slice=as.data.frame(stats$jsfs_2_6[42:52,])
summary(slice)
ncol(slice)#1616 - coastal
sumofcols=colSums(slice) #sum of columns 

#P-val
sum(sumofcols[92:1616]) / sum(sumofcols) #col 92 will correspond to variants at freq. of 91
#proportion of overall SNPs [(sum(sumofcols))] that are above a freq. of 91 
#0.2690355 (for this run) - not significant, unable to detect non-neutrality


#Test plot
index=seq_along(sumofcols)
a=data.frame(sumofcols,index)
summary(a)
ggplot(a, aes(index, sumofcols))+geom_point()+labs(y="Summed frequency")

min_sum=min(sumofcols)
max_sum=max(sumofcols)
lims=c(min_sum,max_sum)


#Cumulative plot

cumsumofcols=cumsum(sumofcols)
cbind(a,cumsumofcols)
cumsumofcols #finding where it levels off - 394
cumsumofcols<394 
ggplot(a, aes(index,cumsumofcols))+
  geom_point(size=1, shape=16, colour="black")+
  geom_vline(xintercept=91, colour="red", linetype="dashed")+
  labs(x="Number of derived mutations at a site", y="Cumulative frequency")+
  annotate("rect",xmin=92, xmax=344, ymin=0, ymax=Inf, alpha=0.4, fill="lightblue")
#ggsave("Cumsum_of_JSFS.tiff",width=5, height=5, units="in",dpi=300)
#ggsave("Cumsum_of_JSFS_shading.tiff",width=5, height=5, units="in",dpi=300)
dev.off()


#Cumulative proportion plot - Figure used in results

cumsumofcols.p=(cumsum(sumofcols))/400 #proportion 
cbind(a,cumsumofcols.p)
cumsumofcols #finding where it levels off - 394
cumsumofcols<394 
ggplot(a, aes(index,cumsumofcols.p))+
  geom_point(size=1, shape=16, colour="black")+
  geom_vline(xintercept=91, colour="red", linetype="dashed")+
  labs(x="< X derived alleles in coastal population group", y="Cumulative proportion of SNPs")+
  annotate("rect",xmin=92, xmax=344, ymin=0, ymax=Inf, alpha=0.4, fill="lightblue")
#ggsave("Cumsum_of_JSFS.tiff",width=5, height=5, units="in",dpi=300)
#ggsave("Cumsum_of_JSFS_shading_proportion.tiff",width=5, height=5, units="in",dpi=300)
dev.off()


