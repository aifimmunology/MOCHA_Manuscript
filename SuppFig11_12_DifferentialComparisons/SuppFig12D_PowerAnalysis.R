################################################################################
## Do a power analysis based on a 
## 2-sample T-test with different 
## sample sizes to showcase power loss
## via downsampling 
################################################################################

### Original N=39
###   - 17 COVID+
###   - 22 COVID- 

prop1 = 22/39
prop2 = 17/39

### Set effect size to 0.5 (medium cohen's effect size)
effectSize = 0.5

### Dowsample from N=39 to N=30
### and keep the # of samples per group
### proportional to original cohort 

power_es05 <- rbindlist(lapply(39:30,
			 function(x){
			 	n1 = round(prop1 * x)
			 	n2 = round(prop2 * x)
			 	
			 	data.frame(
			 		N1= n1,
			 		N2= n2,
			 		TotalN=x,
			 		Power=
			 	
			 	pwr::pwr.t2n.test(n1 = n1,
			 										n2 = n2,
			 										d=effectSize,
			 										sig.level = 0.05
			 										)$power			 	
			 	)
			 }
			 	
			 	)
)

### Calculate Relative Power based on 
### How much power there is at each N
### relative to N=39
power_es05$RelativePower = power_es05$Power / power_es05[TotalN==39]$Power

### create plot and save results 
pdf('Power_Analysis_SuppFig12D.pdf')
ggplot(power_es05,
			 aes(x=TotalN,
			 		 y= RelativePower))+geom_point()+geom_smooth()+
	theme_minimal() + ggtitle("Relative Power of 2-Sample T-test After Downsampling")+
	xlab('Total Sample Size')+
	ylab('Relative Power')+
	theme(text=element_text(size=20))
dev.off()
			 
write.csv(power_es05,
					file='../PowerAnalysis.csv')


