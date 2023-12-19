require(data.table)
require(ggplot2)
dats = fread('cd16_allTiles.csv')
tested_dats = dats[Avg_Intensity_Control > 12 | Avg_Intensity_Case >12| abs(Pct0_Case - Pct0_Control) > 0.5,]

pdf('~/Google Drive/Shared drives/Imm - CompBio - Temp Backup/Samir_workspace/Project2: MOCHA scMACS/Figure 3/bi_modal_pvalues.pdf',
		width=11, height=7)
ggplot(tested_dats,
			 aes(x=P_value,
			 		 ))+geom_histogram(bins=50)+
		theme_minimal()+
		xlab('P values')
dev.off()

write.csv(tested_dats,
					file='~/Google Drive/Shared drives/Imm - CompBio - Temp Backup/Samir_workspace/Project2: MOCHA scMACS/Figure 3/bi_modal_pvals.csv')
