#####################################################################
#####################################################################
require(ggplot2)
require(scMACS)
require(data.table)
require(dplyr)

cell_sampled <- c(5, 10, 25, 50, 100,
                      200, 500,
                      1000, 2000,
                      3000, 5000,
                      10000, 15000, 20000, 30000, 40000, 50000, 75000,
                      100000,110000, 120000, 130000,
                      140000)



interpolateCoefficients <- function(cell_model){
    
    ### Apply model fit based on the
    ### cell abundance. First Part is a loess
    ### Fit, final part is a linear fit,
    ### Middle part is an average between
    ### Loess & Linear fits.

    if (cell_model <= 100000) {

      ## Loess Fit
      Intercept <- stats::predict(
        finalModelObject$Loess$Intercept,
        data.frame(NumCells = cell_model)
      )
      Total <- stats::predict(
        finalModelObject$Loess$Total,
        data.frame(NumCells = cell_model)
      )
      Max <- stats::predict(
        finalModelObject$Loess$Max,
        data.frame(NumCells = cell_model)
      )
      tmpModel <- c(Intercept, Total, Max)
      tmpModel <- as.matrix(tmpModel)
    } else if (cell_model > 100000 & cell_model <= 140000) {

      ## Average of Loess & Linear Fit
      Intercept <- mean(
        stats::predict(
          finalModelObject$Loess$Intercept,
          data.frame(NumCells = cell_model)
        ),
        stats::predict(
          finalModelObject$Linear$Intercept,
          data.frame(NumCells = cell_model)
        )
      )

      Total <- mean(
        stats::predict(
          finalModelObject$Loess$Total,
          data.frame(NumCells = cell_model)
        ),
        stats::predict(
          finalModelObject$Linear$Total,
          data.frame(NumCells = cell_model)
        )
      )
      Max <- mean(
        stats::predict(
          finalModelObject$Loess$Max,
          data.frame(NumCells = cell_model)
        ),
        stats::predict(
          finalModelObject$Linear$Max,
          data.frame(NumCells = cell_model)
        )
      )

      tmpModel <- c(Intercept, Total, Max)
      tmpModel <- as.matrix(tmpModel)
    } else {
      ## Linear Fit
      Intercept <- stats::predict(finalModelObject$Linear$Intercept, data.frame(NumCells = cell_model))
      Total <- stats::predict(finalModelObject$Linear$Total, data.frame(NumCells = cell_model))
      Max <- stats::predict(finalModelObject$Linear$Max, data.frame(NumCells = cell_model))
      tmpModel <- c(Intercept, Total, Max)
      tmpModel <- as.matrix(tmpModel)
    }
    newdata <- data.frame(Ncells = cell_model)
    threshold <-  as.numeric(stats::predict(scMACS::youden_threshold, newdata = newdata))
    tmpModel <- c(tmpModel, threshold)
    return(t(tmpModel))
}                    


### Get Model Values 
model_interpolations <- do.call(rbind, lapply(cell_sampled,
       function(x) 
       cbind(interpolateCoefficients(x),x)
            )
                                         )

colnames(model_interpolations) <- c('Intercept','L1','L2','Threshold','Cells')
model_interpolations = as.data.table(model_interpolations)

setwd('/home/jupyter/MOCHA_Manuscript/Fig1/')
write.csv(model_interpolations,
          file='final_model.csv',
          row.names=F)

model_interpolations <- melt(model_interpolations, measure.vars=c('Intercept','L1','L2','Threshold'))

#### 

pdf('finalModel.pdf')
ggplot(model_interpolations[variable!='Threshold'],
       aes(x=Cells,
           y=value))+geom_line()+geom_point()+facet_wrap(~variable, ncol=1, scales='free_y')+theme_minimal()+
    theme(text=element_text(size=20))+
    xlab('Cells Sampled') + ylab('Coefficient Value')
dev.off()

pdf('threshold.pdf')
ggplot(model_interpolations[variable=='Threshold'],
       aes(x=Cells,
           y=value))+geom_line()+geom_point()+facet_wrap(~variable, ncol=1, scales='free_y')+theme_minimal()+
    theme(text=element_text(size=20))+
    xlab('Cells Sampled') + ylab('Threshold Value')+
    ylim(0,1)
dev.off()