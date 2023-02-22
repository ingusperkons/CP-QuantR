
# Code (version R 4.1.1.)  
# Load libraries ####


library(readxl)
library(stringr)
library(stringi)
library(ggplot2)
library(tidyr)
library(dplyr)
library(knitr)
library(xlsx)
library(ggthemes)
library(broom)
library(ggpmisc)
library(gridExtra)


# Define variables ####

### Define file name for the XLSX export, where the results of the data analysis will be stored


temp.xlsx.filename = "Excel_data_report"


### Define file name of the exported peak integration data (.csv) from Skyline   
#'Note: the Skyline report must contain the following columns:  
#' "Product.Mz" - m/z values of CPs  
#' "Area" - peak areas of CPs  
#' "Replicate" - file names of the instrumental data files  

temp.csv.filename = "Skyline_export"


### Define the volume of the final sample extract (??L)  

Final.extract.volume = 100 #??L

### Define the molecular formula of the internal standard compound  

ISTD.formula = "C'10Cl6H16"


#===============================================================================
## Read external files  ####
#'All files must be placed in the same folder as the .Rmd file of the analysis 
#'To enable reading from the same folder as the .Rmd file is located:    
#'Session --> Set Working Directory --> To Source File Location*
  
### Read sequence and mass-list (both must be placed in the same folder of .R file)

LC.input <- read_xlsx("Reference_sheet.xlsx", sheet = 3)
LC.input = LC.input %>% mutate(ID.name = paste(molecular.formula,"-",round(m.z), sep = "")) %>% select(-m.z,-adduct)
LC.sequence.1M <- read_xlsx("Reference_sheet.xlsx", sheet = 1)

### Read .csv export from SkyLine

store.lc <- read.csv(paste(temp.csv.filename,".csv",sep = ""),sep = ",",header = T)
store.lc = store.lc %>% select(Precursor.Neutral.Formula,Product.Mz,Area,Replicate) %>%
  mutate(m.z = Product.Mz, peak.area = Area) %>% mutate(ID.name = paste(Precursor.Neutral.Formula,"-",round(m.z), sep = "")) %>% select(-Product.Mz,-Area)
store.lc = left_join(store.lc,LC.sequence.1M, by = "Replicate")
store.lc = left_join(store.lc,LC.input, by = "ID.name")
store.lc = store.lc %>%
  ungroup() %>%
  mutate(peak.area = ifelse(peak.area == "#N/A",NA,peak.area)) %>%
  mutate(peak.area = as.numeric(peak.area)) %>%
  group_by(Replicate) %>% mutate(ISTD.area = ifelse(molecular.formula == ISTD.formula,peak.area,0)) %>%
  mutate(ISTD.area = max(ISTD.area, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(molecular.formula = ifelse(molecular.formula == ISTD.formula,"ISTD",molecular.formula)) %>%
  filter(is.na(retention.time) != TRUE) %>%
  select(m.z,retention.time,row.identity,molecular.formula,peak.area,Sample.name, Sample.type, ISTD.area)

### Backup unprocessed input data

store.lc.1M <- store.lc

#===============================================================================

## Prepare data tables for quantification procedure ####
  
### Add carbon chain length and Cl number  
#'Note:  
#'Here only Cl numbers from 4 to 12 Cl atoms per CP homologue group are evaluated  
#'Line mutate(Cl.No = factor(Cl.No, levels = c("4","5","6","7","8","9","10","11","12"))) has to be changed to increase or decrease number of Cl atoms used in the following workflow.


store.lc.1M = store.lc.1M %>%
  ungroup() %>%
  mutate(Chain.length = str_remove(substr(molecular.formula, start=1, stop=3),pattern = "C")) %>%
  mutate(Cl.No = substr(molecular.formula,
                        start=(unlist(str_locate(molecular.formula,"Cl"))+2),
                        stop=str_length(molecular.formula))) %>%
  mutate(Cl.No = factor(Cl.No, levels = c("4","5","6","7","8","9","10","11","12")))


## Export the unprocessed data as an XLSX file  
#'Exported file will contain information about peak identities, peak areas, Cl and C numbers aggregated by sample type (i.e, Sample, Blank, QC, Calibration)  


if (any(store.lc.1M$Sample.type == "Blank")) {
  write.xlsx(spread(data = store.lc.1M[,-8]%>% filter(Sample.type == "Blank"),key = Sample.name, value = peak.area ) %>% arrange(row.identity) %>% as.data.frame(.),
             file=paste(temp.xlsx.filename,"-RAW_export.xlsx",sep=""), sheetName="Blanks", row.names=FALSE,showNA = FALSE)}
if (any(store.lc.1M$Sample.type == "Sample")) {
  write.xlsx(spread(data = store.lc.1M[,-8]%>% filter(Sample.type == "Sample"),key = Sample.name, value = peak.area )%>% arrange(row.identity) %>% as.data.frame(.),
             file=paste(temp.xlsx.filename,"-RAW_export.xlsx",sep=""), sheetName="Samples",append = TRUE, row.names=FALSE,showNA = FALSE)}
if (any(store.lc.1M$Sample.type == "Calibration")) {
  write.xlsx(spread(data = store.lc.1M[,-8]%>% filter(Sample.type == "Calibration"),key = Sample.name, value = peak.area )%>% arrange(row.identity) %>% as.data.frame(.),
             file=paste(temp.xlsx.filename,"-RAW_export.xlsx",sep=""), sheetName="Calibration",append = TRUE, row.names=FALSE,showNA = FALSE)}
if (any(store.lc.1M$Sample.type == "QC")) {
  write.xlsx(spread(data = store.lc.1M[,-8]%>% filter(Sample.type == "QC"),key = Sample.name, value = peak.area )%>% arrange(row.identity) %>% as.data.frame(.),
             file=paste(temp.xlsx.filename,"-RAW_export.xlsx",sep=""), sheetName="QC",append = TRUE, row.names=FALSE,showNA = FALSE)}


### Normalize peak areas using internal standard


store.lc.1M = store.lc.1M %>%
  mutate(Peak.area.ratio = peak.area/ISTD.area)


### Graph internal standard peak areas to check for sample inconsistencies


store.lc.1M %>%
  select(Sample.name,ISTD.area,Sample.type) %>%
  unique() %>%
  ungroup() %>%
  mutate(Mean.area = mean(ISTD.area, na.rm = TRUE)) %>%
  mutate(SD = sd(ISTD.area, na.rm = TRUE)) %>%
  ggplot(.,aes(x = Sample.name, y = ISTD.area, shape = Sample.type))+
  geom_point(size = 2)+
  ylab(label = "Internal standard area")+
  geom_hline(aes(yintercept = Mean.area), colour = "blue")+
  geom_hline(aes(yintercept = Mean.area+2*SD), colour = "red")+
  geom_hline(aes(yintercept = Mean.area-2*SD), colour = "red")+
  geom_hline(aes(yintercept = Mean.area+1*SD), colour = "darkgreen")+
  geom_hline(aes(yintercept = Mean.area-1*SD), colour = "darkgreen")+
  guides(shape = guide_legend(title = "Sample type"))+
  theme_excel()+
  scale_y_continuous(limits = c(0,max(store.lc.1M$ISTD.area*1.2)))+
  theme(panel.background = element_rect(fill = "grey90"),
        panel.grid.major.y = element_line(colour = "grey50"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank())


### Blank assignment and subtraction using a pre-defined assignment of blanks from "Reference_sheet.xlsx"  
#'Blank subtraction is not conducted for Calibration samples and samples that do not have an assigned blank file from the Reference_sheet.xlsx file 


store.lc.blank.values = store.lc.1M %>%
  filter(Sample.type == "Blank") %>%
  group_by(molecular.formula, Sample.name) %>%
  mutate(Peak.area.ratio.mean = mean(Peak.area.ratio, na.rm = TRUE)) %>%
  select(molecular.formula,Peak.area.ratio.mean,Sample.name) %>%
  mutate(Peak.area.ratio.mean = ifelse(is.nan(Peak.area.ratio.mean)==TRUE,0,Peak.area.ratio.mean)) %>%
  unique(.)
store.lc.blank.values = left_join(store.lc.blank.values, read_xlsx("Reference_sheet.xlsx", sheet = 4), by = "Sample.name") %>%
  ungroup() %>%
  select(-Sample.name)

store.lc.1M = store.lc.1M %>%
  left_join(., read_xlsx("Reference_sheet.xlsx", sheet = 4), by = "Sample.name") %>%
  left_join(.,store.lc.blank.values, by = c("molecular.formula","Blank.sample")) %>%
  mutate(Peak.area.ratio.bs = ifelse(Sample.type != "Calibration",Peak.area.ratio-Peak.area.ratio.mean,Peak.area.ratio)) %>%
  mutate(Peak.area.ratio.bs = ifelse(Peak.area.ratio.bs>0,Peak.area.ratio.bs,0)) %>%
  select(-Peak.area.ratio.mean,-Blank.sample)


### Calculate relative distribution among homologue groups


store.lc.1M = 
  store.lc.1M %>%
  filter(row.identity != "000_ISTD") %>%
  group_by(Sample.name, Chain.length) %>%
  mutate(Sum.chain.length.peak.area.ratio.bs = sum(Peak.area.ratio.bs, na.rm = TRUE)) %>%
  mutate(Sum.chain.length.peak.area.ratio = sum(Peak.area.ratio, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Sample.name, Chain.length, Cl.No) %>%
  mutate(Relative.distribution.bs = Peak.area.ratio.bs/Sum.chain.length.peak.area.ratio.bs*100) %>%
  mutate(Relative.distribution = Peak.area.ratio/Sum.chain.length.peak.area.ratio*100)

### Calculate theoretical Cl% degrees of each C~n~Cl~4~~-~~12~  homologoue groups
  
store.lc.1M = store.lc.1M %>%
  mutate(Cl.degree = (as.numeric(as.character(Cl.No))*35.5)/(as.numeric(as.character(Cl.No))*35.5+as.numeric(as.character(Chain.length))*12+(as.numeric(as.character(Chain.length))*2+2)-as.numeric(as.character(Cl.No)))*100) %>%
  mutate(Cl.degree.impact.bs = Cl.degree * Relative.distribution.bs/100) %>%
  mutate(Cl.degree.impact = Cl.degree * Relative.distribution/100) %>%
  group_by(Sample.name,Chain.length) %>%
  mutate(Chain.Cl.degree.impact.bs = sum(Cl.degree.impact.bs, na.rm = TRUE)) %>%
  mutate(Chain.Cl.degree.impact = sum(Cl.degree.impact, na.rm = TRUE)) %>%
  select(-Cl.degree.impact.bs,-Cl.degree.impact.bs)

store.lc.1M = store.lc.1M %>% unique()

#===============================================================================

## Export processed data with peak area ratios and Cl% degrees as an Excel Output  ####
  
### Output in excel by sample type (normalized by internal standard, **blank subtracted**)
  

if (any(store.lc.1M$Sample.type == "Blank")) {
  write.xlsx(spread(data = store.lc.1M[,c(1:4,6:7,9:10,12)]%>% filter(Sample.type == "Blank"),key = Sample.name, value = Peak.area.ratio.bs) %>% arrange(row.identity) %>% as.data.frame(),
             file=paste(temp.xlsx.filename,"-Subtracted_normalized_export.xlsx",sep=""), sheetName="Blanks", row.names=FALSE,showNA = FALSE)}
if (any(store.lc.1M$Sample.type == "Sample")) {
  write.xlsx(spread(data = store.lc.1M[,c(1:4,6:7,9:10,12)]%>% filter(Sample.type == "Sample"),key = Sample.name, value = Peak.area.ratio.bs)%>% arrange(row.identity) %>% as.data.frame(),
             file=paste(temp.xlsx.filename,"-Subtracted_normalized_export.xlsx",sep=""), sheetName="Samples",append = TRUE, row.names=FALSE,showNA = FALSE)}
if (any(store.lc.1M$Sample.type == "Calibration")) {
  write.xlsx(spread(data = store.lc.1M[,c(1:4,6:7,9:10,12)]%>% filter(Sample.type == "Calibration"),key = Sample.name, value = Peak.area.ratio.bs)%>% arrange(row.identity) %>% as.data.frame(),
             file=paste(temp.xlsx.filename,"-Subtracted_normalized_export.xlsx",sep=""), sheetName="Calibration",append = TRUE, row.names=FALSE,showNA = FALSE)}
if (any(store.lc.1M$Sample.type == "QC")) {
  write.xlsx(spread(data = store.lc.1M[,c(1:4,6:7,9:10,12)]%>% filter(Sample.type == "QC"),key = Sample.name, value = Peak.area.ratio.bs)%>% arrange(row.identity) %>% as.data.frame(),
             file=paste(temp.xlsx.filename,"-Subtracted_normalized_export.xlsx",sep=""), sheetName="QC",append = TRUE, row.names=FALSE,showNA = FALSE)}


### Output in excel by sample type (normalized by internal standard, **not blank subtracted**)

if (any(store.lc.1M$Sample.type == "Blank")) {
  write.xlsx(spread(data = store.lc.1M[,c(1:4,6:7,9:10,11)]%>% filter(Sample.type == "Blank"),key = Sample.name, value = Peak.area.ratio) %>% arrange(row.identity) %>% as.data.frame(),
             file=paste(temp.xlsx.filename,"-Non-Subtracted_normalized_export.xlsx",sep=""), sheetName="Blanks", row.names=FALSE,showNA = FALSE)}
if (any(store.lc.1M$Sample.type == "Sample")) {
  write.xlsx(spread(data = store.lc.1M[,c(1:4,6:7,9:10,11)]%>% filter(Sample.type == "Sample"),key = Sample.name, value = Peak.area.ratio)%>% arrange(row.identity) %>% as.data.frame(),
             file=paste(temp.xlsx.filename,"-Non-Subtracted_normalized_export.xlsx",sep=""), sheetName="Samples",append = TRUE, row.names=FALSE,showNA = FALSE)}
if (any(store.lc.1M$Sample.type == "Calibration")) {
  write.xlsx(spread(data = store.lc.1M[,c(1:4,6:7,9:10,11)]%>% filter(Sample.type == "Calibration"),key = Sample.name, value = Peak.area.ratio)%>% arrange(row.identity) %>% as.data.frame(),
             file=paste(temp.xlsx.filename,"-Non-Subtracted_normalized_export.xlsx",sep=""), sheetName="Calibration",append = TRUE, row.names=FALSE,showNA = FALSE)}
if (any(store.lc.1M$Sample.type == "QC")) {
  write.xlsx(spread(data = store.lc.1M[,c(1:4,6:7,9:10,11)]%>% filter(Sample.type == "QC"),key = Sample.name, value = Peak.area.ratio)%>% arrange(row.identity) %>% as.data.frame(),
             file=paste(temp.xlsx.filename,"-Non-Subtracted_normalized_export.xlsx",sep=""), sheetName="QC",append = TRUE, row.names=FALSE,showNA = FALSE)}

#===============================================================================

## Calibration  ####
  
### Create a separate data.frame for calibration  
  

store.lc.calibration = store.lc.1M %>%
  filter(Sample.type == "Calibration") %>%
  select(-Peak.area.ratio.bs,-Sum.chain.length.peak.area.ratio.bs,-Relative.distribution.bs,-Chain.Cl.degree.impact.bs)


### Read calibration information from "Reference_sheet.xlsx"  
#'Information: calibration mixture code, concentration

store.lc.calibration = left_join(x = store.lc.calibration, y= read_xlsx("Reference_sheet.xlsx", sheet = 2), by = "Sample.name")

### Drop unused homologoue groups (e.g., when a calibration mixture doesn't contain MCCPs)  

store.lc.calibration = store.lc.calibration %>%
  mutate(Use.homologue.group = str_detect(as.character(Chain.length.info), as.character(Chain.length))) %>%
  filter(Use.homologue.group == TRUE)

### Extract R-squared values and calculate slope  
#'intercept = 0 (forced to origin)  


temp.R.sq = store.lc.calibration %>%
  ungroup() %>%
  select(Calibration.type,Chain.length,Calibration.level,Sum.chain.length.peak.area.ratio) %>%
  unique() %>%
  group_by(Calibration.type,Chain.length) %>%
  do(InfoList = glance(lm(.,formula = Sum.chain.length.peak.area.ratio~Calibration.level-1))) %>%
  unnest(InfoList) %>%
  select(Calibration.type,Chain.length,r.squared)

temp.Slope = store.lc.calibration %>%
  ungroup() %>%
  select(Calibration.type,Chain.length,Calibration.level,Sum.chain.length.peak.area.ratio) %>%
  unique() %>%
  group_by(Calibration.type,Chain.length) %>%
  do(InfoList = tidy(lm(.,formula = Sum.chain.length.peak.area.ratio~Calibration.level-1))) %>%
  unnest(InfoList) %>%
  select(Calibration.type,Chain.length,estimate) %>%
  mutate(Slope = estimate) %>%
  select(-estimate)

store.lc.calibration = left_join(store.lc.calibration,
                                 left_join(temp.Slope,temp.R.sq,by = c("Calibration.type","Chain.length")),by = c("Calibration.type","Chain.length"))
rm(temp.R.sq,temp.Slope)


### Create calibration plots and export as .pdf  
#If there are vSCCPs, LCCPs or vLCCPs implemented in the workflow, the following for loop has to be adjusted  
#For example, if only SCCPs are used: `for (i in 10:13)` or if vSCCPs/SCCPs/MCCPs are used: `for (i in 9:17)`  

Conc.vs.Resp.plots = list()

for (i in 10:17) {
  Conc.vs.Resp.plots[[i]] =  store.lc.calibration %>%
    filter(Chain.length == i) %>%
    ungroup() %>%
    select(Calibration.level,Sum.chain.length.peak.area.ratio,Calibration.type) %>%
    unique() %>%
    ggplot(., aes(x = Calibration.level, y = Sum.chain.length.peak.area.ratio))+
    geom_point()+
    theme_excel()+
    theme(panel.background = element_rect(fill = "grey90"),
          panel.grid.major.y = element_line(colour = "grey50"))+
    geom_smooth(method = "lm",formula =y~x-1)+
    ylab(label = "Response\n(peak area sum / ISTD peak area)")+
    xlab(label = "Conc., ng/??L")+
    stat_poly_eq(formula = y~x-1,
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                 parse = TRUE, coef.digits = 3, rr.digits = 3)+
    facet_wrap(~Calibration.type, scales = "free_y")

Calibration.main = store.lc.calibration %>%
  ungroup() %>%
  select(Chain.length,Chain.Cl.degree.impact,Calibration.type,Slope,r.squared) %>%
  unique() %>%
  group_by(Chain.length,Calibration.type) %>%
  mutate(Chain.Cl.degree.impact = mean(Chain.Cl.degree.impact, na.rm = TRUE)) %>%
  unique() %>%
  group_by(Chain.length) %>% 
  mutate(Min.Cl.degree = min(Chain.Cl.degree.impact)) %>%
  mutate(Max.Cl.degree = max(Chain.Cl.degree.impact))
}


Chlorine.vs.Resp.plots = list()
for (i in 10:17) {
  Chlorine.vs.Resp.plots[[i]] = Calibration.main %>%
    filter(Chain.length == i) %>%
    ungroup() %>%
    ggplot(., aes(x = Chain.Cl.degree.impact, y = Slope))+
    geom_point()+
    theme_light()+
    geom_smooth(method = "lm",formula =y~x, colour = "red")+
    ggtitle(label = paste("Homologue group with",i,"carbons"))+
    ylab(label = "Response factor")+
    xlab(label = "Calculated chlorine content, %")+
    stat_poly_eq(formula = y~x,
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                 parse = TRUE, coef.digits = 3, rr.digits = 3)
}

pdf("Calibration_report.pdf", width = 6, height = 8)
for (i in 10:17) {
  grid.arrange(Chlorine.vs.Resp.plots[[i]], Conc.vs.Resp.plots[[i]], ncol=1)
}
dev.off()


### Fetch R-squared and slope values from calibration and join with main dataset


temp.R.sq = Calibration.main %>%
  group_by(Chain.length) %>%
  do(InfoList = glance(lm(.,formula = Slope~Chain.Cl.degree.impact))) %>%
  unnest(InfoList) %>%
  select(Chain.length,r.squared)

temp.Slope = Calibration.main %>%
  group_by(Chain.length) %>%
  do(InfoList = tidy(lm(.,formula = Slope~Chain.Cl.degree.impact))) %>%
  unnest(InfoList)%>%
  select(Chain.length,term,estimate) %>%
  mutate(term = ifelse(term == "Chain.Cl.degree.impact","Slope","Intercept"))%>%
  spread(term,estimate)

Calibration.main = left_join(Calibration.main %>% select(-Slope,-r.squared),
                             left_join(temp.Slope,temp.R.sq,by = c("Chain.length")),by = c("Chain.length"))
rm(temp.R.sq,temp.Slope)

Calibration.main = Calibration.main %>%
  ungroup() %>%
  select(-Chain.Cl.degree.impact,-Calibration.type) %>%
  unique()









#===============================================================================

## Quantification of samples  ####
  
### Create data frame for results and join with sample and calibration data 
  

Results.sheet = store.lc.1M %>% filter(Sample.type != "Calibration") %>%
  select(Sample.name,Sample.type,Chain.length,Sum.chain.length.peak.area.ratio,Sum.chain.length.peak.area.ratio.bs,
         Chain.Cl.degree.impact, Chain.Cl.degree.impact.bs) %>%
  unique()

Results.sheet = left_join(Results.sheet,Calibration.main, by = "Chain.length") %>%
  mutate(CP.type = ifelse(Chain.length %in% c(10,11,12,13),"SCCPs",ifelse(Chain.length %in% c(14,15,16,17),"MCCPs","LCCPs")))


### Quantify samples  

#'When calculated Cl% degree in sample is lower than min(Cl%) degree in the calibration: Cl%, sample ==> min(Cl%), calibration  
#'When calculated Cl% degree in sample is higher than max(Cl%) degree in the calibration: Cl%, sample ==> max(Cl%), calibration  
#'This is done to avoid unrealistic response factors that fall outside the Cl% range of calibration mixtures  


Results.sheet = Results.sheet %>%
  ungroup()%>%
  mutate(Chain.Cl.degree.f = ifelse(Chain.Cl.degree.impact<Min.Cl.degree,Min.Cl.degree,Chain.Cl.degree.impact))%>%
  mutate(Chain.Cl.degree.f.bs = ifelse(Chain.Cl.degree.impact.bs<Min.Cl.degree,Min.Cl.degree,Chain.Cl.degree.impact.bs)) %>%
  mutate(Chain.Cl.degree.f = ifelse(Chain.Cl.degree.f>Max.Cl.degree,Max.Cl.degree,Chain.Cl.degree.f))%>%
  mutate(Chain.Cl.degree.f.bs = ifelse(Chain.Cl.degree.f.bs>Max.Cl.degree,Max.Cl.degree,Chain.Cl.degree.f.bs)) %>%
  mutate(Chain.RF = (Intercept+Slope*Chain.Cl.degree.f)) %>%
  mutate(Chain.RF.bs = (Intercept+Slope*Chain.Cl.degree.f.bs)) %>%
  mutate(Chain.RF = ifelse(Chain.RF<1,1,Chain.RF)) %>%
  mutate(Chain.RF.bs = ifelse(Chain.RF.bs<1,1,Chain.RF.bs)) %>%
  mutate(Conc.nguL = Sum.chain.length.peak.area.ratio/Chain.RF) %>%
  mutate(Conc.nguL.bs = Sum.chain.length.peak.area.ratio.bs/Chain.RF.bs) %>%
  group_by(Sample.name,CP.type) %>%
  mutate(Sum.CPs.nguL = sum(Conc.nguL, na.rm = TRUE)) %>%
  mutate(Sum.CPs.nguL.bs = sum(Conc.nguL.bs, na.rm = TRUE))









#===============================================================================

## Export quantification results  ####
  
### Export blank subtracted results 
#' The results are shown as ng per extract.
#' If you wish to obtain final concentration in the extract, this value has to be divided by sample volume (??L)  


rbind(Results.sheet %>%
        ungroup() %>%
        select(Sample.name, Chain.length, Conc.nguL.bs),
      Results.sheet %>%
        ungroup() %>%
        select(Sample.name, CP.type, Sum.CPs.nguL.bs) %>%
        mutate(Chain.length = CP.type) %>%
        mutate(Conc.nguL.bs = Sum.CPs.nguL.bs) %>%
        select(-CP.type,-Sum.CPs.nguL.bs) %>%
        unique()) %>%
  ungroup() %>%
  mutate(Conc.nguL.bs = Conc.nguL.bs*Final.extract.volume) %>%
  spread(key = Sample.name, value = Conc.nguL.bs) %>%
  arrange(Chain.length) %>% as.data.frame() %>%
  write.xlsx(.,file=paste(temp.xlsx.filename,"-Concentration_results.xlsx",sep=""),
             sheetName="Blank subtract", row.names=FALSE,showNA = FALSE)


### Export results without blank subtraction  
#' The results are shown as ng per extract. If you wish to obtain final concentration in the extract, this value has to be divided by sample volume (??L)  


rbind(Results.sheet %>%
        ungroup() %>%
        select(Sample.name, Chain.length, Conc.nguL),
      Results.sheet %>%
        ungroup() %>%
        select(Sample.name, CP.type, Sum.CPs.nguL) %>%
        mutate(Chain.length = CP.type) %>%
        mutate(Conc.nguL = Sum.CPs.nguL) %>%
        select(-CP.type,-Sum.CPs.nguL) %>%
        unique()) %>%
  ungroup() %>%
  mutate(Conc.nguL = Conc.nguL*Final.extract.volume) %>%
  spread(key = Sample.name, value = Conc.nguL) %>%
  arrange(Chain.length) %>% as.data.frame() %>%
  write.xlsx(.,file=paste(temp.xlsx.filename,"-Concentration_results.xlsx",sep=""),
             sheetName="No blank subtract", row.names=FALSE,showNA = FALSE)








#===============================================================================

## Visualization of results  ####
  
### Visualization of concentrations  
  
#'Generate a list of sample names  


Sample.list = Results.sheet$Sample.name %>% unique(.) %>% as.character(.)
Sample.list

### Select samples for visualization  
#'Examples:  
#'Sample.subset = Sample.list select all samples  
#'Sample.subset = Sample.list[1:10] select samples from 1 to 10  
#'Sample.subset = Sample.list[c(1:10,14,17)] select samples from 1 to 10 as well as sample 14 and 17  

Sample.subset = Sample.list[1:6]

### Concentration of SCCPs, MCCPs and LCCPs (blank subtracted)  

Results.sheet %>%
  ungroup() %>%
  filter(Sample.name %in% Sample.subset) %>%
  select(Sample.name,CP.type,Sum.CPs.nguL.bs) %>%
  unique(.) %>%
  ggplot(., aes(x=Sample.name, y=Sum.CPs.nguL.bs, fill = CP.type))+
  geom_bar(stat = "identity",position = "dodge", colour = "grey10")+
  scale_fill_brewer(type = "qual",palette =3)+
  ylab("ng/??L")+
  guides(fill = guide_legend(title = "Type"))+
  theme_clean()+
  theme(axis.title.x = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90,hjust =1, vjust=0.5))


### Concentration of SCCPs, MCCPs and LCCPs (without blank subtraction)  


Results.sheet %>%
  ungroup() %>%
  filter(Sample.name %in% Sample.subset) %>%
  select(Sample.name,CP.type,Sum.CPs.nguL) %>%
  unique(.) %>%
  ggplot(., aes(x=Sample.name, y=Sum.CPs.nguL, fill = CP.type))+
  geom_bar(stat = "identity",position = "dodge", colour = "grey10")+
  scale_fill_brewer(type = "qual",palette =3)+
  ylab("ng/??L")+
  guides(fill = guide_legend(title = "Type"))+
  theme_clean()+
  theme(axis.title.x = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90,hjust =1, vjust=0.5))


### Concentration of CP homologue groups (by chain length, blank subtracted)  


Results.sheet %>%
  ungroup() %>%
  filter(Sample.name %in% Sample.subset) %>%
  select(Sample.name,Chain.length,Conc.nguL.bs) %>%
  unique(.) %>%
  ggplot(., aes(x=Sample.name, y=Conc.nguL.bs, fill = Chain.length))+
  geom_bar(stat = "identity", colour = "grey10")+
  scale_fill_brewer(type = "qual",palette =3)+
  ylab("ng/??L")+
  guides(fill = guide_legend(title = "Type"))+
  theme_clean()+
  theme(axis.title.x = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90,hjust =1, vjust=0.5))


### Concentration of CP homologue groups (by chain length, without blank subtraction)  


Results.sheet %>%
  ungroup() %>%
  filter(Sample.name %in% Sample.subset) %>%
  select(Sample.name,Chain.length,Conc.nguL) %>%
  unique(.) %>%
  ggplot(., aes(x=Sample.name, y=Conc.nguL, fill = Chain.length))+
  geom_bar(stat = "identity", colour = "grey10")+
  scale_fill_brewer(type = "qual",palette =3)+
  ylab("ng/??L")+
  guides(fill = guide_legend(title = "Type"))+
  theme_clean()+
  theme(axis.title.x = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90,hjust =1, vjust=0.5))







#===============================================================================

## Visualize homologue profiles  ####
  
#'Enter sample name/-s (up to 5 samples):  
  

Sample.name.select = c("Lard B","Lard D","Pork 1")


### Visualize peak realative distribution within each Cn CP homologue groups (%, blank subtracted)  

store.lc.1M %>%
  ungroup() %>%
  filter(Sample.name %in% Sample.name.select) %>%
  mutate(Cl.No = as.numeric(Cl.No)) %>%
  filter(Cl.No %in% c(4,5,6,7,8,9,10,11,12)) %>%
  select(Sample.name,Chain.length,Cl.No,Peak.area.ratio.bs) %>%
  mutate(Peak.area.ratio.bs = ifelse(is.na(Peak.area.ratio.bs)==TRUE,0,Peak.area.ratio.bs)) %>%
  mutate(Profile.total = Peak.area.ratio.bs/sum(Peak.area.ratio.bs, na.rm = TRUE)*100) %>%
  group_by(Chain.length) %>%
  mutate(Profile.chain = Peak.area.ratio.bs/sum(Peak.area.ratio.bs, na.rm = TRUE)*100) %>%
  ggplot(., aes(x=Cl.No, y=Profile.chain, fill = Sample.name))+
  geom_bar(stat = "identity", position = "dodge", colour = "grey10")+
  guides(fill = guide_legend(title = "Sample name"))+
  scale_fill_brewer(type = "qual",palette =3)+
  facet_wrap(~paste("C",Chain.length, sep =""), scales = "free_x")+
  theme_clean()+
  theme(axis.title = element_blank(), legend.position = "bottom")

### Visualize peak realative distribution within each Cn CP homologue groups (%, without blank subtraction)  

store.lc.1M %>%
  ungroup() %>%
  filter(Sample.name %in% Sample.name.select) %>%
  mutate(Cl.No = as.numeric(Cl.No)) %>%
  filter(Cl.No %in% c(4,5,6,7,8,9,10,11,12)) %>%
  select(Sample.name,Chain.length,Cl.No,Peak.area.ratio) %>%
  mutate(Peak.area.ratio = ifelse(is.na(Peak.area.ratio)==TRUE,0,Peak.area.ratio)) %>%
  mutate(Profile.total = Peak.area.ratio/sum(Peak.area.ratio, na.rm = TRUE)*100) %>%
  group_by(Chain.length) %>%
  mutate(Profile.chain = Peak.area.ratio/sum(Peak.area.ratio, na.rm = TRUE)*100) %>%
  ggplot(., aes(x=Cl.No, y=Profile.chain, fill = Sample.name))+
  geom_bar(stat = "identity", position = "dodge", colour = "grey10")+
  guides(fill = guide_legend(title = "Sample name"))+
  scale_fill_brewer(type = "qual",palette =3)+
  facet_wrap(~paste("C",Chain.length, sep =""), scales = "free_x")+
  theme_clean()+
  theme(axis.title = element_blank(), legend.position = "bottom")


### Visualize peak area distribution (absolute peak areas, blank subtracted)  


store.lc.1M %>%
  ungroup() %>%
  filter(Sample.name %in% Sample.name.select) %>%
  mutate(Cl.No = as.numeric(Cl.No)) %>%
  filter(Cl.No %in% c(4,5,6,7,8,9,10,11,12)) %>%
  select(Sample.name,Chain.length,Cl.No,Peak.area.ratio.bs) %>%
  mutate(Peak.area.ratio.bs = ifelse(is.na(Peak.area.ratio.bs)==TRUE,0,Peak.area.ratio.bs)) %>%
  mutate(Profile.total = Peak.area.ratio.bs) %>%
  group_by(Chain.length) %>%
  mutate(Profile.chain = Peak.area.ratio.bs) %>%
  ggplot(., aes(x=Cl.No, y=Profile.chain, fill = Sample.name))+
  geom_bar(stat = "identity", position = "dodge", colour = "grey10")+
  guides(fill = guide_legend(title = "Sample name"))+
  scale_fill_brewer(type = "qual",palette =3)+
  facet_wrap(~paste("C",Chain.length, sep =""), scales = "free_x")+
  theme_clean()+
  theme(axis.title = element_blank(), legend.position = "bottom")


### Visualize peak area distribution (absolute peak areas, without blank subtraction)  

store.lc.1M %>%
  ungroup() %>%
  filter(Sample.name %in% Sample.name.select) %>%
  mutate(Cl.No = as.numeric(Cl.No)) %>%
  filter(Cl.No %in% c(4,5,6,7,8,9,10,11,12)) %>%
  select(Sample.name,Chain.length,Cl.No,Peak.area.ratio) %>%
  mutate(Peak.area.ratio = ifelse(is.na(Peak.area.ratio)==TRUE,0,Peak.area.ratio)) %>%
  mutate(Profile.total = Peak.area.ratio) %>%
  group_by(Chain.length) %>%
  mutate(Profile.chain = Peak.area.ratio) %>%
  ggplot(., aes(x=Cl.No, y=Profile.chain, fill = Sample.name))+
  geom_bar(stat = "identity", position = "dodge", colour = "grey10")+
  guides(fill = guide_legend(title = "Sample name"))+
  scale_fill_brewer(type = "qual",palette =3)+
  facet_wrap(~paste("C",Chain.length, sep =""), scales = "free_x")+
  theme_clean()+
  theme(axis.title = element_blank(), legend.position = "bottom")


