library(MSDtracking)

# BRCA2 HU ----------------------------------------------------------------
all_data <- list()
datasets <- c("D:/Stack/Genetics/180110 BRCA2 WTG10_dDdC_F4",
              "O:/Maarten/Genetics/170711 BRCA2 dDBD CTD Halo HU",
              "O:/Maarten/Genetics/170523 BRCA2 tracking HU",
              "D:/Stack/Genetics/180502 BRCA2 Halo dDBD+CTD HU/",
              "D:/Stack/Genetics/180606 BRCA2 WT dCTD Halo/")

all_data[["WT -HU"]] <- merge_dataset(datasets[c(2,3,4,5)],c("WT G10 noHU","WT E10 noHU","WT E10 noHU","WT G10 -HU","WT G10 -HU")[c(2,3,4,5)],con_name="WT -HU")

all_data[["WT +HU"]] <- merge_dataset(datasets[c(2,3,4,5)],c("WT G10 +HU 2h","WT E10 HU","WT E10 HU","WT G10 +HU","WT G10 +HU")[c(2,3,4,5)],con_name="WT +HU")

all_data[["dDBDdCTD -HU"]] <- merge_dataset(datasets[c(2,3,4)],c("dDBDdCTD F4 noHU","dDBDdCTD G9 noHU","dDBD G9 noHU","dDBDdCTD F4 -HU")[c(2,3,4)],con_name="dDBDdCTD -HU")

all_data[["dDBDdCTD +HU"]] <- merge_dataset(datasets[c(2,3,4)],c("dDBDdCTD F4 +HU 2h","dDBDdCTD G9 HU","dDBD G9 HU","dDBDdCTD F4 +HU")[c(2,3,4)],con_name="dDBDdCTD +HU")

all_data[["dDBD -HU"]] <- merge_dataset(datasets[2],c("dDBD+CTD A4 noHU"),con_name="dDBD -HU")

all_data[["dDBD +HU"]] <- merge_dataset(datasets[2],c("dDBD+CTD A4 HU"),con_name="dDBD +HU")

all_data[["dCTD -HU"]] <- merge_dataset(datasets[5],c("dCTD A2 -HU"),con_name="dCTD -HU")

all_data[["dCTD +HU"]] <- merge_dataset(datasets[5],c("dCTD A2 +HU"),con_name="dCTD +HU")

save(all_data,file = file.path("D:/Stack/Genetics/Tracking data merged/","BRCA2_HU_all_data.Rdata"))
msd_histogram(all_data,directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="BRCA2_HU",threshold=0.05)
msd_histogram(all_data[c(1,2)],directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="BRCA2_HU_WT",threshold=0.05)
msd_histogram(all_data[c(3,4)],directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="BRCA2_HU_dDBDdCTD",threshold=0.05)
msd_histogram(all_data[c(5,6)],directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="BRCA2_HU_dDBD",threshold=0.05)
msd_histogram(all_data[c(1,2,3,4)],directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="BRCA2_HU+dDBDdCTD_WT",threshold=0.05)

# BRCA2 MMC ---------------------------------------------------------------
all_data <- list()
datasets <- c("F://170517 BRCA2 tracking MMC",
              "F://170623 BRCA2 MMC tracking",
              "D:/Stack/Genetics/170809 tracking MMC WT_dDBD_+CTD")
all_data[["WT -MMC"]] <- merge_dataset(datasets[c(1,2,3)],c("WT E10 5nM noMMC","WT E10 noMMC2h","WT E10 noMMC")[c(1,2,3)],con_name="WT -MMC")

all_data[["WT +MMC"]] <- merge_dataset(datasets[c(1,2,3)],c("WT E10 5nM 1ug MMC","WT E10 MMC","WT E10 MMC")[c(1,2,3)],con_name="WT +MMC")

all_data[["dDBDdCTD -MMC"]] <- merge_dataset(datasets[c(1,2,3)],c("dDBDdCTD G9 5nM noMMC","dDBDcCTD G9 noMMC","dDBDdCTD G9 noMMC")[c(1,2,3)],con_name="dDBDdCTD -MMC")

all_data[["dDBDdCTD +MMC"]] <- merge_dataset(datasets[c(1,2,3)],c("dDBDdCTD G9 5nM 1ug MMC","dDBDdCTD G9 MMC","dDBDdCTD G9 MMC")[c(1,2,3)],con_name="dDBDdCTD +MMC")

all_data[["dDBD -MMC"]] <- merge_dataset("D:/Stack/Genetics/170809 tracking MMC WT_dDBD_+CTD",c("dDBD A4 noMMC"),con_name="dDBD -MMC")

all_data[["dDBD +MMC"]] <- merge_dataset("D:/Stack/Genetics/170809 tracking MMC WT_dDBD_+CTD",c("dDBD A4 MMC"),con_name="dDBD +MMC")

save(all_data,file = file.path("D:/Stack/Genetics/Tracking data merged/","BRCA2_MMC_all_data.Rdata"))
msd_histogram(all_data,directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="MMC",threshold=0.05)
msd_histogram(all_data[c(1,2)],directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="BRCA2_MMC_WT",threshold=0.05)
msd_histogram(all_data[c(3,4)],directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="BRCA2_MMC_dDBDdCTD",threshold=0.05)
msd_histogram(all_data[c(1,2,3,4)],directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="BRCA2_MMC_WT+dDBDdCTD",threshold=0.05)
msd_histogram(all_data[c(5,6)],directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="BRCA2_MMC_dDBD",threshold=0.05)

# BRCA2 IR ----------------------------------------------------------------

all_data <- list()
datasets <- c("D:/Stack/Genetics/170510 particle tracking BRCA2/","F://180117 BRCA2 dDdC IR tracking")
all_data[["WT -IR"]] <- merge_dataset(datasets,c("BRCA2 WT Halo E10 noIR","WT G10 -IR"),con_name="WT -IR")

all_data[["WT +IR"]] <- merge_dataset(datasets,c("BRCA2 WT Halo E10 5gyIR 2h","WT G10 +IR"),con_name="WT +IR")

all_data[["dDBDdCTD -IR"]] <- merge_dataset(datasets,c("BRCA2 dDBDdCTD Halo G9 noIR","dDBDdCTD F4 -IR"),con_name="dDBDdCTD -IR")

all_data[["dDBDdCTD +IR"]] <- merge_dataset(datasets,c("BRCA2 dDBDdCTD Halo G9 5gyIR 2h","dDBDdCTD F4 +IR"),con_name="dDBDdCTD +IR")

save(all_data,file = file.path("D:/Stack/Genetics/Tracking data merged/","BRCA2_IR_all_data.Rdata"))
load(file.path("D:/Stack/Genetics/Tracking data merged/","BRCA2_IR_all_data.Rdata"))

msd_histogram(all_data,directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="BRCA2_IR",threshold=0.05)
msd_histogram(all_data[c(1,2)],directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="BRCA2_WT_IR",threshold=0.05)
msd_histogram(all_data[c(3,4)],directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="BRCA2_dDBDdCTD_IR",threshold=0.05)


# ##PALB2 HU --------------------------------------------------------------
all_data <- list()
stop("TO DO")
datasets <- c("F://170517 BRCA2 tracking MMC","F://170623 BRCA2 MMC tracking","D:/Stack/Genetics/170809 tracking MMC WT_dDBD_+CTD")
all_data[["WTnoMMC"]] <- merge_dataset(datasets,c("WT E10 5nM noMMC","WT E10 noMMC2h","WT E10 noMMC"),con_name="WT -MMC")

all_data[["WTMMC"]] <- merge_dataset(datasets,c("WT E10 5nM 1ug MMC","WT E10 MMC","WT E10 MMC"),con_name="WT +MMC")


# PALB2 MMC ---------------------------------------------------------------
all_data <- list()
datasets <- c("D:/Stack/Genetics/171024 PALB2-Halo tracking")
all_data[["PALB2 -MMC"]] <- merge_dataset(datasets,c("PALB2 F6 mock"),con_name="PALB2 -MMC")

all_data[["PALB2 +MMC"]] <- merge_dataset(datasets,c("PALB2 F6 MMC 2h"),con_name="PALB2 +MMC")

msd_histogram(all_data,directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="PALB2 MMC",threshold=0.05)


# Halo-BRCA2 WT d1-40 -----------------------------------------------------
all_data <- list()
datasets <- c("D:/Stack/Genetics/171130 d1-40 Halo-BRCA2","D:/Stack/Genetics/171130 N-BRCA2 Halo H8 HU")
all_data[["Halo-BRCA2"]] <- merge_dataset(datasets[2],c("-HU"),con_name="Halo-BRCA2")

all_data[["Halo-BRCA2 d1-40"]] <- merge_dataset(datasets[1],c("F9"),con_name="Halo-BRCA2 d1-40")

msd_histogram(all_data,directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="HaloTag BRCA2 d1-40",threshold=0.05)




##Old code
# directory <- "D:/Stack/Genetics/171208 N-BRCA2 Halo H8 MMC/"
# load(file.path(directory,"msd_fit_all.Rdata"))
# msd_fit_all$`-MMC`$cellID <- paste0('171208_',msd_fit_all$`-MMC`$cellID)
# msd_fit_all$`+MMC`$cellID <- paste0('171208_',msd_fit_all$`+MMC`$cellID)
#
# all$`-MMC` <- msd_fit_all$`-MMC`
# all$`+MMC` <- msd_fit_all$`+MMC`
#
# directory <- "D:/Stack/Genetics/171201 N-BRCA2 Halo H8 MMC/"
# load(file.path(directory,"msd_fit_all.Rdata"))
# msd_fit_all$`-MMC`$cellID <- paste0('171201_',msd_fit_all$`-MMC`$cellID)
# msd_fit_all$`+MMC`$cellID <- paste0('171201_',msd_fit_all$`+MMC`$cellID)
#
# all$`-MMC` <- rbind(all$`-MMC`,msd_fit_all$`-MMC`)
# all$`+MMC` <- rbind(all$`+MMC`,msd_fit_all$`+MMC`)
#
# directory <- "D:/Stack/Genetics/171117 N-BRCA2 MMC/"
# load(file.path(directory,"msd_fit_all.Rdata"))
# msd_fit_all$`H8 -MCC`$cellID <- paste0('171117_',msd_fit_all$`H8 -MCC`$cellID)
# msd_fit_all$`H8 +MCC`$cellID <- paste0('171117_',msd_fit_all$`H8 +MCC`$cellID)
#
# all$`-MMC` <- rbind(all$`-MMC`,msd_fit_all$`H8 -MCC`)
# all$`+MMC` <- rbind(all$`+MMC`,msd_fit_all$`H8 +MCC`)
#
# all$`-MMC`$condition <- "-MMC"
# all$`+MMC`$condition <- "+MMC"
# msd_fit_all <-all



# H2B/NLS ---------------------------------------------------------------------
all_data <- list()
#all_data[["H2B"]] <- merge_dataset(datasets[1],c("H2B"),con_name="H2B")
load("Y:/Maarten/Genetics/170809 tracking MMC WT_dDBD_+CTD/H2B/msd_fit.Rdata")
all_data[["H2B"]] <- ldply(tracks)
all_data[["H2B"]]$condition <- "H2B"

load("Y:/Maarten/Genetics/170512 BRCA2 slow tracking/Halo-NLS/30 ms/msd_fit.Rdata")
all_data[["NLS"]] <- ldply(tracks)
all_data[["NLS"]]$condition <- "NLS"
msd_histogram(all_data,directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="H2B-NLS",threshold=0.05)


#BRC3 preliminary
all_data <- list()
datasets <- c("D:/Stack/Genetics/180502 BRCA2 Halo dDBD+CTD HU","D:/Stack/Genetics/180502 BRCA2 Halo WT BRC3-F2A-GFP",
              "D:/Stack/Genetics/180426 BRCA2 Halo dDBD+CTD HU","D:/Stack/Genetics/180426 BRC3-F2A-GFP/WT G10 stable #2")
all_data[["WT -BRC3"]] <- merge_dataset(datasets[c(1,3)],c("WT G10 -HU","BRCA2 WT G10 -HU"),con_name="WT -BRC3")

all_data[["WT +BRC3"]] <- merge_dataset(datasets[c(2,4)],c("+BRC3","+BRC3"),con_name="WT +BRC3")
msd_histogram(all_data,directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="BRC3",threshold=0.05)

