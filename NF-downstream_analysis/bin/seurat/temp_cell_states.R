cluster_res = list('hh4' = 0.3, 'hh5' = 0.6)

cluster_order = list('hh4' = c(2,0,1,3))

cell_classification_quantile = list( 'hh4' = 0.6, 'hh5' = 0.6)

hh5_cell_type_markers = list(early_caudal_neural = c('GBX2', 'SP5', 'HOXB1', 'CDX2', 'SOX2', 'SOX21', 'SOX3'),
                             early_neural_plate = c('OTX2', 'SOX2', 'SOX21', 'SOX3'),
                             early_pNPB = c("PAX7", "MSX1", "GBX2", 'SP5', "DLX5", "DLX6", "TFAP2A", "TFAP2C", "PRDM1", "SOX2", 'SOX3', "SOX21"),
                             early_aNPB = c("SIX3", "OTX2", "DLX5", "DLX6", "TFAP2A", "TFAP2C", "PRDM1", "SOX2", 'SOX3', "SOX21"),
                             early_NPB = c("DLX5", "DLX6", "TFAP2A", "TFAP2C", "PRDM1", "SOX2", 'SOX3', "SOX21"),
                             early_NNE = c('ASTL', "DLX5", "DLX6", 'TFAP2A', "TFAP2C", "GATA2", "GATA3", "EPAS1"))

hh4_cell_type_markers = list(  node = c('EOMES', 'ADMP', 'CHRD', 'TBX6'),
                               early_neural = c("SOX2", "SOX3", 'OTX2', 'EPCAM', 'MAFA', 'FRZB', "YEATS4", 'SOX13', 'SOX11'), # many from trevers 2021 / katherine thesis
                               early_border = c("SOX2", "SOX3", 'OTX2', 'EPCAM', 'MAFA', 'FRZB', "YEATS4", 'SOX13', 'SOX11', "DLX5", "DLX6", "GATA2", "GATA3"),
                               early_non_neural = c("DLX5", "DLX6", "GATA2", "GATA3"),
                               extra_embryonic = c('VGLL1', 'GRHL3', 'GATA2', 'GATA3'))

hh6_cell_type_markers = list(  non_neural = c("EPAS1", "MSX2", "GATA2", "GATA3", "GRHL3", "EPAS1", "MSX2"),#krt19 keith mclaren 2003
                               NPB = c("SIX1", "EYA2", "DLX5", "DLX6", "TFAP2B", "TFAP2A", "TFAP2C", "PRDM1", "SOX2", "SOX3", 'SOX21', "MSX2"),
                               pNPB = c("PAX7", "MSX1", "GBX2", "SIX1", "EYA2", "DLX5", "DLX6", "TFAP2B", "TFAP2A", "TFAP2C", "PRDM1", "SOX2", "SOX3",'SOX21', "MSX2"), # check ZIC1 expression
                               aNPB = c("SIX3", "PAX6", "OTX2", "SIX1", "EYA2", "DLX5", "DLX6", "TFAP2B", "TFAP2A", "TFAP2C", "PRDM1", "SOX2", "SOX3",'SOX21', "MSX2"),
                               early_aPPR = c("SIX1", "EYA2", "DLX3", "DLX5", "DLX6", "PRDM1", "SIX3", "PAX6", "HESX1", "OTX2", "TFAP2A"), # PNOC, SSTR5
                               early_pPPR = c("SIX1", "EYA2", "DLX3", "DLX5", "DLX6", "PRDM1", "GBX2", "TFAP2A"), # FOXI3
                               early_PPR = c("SIX1", "EYA2", "DLX3", "DLX5", "DLX6", "PRDM1", "TFAP2A"), # TFAPs
                               neural_progenitors = c("SOX2", "SOX21", "LMO1", "ZEB2", "SOX1", "SOX3", "FRZB"),
                               a_neural_progenitors = c("OTX2", "SIX3", "HESX1", "SOX2", "SOX21", "LMO1", "ZEB2", "SOX3", "FRZB"),
                               p_neural_progenitors = c("GBX2", "SOX2", "SOX21", "LMO1", "ZEB2", "SOX3", "FRZB"),
                               early_hindbrain = c("GBX2", "HOXA2", "HOXA3", "HOXB2", "KROX20", "SOX2", "SOX21", "LMO1", "ZEB2", "SOX1", "SOX3", "FRZB"),
                               early_midbrain = c("WNT4", "PAX2", "FGF8", "WNT1", "OTX2", "SOX2", "SOX21", "LMO1", "ZEB2", "SOX1", "SOX3", "FRZB"),
                               early_forebrain = c("PAX6" , "SIX3", "OTX2", "SOX2", "SOX21", "LMO1", "ZEB2", "SOX1", "SOX3", "FRZB")) #TLX1)