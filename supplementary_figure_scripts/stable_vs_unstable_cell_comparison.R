#analyze QC differences between stable and unstably assigned cells 


ileum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_ileum_cell_annotation_udpate_11_26_24.rds")

# Extract metadata
meta_stable_filt <- ileum@meta.data

# Compute mean mitochondrial percentage per sample
mean_mt_per_sample_stable_filt <- meta_stable_filt %>%
  group_by(orig.ident) %>%  # replace with your sample column name
  summarise(mean_percent_mt_stable_filt = mean(percent.mt, na.rm = TRUE))

#merge mean mt per sample 
merge_mt_percentage_filt <- merge(mean_mt_per_sample_unstable, mean_mt_per_sample_stable_filt, by = "orig.ident")


#t test
t.test(merge_mt_percentage_filt$mean_percent_mt_unstable, merge_mt_percentage_filt$mean_percent_mt_stable_filt, paired = TRUE)

#t = 3.6884, df = 26, p-value = 0.001048; mean difference: 1.262635

#nFeature_RNA

# Compute mean nFeature per sample
mean_feature_per_sample_stable_filt <- meta_stable_filt %>%
  group_by(orig.ident) %>%  # replace with your sample column name
  summarise(mean_percent_feature_stable_filt = mean(nFeature_RNA, na.rm = TRUE))


# Compute mean nFeature per sample
mean_feature_per_sample_unstable <- meta_unstable %>%
  group_by(orig.ident) %>%  # replace with your sample column name
  summarise(mean_percent_feature_unstable = mean(nFeature_RNA, na.rm = TRUE))


#merge mean nfeature per sample 
merge_feature_filt <- merge(mean_feature_per_sample_unstable, mean_feature_per_sample_stable_filt, by = "orig.ident")


#t test
t.test(merge_feature_filt$mean_percent_feature_unstable, merge_feature_filt$mean_percent_feature_stable_filt, paired = TRUE)

#t = -3.4353, df = 26, p-value = 0.001998; mean difference -188.5721 

#nCount_RNA

# Compute mean nCount per sample
mean_count_per_sample_stable_filt <- meta_stable_filt %>%
  group_by(orig.ident) %>%  # replace with your sample column name
  summarise(mean_count_stable_filt = mean(nCount_RNA, na.rm = TRUE))


# Compute mean nCount per sample
mean_count_per_sample_unstable <- meta_unstable %>%
  group_by(orig.ident) %>%  # replace with your sample column name
  summarise(mean_count_unstable = mean(nCount_RNA, na.rm = TRUE))


#merge mean nCount per sample 
merge_count_filt <- merge(mean_count_per_sample_unstable, mean_count_per_sample_stable_filt, by = "orig.ident")

#t test
t.test(merge_count_filt$mean_count_unstable, merge_count_filt$mean_count_stable_filt, paired = TRUE)

#t = -5.0285, df = 26, p-value = 3.117e-05; mean difference: -1148.499



#-----------------------------------------------------------------------------#

#colon 

#colon 
colon <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_colon_cell_annotation_11_26_24.rds")

filt_barcodes <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/colon/low_res/barcodes/filtered_low_res_barcodes_10_29_24.csv', what = "", sep = ",", skip = 1)

barcodes <- colnames(colon)

unstable <- setdiff(barcodes, filt_barcodes)

unstable_cells <- subset(colon, cells = unstable)

stable_cells <- subset(colon, cells = filt_barcodes)

#compute mean mt percentage per sample - stably assigned cell

# Extract metadata
meta_stable <- stable_cells@meta.data

# Compute mean mitochondrial percentage per sample
mean_mt_per_sample_stable <- meta_stable %>%
  group_by(orig.ident) %>%  # replace with your sample column name
  summarise(mean_percent_mt_stable = mean(percent.mt, na.rm = TRUE))




#t test
t.test(merge_mt_percentage_filt$mean_percent_mt_unstable, merge_mt_percentage_filt$mean_percent_mt_stable_filt, paired = TRUE)

#t = 7.1638, df = 35, p-value = 2.351e-08; mean difference: 1.797277

#nFeature_RNA

# Compute mean nFeature per sample
mean_feature_per_sample_stable_filt <- meta_stable_filt %>%
  group_by(orig.ident) %>%  # replace with your sample column name
  summarise(mean_feature_stable_filt = mean(nFeature_RNA, na.rm = TRUE))

# Compute mean nFeature per sample
mean_feature_per_sample_unstable <- meta_unstable %>%
  group_by(orig.ident) %>%  # replace with your sample column name
  summarise(mean_feature_unstable = mean(nFeature_RNA, na.rm = TRUE))


#merge mean nFeature per sample 
merge_feature_filt <- merge(mean_feature_per_sample_unstable, mean_feature_per_sample_stable_filt, by = "orig.ident")

#t test
t.test(merge_feature_filt$mean_feature_unstable, merge_feature_filt$mean_feature_stable_filt, paired = TRUE)

#t = -10.829, df = 35, p-value = 1.014e-12; mean difference: -325.6369

#nCount_RNA

# Compute mean nCount per sample
mean_count_per_sample_stable_filt <- meta_stable_filt %>%
  group_by(orig.ident) %>%  # replace with your sample column name
  summarise(mean_count_stable_filt = mean(nCount_RNA, na.rm = TRUE))

# Compute mean nCount per sample
mean_count_per_sample_unstable <- meta_unstable %>%
  group_by(orig.ident) %>%  # replace with your sample column name
  summarise(mean_count_unstable = mean(nCount_RNA, na.rm = TRUE))


#merge mean nCount per sample 
merge_count_filt <- merge(mean_count_per_sample_unstable, mean_count_per_sample_stable_filt, by = "orig.ident")

#t test
t.test(merge_count_filt$mean_count_unstable, merge_count_filt$mean_count_stable_filt, paired = TRUE)

#t = -11.191, df = 35, p-value = 4.12e-13; mean difference: -2312.415


#-----------------------------------------------------------------------------#

#rectum 

rectum <- readRDS(file = "/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/seurat_object/helm_batch1_13_rectum_cell_annotation_12_8_24.rds")

#load in stably assigned cells 
filtered_values <- scan('/Users/swashburn30/Desktop/Helmsley_Project/Batch1_13/stability_assessment/rectum/low_res/barcodes/filtered_low_res_barcodes_11_4_24.csv', what = "", sep = ",", skip = 1)

barcodes <- colnames(rectum)

unstable <- setdiff(barcodes, filtered_values)

unstable_cells <- subset(rectum, cells = unstable)

stable_cells <- subset(rectum, cells = filtered_values)

#compute mean mt percentage per sample - stably assigned cell

# Extract metadata
meta_stable <- stable_cells@meta.data

# Compute mean mitochondrial percentage per sample
mean_mt_per_sample_stable <- meta_stable %>%
  group_by(orig.ident) %>%  # replace with your sample column name
  summarise(mean_percent_mt_stable = mean(percent.mt, na.rm = TRUE))



#compute mean mt percentage per sample - unstably assigned cell 

# Extract metadata
meta_unstable <- unstable_cells@meta.data

# Compute mean mitochondrial percentage per sample
mean_mt_per_sample_unstable <- meta_unstable %>%
  group_by(orig.ident) %>%  # replace with your sample column name
  summarise(mean_percent_mt_unstable = mean(percent.mt, na.rm = TRUE))

mean_mt_per_sample_unstable


#merge mean mt per sample 
merge_mt_percentage_filt <- merge(mean_mt_per_sample_unstable, mean_mt_per_sample_stable_filt, by = "orig.ident")


#t test
t.test(merge_mt_percentage_filt$mean_percent_mt_unstable, merge_mt_percentage_filt$mean_percent_mt_stable_filt, paired = TRUE)

#t = 3.7372, df = 31, p-value = 0.0007539; mean difference: 1.269651

#statistically different 

#nFeature_RNA

# Compute mean nFeature per sample
mean_feature_per_sample_stable_filt <- meta_stable_filt %>%
  group_by(orig.ident) %>%  # replace with your sample column name
  summarise(mean_feature_stable_filt = mean(nFeature_RNA, na.rm = TRUE))


mean_feature_per_sample_unstable <- meta_unstable %>%
  group_by(orig.ident) %>%  # replace with your sample column name
  summarise(mean_feature_unstable = mean(nFeature_RNA, na.rm = TRUE))

#merge mean nFeature per sample 
merge_feature_filt <- merge(mean_feature_per_sample_unstable, mean_feature_per_sample_stable_filt, by = "orig.ident")


#t test
t.test(merge_feature_filt$mean_feature_unstable, merge_feature_filt$mean_feature_stable_filt, paired = TRUE)

#t = -5.7906, df = 31, p-value = 2.229e-06; mean difference: -292.8111 

#nCount_RNA

# Compute mean nCount per sample
mean_count_per_sample_stable_filt <- meta_stable_filt %>%
  group_by(orig.ident) %>%  # replace with your sample column name
  summarise(mean_count_stable_filt = mean(nCount_RNA, na.rm = TRUE))


mean_count_per_sample_unstable <- meta_unstable %>%
  group_by(orig.ident) %>%  # replace with your sample column name
  summarise(mean_count_unstable = mean(nCount_RNA, na.rm = TRUE))

#merge mean nCount per sample 
merge_count_filt <- merge(mean_count_per_sample_unstable, mean_count_per_sample_stable_filt, by = "orig.ident")


#t test
t.test(merge_count_filt$mean_count_unstable, merge_count_filt$mean_count_stable_filt, paired = TRUE)

#t = -7.1768, df = 31, p-value = 4.55e-08; mean difference: -1640.899  



