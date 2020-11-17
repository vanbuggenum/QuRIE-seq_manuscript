# I will create datatable with all filenames, custom samplenames(sampleIDs), and sample info seperated in multiple columns

# read all file locations from data/raw/ (dirsample), list original sample name (sample_origname), manually add if protein or RNA data (dataorigin)
samples_df <- data.frame(dirsample = paste0("data/raw/", list.files(path = "data/raw/")), sample_origname = list.files(path = "data/raw/"),dataorigin = c(rep(c("RNA","protein"),11)))

# Extract sampleid_brief from original sample_origname, manually add formatted condition per sample (usefull later for ordered orig.ident etc in seurat)
samples_df <- data.frame(samples_df,
                         sampleid_brief = separate(samples_df,col = "sample_origname", c("A", NA))[,"A"],
                         condition = c(rep(c("060.aIg.contr","180.aIg.contr","180.aIg.ibr","000.aIg.contr","002.aIg.contr","004.aIg.contr","006.aIg.contr","006.aIg.ibr","010.aIg.contr","030.aIg.contr","030.aIg.ibr"), each = 2)))
# Extract from condition metadata per sample (time, stimulus, inhibitor columns) also usefull for later seurat metadata, add sampleIDs that will be used in R
samples_df <- data.frame(samples_df,
                         time = separate(samples_df,col = "condition", c("A", "B", "C"))[,"A"],
                         stimulus = separate(samples_df,col = "condition", c("A", "B", "C"))[,"B"],
                         inhibitor = separate(samples_df,col = "condition", c("A", "B", "C"))[,"C"]) %>%
  mutate(sampleIDs = paste0(dataorigin, "_", condition))

write.csv2(samples_df, file = "output/metadata.csv")

RNAsamples <- filter(samples_df, str_detect(sampleIDs, "RNA"))
proteinsamples <-    filter(samples_df, str_detect(sampleIDs, "protein"))

# A for loop that reads all files defined in the samples dataframes. Including renaming of cellIDs:

##for proteins
for (i in 1:nrow(proteinsamples)){
  assign(proteinsamples$sampleIDs[i], read.table(as.character(proteinsamples$dirsample[i]), header=T, row.names=1, sep="\t") )

  assign(proteinsamples$sampleIDs[i],
         get(proteinsamples$sampleIDs[i])  %>%
           rename_all(funs(str_replace(., ".*prottrim_fully_", ""))) %>%
           rename_all(funs(str_replace(., ., paste0(proteinsamples$condition[i], "_",.)))) %>%
           assign(proteinsamples$sampleIDs[i],.)
  )

}


##for RNA (needs transform cols to rows)
for (i in 1:nrow(RNAsamples)) {
  assign(RNAsamples$sampleIDs[i], as.data.frame(t(read.table(as.character(RNAsamples$dirsample[i]), header=T, row.names=1, sep="\t"))) )

  assign(RNAsamples$sampleIDs[i],
         get(RNAsamples$sampleIDs[i])  %>%
           # rename_all(funs(str_replace(., ".*prottrim_fully_", ""))) %>% # remove extra info
           rename_all(funs(str_replace(., ., paste0(RNAsamples$condition[i], "_",.)))) # add sampleconditionbeforecellID
  )
}

## Combine imported counts into 2 tables & remove nonmeasured genes

all_RNA <- data.frame(row.names = 1:nrow(RNA_000.aIg.contr))
for(i in 1:nrow(RNAsamples)) {all_RNA <- cbind(all_RNA,get(RNAsamples$sampleIDs[i]))}
rownames(all_RNA) <- rownames(get(RNAsamples$sampleIDs[i]))

all_RNA <- all_RNA[!!rowSums(abs(all_RNA[-c(1:2)])),] # remove notmeasured genes

all_protein <- data.frame(row.names = 1:nrow(protein_000.aIg.contr))
for(i in 1:nrow(proteinsamples)) {all_protein <- cbind(all_protein,get(proteinsamples$sampleIDs[i]))}
rownames(all_protein) <- rownames(get(proteinsamples$sampleIDs[i]))

## Filter common cells for full dataset
commoncellIDs <- intersect(colnames(all_RNA), colnames(all_protein))
common_RNA <- all_RNA[, commoncellIDs]

common_protein <- all_protein[, commoncellIDs]

## Create seurat object with:
## genes detected in at least 100 cells,
## and cells with at least 100 genes

seu_combined <- CreateSeuratObject(counts = common_RNA, project = "common_cells", assay = "RNA", min.cells = 100, min.features = 100)
seu_combined[["PROT"]] <- CreateAssayObject(counts = common_protein)

## Prepare metadata information for seuratobject
meta.all <- select(RNAsamples, c("condition", "time", "stimulus", "inhibitor")) %>%
  mutate(orig.ident = condition)
meta.all <- left_join(seu_combined[[]], meta.all,by = "orig.ident")
rownames(meta.all) <- rownames(seu_combined[[]])
## Add metadata to seuratobject
seu_combined <- AddMetaData(seu_combined, metadata = meta.all)

## Calculate metadata properties (percentages features)
seu_combined[["percent.mt"]] <- PercentageFeatureSet(seu_combined, pattern = "^MT\\.", assay = "RNA")
seu_combined[["percent.rb"]] <- PercentageFeatureSet(seu_combined, pattern = "^RPL|^RPS", assay = "RNA")

seu_combined[["percent.Ig"]] <- PercentageFeatureSet(seu_combined, pattern = "^Ig", assay = "PROT")
seu_combined[["percent.HisH3"]] <- PercentageFeatureSet(seu_combined, pattern = "^His", assay = "PROT")
seu_combined[["percent.phospho.cum"]] <- PercentageFeatureSet(seu_combined, pattern = "^p-", assay = "PROT")

saveRDS(seu_combined, file = "output/seu_combined_raw.rds")
