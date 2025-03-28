suppressMessages(library(data.table))
# suppressMessages(library(VennDiagram))

run_all <- function(args){
  var_file <- args[1]
  output_file <- args[2]
  min_variant_frequency <- as.numeric(args[3]) / 100
  min_callers_threshold <- as.numeric(args[4])
  min_var_reads_threshold <- as.numeric(args[5])
  tmp_dir <- args[6]

  if(file.info(var_file)$size != 0){
    vcf <- vcfR::read.vcfR(var_file,verbose = F)

    var_tab <- data.table(chrom = vcf@fix[,"CHROM"]
                          ,position = as.integer(vcf@fix[,"POS"])
                          ,reference = vcf@fix[,"REF"]
                          ,alternative = vcf@fix[,"ALT"]
                          ,filter = vcf@fix[,"FILTER"]
                          ,genotype = vcfR::extract.gt(vcf)[,1]
                          ,var_reads = vcfR::extract.gt(vcf,element = "AD",as.numeric = T)[,1]
                          ,coverage_depth = vcfR::extract.gt(vcf,element = "DP",as.numeric = T)[,1])

    var_tab[is.na(coverage_depth),coverage_depth := vcfR::extract.gt(vcf,element = "DPI",as.numeric = T)[,1][!is.na(vcfR::extract.gt(vcf,element = "DPI",as.numeric = T)[,1])]]

    # ADD INFO ABOUT CALLERS CALLERS
    var_tab[!is.na(vcfR::extract.gt(vcf,element = "GQ")[,1]),caller := "haplotypecaller"]
    var_tab[!is.na(vcfR::extract.gt(vcf,element = "ALD")[,1]),caller := "vardict"]
    var_tab[!is.na(vcfR::extract.gt(vcf,element = "GQX")[,1]),caller := "strelka"]
    var_tab[!is.na(vcfR::extract.gt(vcf,element = "ABQ")[,1]),caller := "varscan"]

    #add and correct info
    var_tab[caller == "haplotypecaller",filter := "PASS"]
    var_tab[caller != "varscan",var_reads := coverage_depth - var_reads]
    var_tab[caller == "vardict",var_reads := vcfR::extract.gt(vcf,element = "VD",as.numeric = T)[which(var_tab$caller == "vardict"),1]]
    var_tab[var_reads > coverage_depth,var_reads := coverage_depth]
    if(any(var_tab$caller == "varscan",na.rm = T)){
      ADF <- vcfR::extract.gt(vcf,element = "ADF",as.numeric = T)[var_tab$caller == "varscan"]
      RDF <- vcfR::extract.gt(vcf,element = "RDF",as.numeric = T)[var_tab$caller == "varscan"]
      ADR <- vcfR::extract.gt(vcf,element = "ADR",as.numeric = T)[var_tab$caller == "varscan"]
      RDR <- vcfR::extract.gt(vcf,element = "RDR",as.numeric = T)[var_tab$caller == "varscan"]

      var_tab[var_tab$caller == "varscan",strand_bias := boot::inv.logit(ADF * (ADR + RDR) / (ADR * (ADF + RDF)))]

    } else {
      var_tab[,strand_bias := NA]
    }

    var_tab[,index := seq_along(var_tab$chrom)]


    #FILTER WIERD CALLS
    var_tab <- unique(var_tab,by = c("chrom","position","reference","alternative","caller"))
    var_tab <- var_tab[coverage_depth > 0]
    var_tab <- var_tab[!(is.na(alternative) | is.na(genotype) | is.na(filter))]
    var_tab <- var_tab[alternative != "*"]



    #ADD PASS COLLUMN
    var_tab[,is_pass := filter == "PASS" | filter == "f25"]
    var_tab[,is_pass := is_pass & var_reads >= min_var_reads_threshold]

    # GET CALLING STATISTICS
    stat_tab <- var_tab[,list(raw_variants = .N,pass_filter = sum(is_pass),ratio = round(sum(is_pass) / .N,3)),by = c("caller")]

    fwrite(stat_tab,file = gsub(".tsv$",".variant_stats.tsv",output_file),sep = "\t",col.names = T)

    if(length(unique(var_tab$caller)) > 1){
      caller_types <- unique(var_tab$caller)
      venn_list <- lapply(caller_types,function(x) var_tab[caller == x,paste(chrom,position,alternative,sep = "_")])
      venn_raw <- venn.diagram(venn_list,NULL,disable.logging = TRUE,
                               ,imagetype = "tiff"
                               ,fill=c("red", "green","blue","yellow")[seq_along(caller_types)]
                               ,alpha=c(0.5,0.5,0.5,0.5)[seq_along(caller_types)]
                               ,cex = 1.5
                               ,cat.fontface=2
                               ,category.names=caller_types,main = "Variants all")

      venn_raw_rel <- venn.diagram(venn_list,NULL,print.mode = "percent"
                                   ,imagetype = "tiff"
                                   ,fill=c("red", "green","blue","yellow")[seq_along(caller_types)]
                                   ,alpha=c(0.5,0.5,0.5,0.5)[seq_along(caller_types)]
                                   ,cex = 1.5
                                   ,cat.fontface=2
                                   ,category.names=caller_types,main = "Variants all relative")


      venn_list_pass <- lapply(caller_types,function(x) var_tab[caller == x & is_pass == T,paste(chrom,position,alternative,sep = "_")])
      venn_pass <- venn.diagram(venn_list_pass,NULL
                                ,imagetype = "tiff"
                                ,fill=c("red", "green","blue","yellow")[seq_along(caller_types)]
                                ,alpha=c(0.5,0.5,0.5,0.5)[seq_along(caller_types)]
                                ,cex = 1.5
                                ,cat.fontface=2
                                ,category.names=caller_types,main = "Variants pass")

      venn_pass_rel <- venn.diagram(venn_list_pass,NULL,print.mode = "percent"
                                    ,imagetype = "tiff"
                                    ,fill=c("red", "green","blue","yellow")[seq_along(caller_types)]
                                    ,alpha=c(0.5,0.5,0.5,0.5)[seq_along(caller_types)]
                                    ,cex = 1.5
                                    ,cat.fontface=2
                                    ,category.names=caller_types,main = "Variants pass relative")

      curent_wd <- getwd()
      setwd(tmp_dir)

      pdf(gsub(".tsv$",".variant_stats.pdf",output_file))
      grid.draw(venn_raw)
      grid.newpage()
      grid.draw(venn_pass)
      grid.newpage()
      grid.draw(venn_raw_rel)
      grid.newpage()
      grid.draw(venn_pass_rel)
      dev.off()

      file.remove(list.files(".",pattern = "VennDiagram.*.log"))

      setwd(curent_wd)
    }

    #filter all not pass vars
    var_tab <- var_tab[is_pass == T]
    var_tab[,is_pass := NULL]

    #join variants from callers and filter for min_callers_threshold,min_variant_frequency and min_variant_frequency
    var_tab[,call_info := paste0(var_reads,",",coverage_depth)]
    setorder(var_tab,caller)
    var_tab <- var_tab[,c("caller_count", "callers","all_callers_info") := list(.N,paste(caller,collapse = ","),paste(call_info,collapse = ";")),by = c("chrom","position","reference","alternative")]
    setorder(var_tab,-var_reads)
    var_tab <- unique(var_tab,by = c("chrom","position","reference","alternative"))

    var_tab <- var_tab[caller_count >= min_callers_threshold,]
    var_tab <- var_tab[ var_reads / coverage_depth >= min_variant_frequency,]
    var_tab <- var_tab[ var_reads >= min_variant_frequency,]

    #print vcf with just filtered variants
    vcf_out <- vcf[var_tab$index]
    vcfR::write.vcf(vcf_out,gsub(".tsv$",".vcf",output_file))

    #normalize deletions and insertions style
    if(nrow(var_tab) > 0){
      var_tab[,new_ref := reference]
      var_tab[sapply(seq_along(reference),function(x) grepl(reference[x],alternative[x])),new_ref := "-"]
      var_tab[sapply(seq_along(reference),function(x) grepl(alternative[x],reference[x])),new_ref := stringi::stri_sub(new_ref,from = nchar(alternative) + 1)]

      var_tab[,new_alt := alternative]
      var_tab[sapply(seq_along(reference),function(x) grepl(alternative[x],reference[x])),new_alt := "-"]
      var_tab[sapply(seq_along(reference),function(x) grepl(reference[x],alternative[x])),new_alt := stringi::stri_sub(new_alt,from = nchar(reference) + 1)]

      var_tab[,new_pos := position]
      var_tab[sapply(seq_along(reference),function(x) grepl(reference[x],alternative[x])),new_pos := new_pos + nchar(reference)]
      var_tab[sapply(seq_along(reference),function(x) grepl(alternative[x],reference[x])),new_pos := new_pos + nchar(alternative)]

      var_tab[,var_name := paste0(chrom,"_",new_pos,"_",new_ref,"/",new_alt)]
    } else {
      var_tab[,var_name := chrom]
    }

    #print as table
    var_tab[,variant_freq := round(var_reads / coverage_depth,5)]
    var_tab <- var_tab[,.(var_name,genotype,variant_freq,coverage_depth,caller_count,callers,strand_bias)]
    write.table(var_tab,file = output_file,sep = "\t",row.names = F,col.names = T,quote = F,na = "")
  } else {
    cols <- c("var_name","genotype","variant_freq","coverage_depth","caller_count","callers","strand_bias")
    var_tab <- data.table(matrix(vector(), nrow = 0, ncol = length(cols), dimnames = list(NULL, cols)))
    write.table(var_tab,file = output_file,sep = "\t",row.names = F,col.names = T,quote = F,na = "")
    system(paste0("touch ",gsub(".tsv$",".vcf",output_file)))
  }


}

# develop and test
# setwd("/Volumes/lamb/shared/S3acgt/sequia/11598__germline_small_var_call__var_call__230908/")
# args <- c("germline_varcalls/M485_C3_E5_2/merged/raw_calls.vcf","germline_varcalls/M485_C3_E5_2.final_variants.tsv","0","1","20","/tmp")

#run as Rscript
# 
# args <- commandArgs(trailingOnly = T)
# run_all(args)
