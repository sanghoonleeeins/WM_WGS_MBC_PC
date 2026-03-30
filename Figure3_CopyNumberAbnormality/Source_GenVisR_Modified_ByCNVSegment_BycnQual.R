library(ggplot2)
library(circlize) # rand_color(n=4)  To choose random colors.

dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir); print(dir) # 

multi_selectOut <- function(data, plot, out="plot", draw="FALSE")
{
  # Decide what to output
  if(toupper(out) == "DATA")
  {
    return(data)
  } else if(toupper(out) == "PLOT" & isTRUE(draw)) {
    return(grid::grid.draw(plot))
  } else if(toupper(out) == "PLOT" & !isTRUE(draw)) {
    return(plot)
  } else if(toupper(out) == "GROB" & isTRUE(draw)) {
    return(plot)
  } else if(toupper(out) == "GROB" & !isTRUE(draw)) {
    return(ggplot2::ggplotGrob(plot))
  } else {
    warning("Did not recognize input to out...")
    if(isTRUE(draw))
    {
      return(grid::grid.draw(plot))
    } else {
      return(plot)
    }
  }
}


cnFreq_disjoin <- function(x){

  # create the Granges object for the data
  x <- GenomicRanges::GRanges(seqnames=x$chromosome,
                              ranges=IRanges::IRanges(start=x$start, end=x$end),
                              "sample"=x$sample, "segmean"=x$segmean)

  # disjoin with grange, get a mapping of meta columns and expand it
  disJoint_x <- GenomicRanges::disjoin(x, with.revmap=TRUE)
  revmap <- GenomicRanges::mcols(disJoint_x)$revmap
  disJoint_x <- rep(disJoint_x, lengths(revmap))


  # exract the meta columns and map them back to the disJoint GRanges object
  sample <- unlist(IRanges::extractList(GenomicRanges::mcols(x)$sample, revmap))
  segmean <- unlist(IRanges::extractList(GenomicRanges::mcols(x)$segmean, revmap))
  GenomicRanges::mcols(disJoint_x)$sample <- sample
  GenomicRanges::mcols(disJoint_x)$segmean <- segmean

  # convert the GRanges Object back to a data frame
  disJoint_x <- as.data.frame(disJoint_x)[,c("seqnames", "start", "end", "width",
                                             "sample", "segmean")]
  colnames(disJoint_x) <- c("chromosome", "start", "end", "width", "sample", "segmean")
  return(disJoint_x)
}


multi_chrBound <- function(x) {
  # Check that input has size
  if(nrow(x) < 1)
  {
    memo <- paste0("input has 0 rows, it is possible that the UCSC",
                   " MySQL query has failed")
    stop(memo)
  }

  # Extract the columns needed
  data <- x[,c('chrom' ,'chromStart' , 'chromEnd')]

  # Obtain max for each chromosome
  maxChrom <- stats::aggregate(chromEnd ~ chrom, data=data, max)
  maxChrom <- cbind(maxChrom, maxChrom[,2])
  colnames(maxChrom) <- c('chromosome', 'start', 'end')

  # Obtain min for each chromosome
  minChrom <- stats::aggregate(chromStart ~ chrom, data=data, min)
  minChrom <- cbind(minChrom, minChrom[,2])
  colnames(minChrom) <- c('chromosome', 'start', 'end')

  # bind all the data together
  data <- rbind(maxChrom, minChrom)

  return(data)
}



cnFreq_qual <- function(x) {
  # Check that x is a data frame
  if(!is.data.frame(x)){
    memo <- paste0("Did not detect a data frame in argument supplied",
                   " to x... attempting to coerce")
    warning(memo)
    x <- as.data.frame(x)
    x <- droplevels(x)
  }

  # Check that x has at least 1 row
  if(nrow(x) < 1){
    memo <- paste0("x needs at least one row")
    stop(memo)
  }

  # remove any NA values in the data
  if(any(is.na(x))){
    na_rows_removed <- nrow(x) - nrow(na.omit(x))
    memo <- paste0("Removing ", na_rows_removed, " rows containing NA values")
    message(memo)
    x <- na.omit(x)
  }

  if(all(c('chromosome', 'start','end', 'segmean', 'sample') %in% colnames(x))){

    # make sure columns are of the correct type
    x$chromosome <- as.factor(x$chromosome)
    x$start <- as.integer(as.character(x$start))
    x$end <- as.integer(as.character(x$end))
    x$segmean <- as.numeric(as.character(x$segmean))
    x$sample <- as.factor(x$sample)

    # make sure windows are consistent if not disjoin them
    tmp <- split(x, x$sample)
    tmp_vec <- tmp[[1]]$end
    if(any(!unlist(sapply(tmp, function(x) x[,"end"] %in% tmp_vec), use.names=F))){
      memo <- paste0("Did not detect identical genomic segments for all samples",
                     " ...Performing disjoin operation")
      message(memo)

      # here we split the DF up in an attempt to avoid complaints that lists are to large
      x <- split(x, f=x$chromosome)
      x <- lapply(x, cnFreq_disjoin)


      x <- do.call(rbind, x)
    }
    rm(tmp)
    rm(tmp_vec)
  } else {
    memo <- paste0("Did not detect correct columns in argument supplied",
                   " to x!")
    stop(memo)
  }

  # Check chromosome column in x
  if(!all(grepl("^chr", x$chromosome))){
    memo <- paste0("Did not detect the prefix \"chr\" in the chromosome",
                   " column of x... adding prefix")
    message(memo)
    x$chromosome <- paste0("chr", x$chromosome)
    x$chromosome <- as.factor(x$chromosome)
  } else if(all(grepl("^chr", x$chromosome))) {
    memo <- paste0("Detected \"chr\" in the chromosome column of x...",
                   " proceeding")
    message(memo)
  } else {
    memo <- paste0("Detected unknown or mixed prefixes in the chromosome ",
                   " column of x, should either have a chr prefix or ",
                   "none at all!")
    stop(memo)
  }

  return(x)
}


cnData_Chr6qOnly <- readRDS("cnData_Chr6qOnly.rds")
x=cnData_Chr6qOnly; WMSubtype="AllType"
cnData_NoWM91WM93WM94_MBCLike <- readRDS("cnData_NoWM91WM93WM94_MBCLike.rds")

x=cnData_NoWM91WM93WM94_MBCLike; WMSubtype="MBCLike"    #  x=cnData_NoWM91WM93WM94_PCLike; WMSubtype="PCLike"  # x=cnData_NoWM91WM93WM94_NotAvail; WMSubtype="NotAvail"
CN_low_cutoff=1.5;CN_high_cutoff=2.5;plot_title=NULL; CN_Loss_colour="#002EB8";CN_Gain_colour="#A30000";
x_title_size=12; y_title_size=12; facet_lab_size=10; plotType="proportion"; out="plot"
genome="hg38"; plotChr = "chr6";

################################### ++++++++++++++++++++++  cnFreq() 
cnFreq_Modified <- function (x, CN_low_cutoff=1.5, CN_high_cutoff=2.5, plot_title=NULL,  CN_Loss_colour="#002EB8", CN_Gain_colour="#A30000", 
                             x_title_size=12,  y_title_size=12, facet_lab_size=10, plotLayer, plotType="proportion", 
                             out="plot", genome="hg38",plotChr=plotChr,  WMSubtype, ChrNumb, CentromereStart, ChrStartPosition, ChrEndPosition) {
    ## I added this line
    OriginalX <- x; dim(OriginalX) # 32 17 # NotAvail: 18 17  ## All Chr01q: 159 17
  
    x <- cnFreq_qual(x); dim(x) # 11461 6 # PC: 262 6  #  length(unique(x$sample)): 47  Nothing lost. 
    cnFreqQual <- x
  
    samples <- unique(x$sample); length(samples) # 47 # [1] 01-190 WM33   WM70   WM47
    gainFreq <- function(x) {
      length(x[x >= CN_high_cutoff])
    }
    # gainFrequency <- aggregate(segmean ~ chromosome + start +  end, data = x, gainFreq)$segmean
    gainFrequency <- aggregate(segmean ~ chromosome + start +  end, data = x, gainFreq)$segmean
    lossFreq <- function(x) {
      length(x[x <= CN_low_cutoff])
    }
    lossFrequency <- aggregate(segmean ~ chromosome + start + end, data = x, lossFreq)$segmean
    x <- aggregate(segmean ~ chromosome + start + end, data = x, length); dim(x) # 63 4
    colnames(x)[which(colnames(x) %in% "segmean")] <- "sampleFrequency"
    x$gainFrequency <- gainFrequency
    x$lossFrequency <- lossFrequency
    if (max(x$sampleFrequency) > length(samples)) {
        memo <- paste0("Detected additional sample rows after disjoin operation",
                       " typically this indicates coordinates are 0-based, please convert",
                       " coordinates to 1-base for accurate results")
        warning(memo)
    }
    x$gainProportion <- x$gainFrequency/length(samples)
    x$lossProportion <- x$lossFrequency/length(samples); dim(x) # 63 8
    preloaded <- c("hg38", "hg19", "mm10", "mm9", "rn5")
    if (any(genome == preloaded)) {
        message("genome specified is preloaded, retrieving data...")
        UCSC_Chr_pos <- GenVisR::cytoGeno[GenVisR::cytoGeno$genome == genome, ]
        UCSC_Chr_pos <- multi_chrBound(UCSC_Chr_pos)
    } else {
        memo <- paste0("attempting to query UCSC mySQL database for chromosome",
                       " positions")
        message(memo)
        cyto_data <- suppressWarnings(multi_cytobandRet(genome))
        UCSC_Chr_pos <- multi_chrBound(cyto_data)
    }
  
    if (nrow(UCSC_Chr_pos) < 1) {
        memo <- paste0("did not recognize genome ", genome, ", plotting provided data and ignoring chromosome ",
                       "boundaries! Output could be decieving!")
        warning(memo)
    }
  
    y <- x %>% dplyr::mutate(chromosome=gsub("chr","", chromosome))
  
    dummy_data <- lapply(unique(y$sample), function(sample, chr_pos) cbind(chr_pos,  sample), UCSC_Chr_pos);
  
    dummy_data <- do.call("rbind", dummy_data)
    chr_order <- gtools::mixedsort(unique(dummy_data$chromosome))
    dummy_data$chromosome <- factor(dummy_data$chromosome, levels=chr_order)
    if (!is.null(plotChr)) {
        if (any(!plotChr %in% dummy_data$chromosome)) {
            missingChr <- plotChr[!plotChr %in% dummy_data$chromosome]
            plotChr <- plotChr[!plotChr %in% missingChr]
            memo <- paste0("The following chromosomes: ", toString(missingChr),
                           ", could not be found! Valid chromosomes are: ",
                           toString(unique(dummy_data$chromosome)))
            warning(memo)
        }
        dummy_data <- dummy_data[dummy_data$chromosome %in% plotChr,]
        dummy_data$chromosome <- factor(dummy_data$chromosome, levels = plotChr)
    
        x <- x[x$chromosome %in% plotChr, ]
        x$chromosome <- factor(x$chromosome, levels = plotChr)
    }

    ## p1 is barplot.   centromere vertical lines are made in cnFreq_buildMain_Modified()
    p1 <- cnFreq_buildMain_Modified(x=x, plotType, dummy_data = dummy_data, ############# <<<<<<======== I modified cnFreq_buildMain_Modified() 
                                        plot_title = plot_title, CN_low_colour = CN_Loss_colour,
                                        CN_high_colour = CN_Gain_colour, x_lab_size = x_title_size,
                                        y_lab_size = y_title_size, facet_lab_size = facet_lab_size,
                                        plotLayer = plotLayer, OriginalX=OriginalX, WMSubtype=WMSubtype, ChrNumb, 
                                        CentromereStart, ChrStartPosition, ChrEndPosition)  #### <<<<<<<<========= I added "WMSubtype=WMSubtype"
        # output <- multi_selectOut(data=list(data = x), plot = p1, out = out)
        # return(output)
    return(p1)
       
}


######### ++++++++++++ ######### ++++++++++++ ######### ++++++++++++ ######### ++++++++++++ ######### ++++++++++++ ######### ++++++++++++ ######### ++++++++++++ 
######### ++++++++++++ ######### ++++++++++++ ######### ++++++++++++ ######### ++++++++++++ ######### ++++++++++++ ######### ++++++++++++ ######### ++++++++++++ 

CN_low_colour=CN_Loss_colour; CN_high_colour=CN_Gain_colour; x_lab_size=x_title_size; y_lab_size=y_title_size


cnFreq_buildMain_Modified <- function(x, plotType, dummy_data, plot_title=NULL, CN_low_colour='#002EB8', CN_high_colour='#A30000',
                                      x_lab_size=12, y_lab_size=12, facet_lab_size=10, plotLayer=NULL, OriginalX=OriginalX,
                                      WMSubtype, ChrNumb, CentromereStart, ChrStartPosition, ChrEndPosition) {
  # Transform losses to be negative values for plotting purposes
  x$lossFrequency <- -1*x$lossFrequency
  x$lossProportion <- -1*x$lossProportion

  # Define parameters of plot
  theme <- ggplot2::theme(strip.text.x=element_text(size=facet_lab_size),
                 axis.text.x=ggplot2::element_blank(),
                 axis.ticks.x=ggplot2::element_blank(),
                 legend.position='right',
                 axis.title.x=ggplot2::element_text(size=x_lab_size, face='bold'),
                 axis.title.y=ggplot2::element_text(size=y_lab_size, face='bold'),
                 panel.grid.major.x=ggplot2::element_blank(),
                 panel.grid.minor.x=ggplot2::element_blank())
  facet <- facet_grid(. ~ chromosome, scales='free', space='free')
  xlabel <- xlab('Chromosomes')

  # Choose whether to plot aesthetics for proportion or frequency
  if(grepl("^PROP", plotType, ignore.case=TRUE)){
      ylabel <- ylab("Proportion of Copy Number Gains/Losses")
      ymax <- 1
      x$gain <- x$gainProportion
      x$loss <- x$lossProportion
  } else if(grepl("^FREQ", plotType, ignore.case=TRUE)){
      ylabel <- ylab("Frequency of Copy Number Gains/Losses")
      ymax <- max(as.numeric(as.character(x$sampleFrequency)), na.rm=TRUE)
      x$gain <- x$gainFrequency
      x$loss <- x$lossFrequency
  } else {
      memo <- paste0("did not recognize plotType ", plotType,
                     ", please specify one of \"proportion\" or \"frequency\"")
      stop(memo)
  }

  # Define the initial plot
  # p1 <- ggplot(data=dummy_data, mapping=aes_string(xmin='start', xmax='end', ymin=-1*ymax, ymax=ymax)) + geom_rect(alpha=0) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
  # p1 <- ggplot(data=dummy_data, mapping=aes(xmin=start, xmax=end, ymin=-1*ymax, ymax=ymax)) + geom_rect(alpha=0) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
  if(ChrNumb %in% c("Chr6q","Chr8p","Chr13q","Chr17p","Chr21q") ) {       #   <<<<<==== This is for LOSS either p or q arm
          p1 <- ggplot(data=dummy_data, mapping=aes(xmin=ChrStartPosition, xmax=ChrEndPosition, ymin=-1*ymax, ymax=0.75*ymax)) + geom_rect(alpha=0) + # 172805979  ####### <<<<<<===== I revised here. 
                  scale_x_continuous(expand=c(0,0), n.breaks=25) + scale_y_continuous(expand=c(0,0))  # xmin=59800001 # xmax=170805979
          InitialPlot <- p1  ### <<<<<<<=========   This is to get an empty plot frame and draw CNV segment horizontal lines
  } else if (ChrNumb %in% c("Chr1q","Chr6p","Chr8q","Chr9q","Chr13p","Chr17q","Chr21p") ) {   ###### <<<<<==== This is for GAIN either in p or q arm 
          print(paste0("ChrNumb: ", ChrNumb))
          p1 <- ggplot(data=dummy_data, mapping=aes(xmin=ChrStartPosition, xmax=ChrEndPosition, ymin=-0.5*ymax, ymax=ymax)) + geom_rect(alpha=0) +   ####### <<<<<<===== I revised here.  ### Chr1q: xmin=1.20*10^8, xmax=248956422
                  scale_x_continuous(expand=c(0,0), n.breaks=25) + scale_y_continuous(expand=c(0,0))  ### xmin should be lower than "y" data 'start' position in line 326
          InitialPlot <- p1  ### <<<<<<<=========   This is to get an empty plot frame and draw CNV segment horizontal lines
          # print(paste0("ChrNumb done: ", ChrNumb))
  } else if (ChrNumb %in% c("Chr1","Chr2","Chr3","Chr4", "Chr5",    "Chr6","Chr7","Chr8","Chr9","Chr10",     "Chr11","Chr12","Chr13","Chr14","Chr15",
                            "Chr16","Chr17","Chr18","Chr19","Chr20",    "Chr21","Chr22"))  {
         p1 <- ggplot(data=dummy_data, mapping=aes(xmin=ChrStartPosition, xmax=ChrEndPosition, ymin=-0.5*ymax, ymax=ymax)) + geom_rect(alpha=0) +   ####### <<<<<<===== I revised here.  
                  scale_x_continuous(expand=c(0,0), n.breaks=25) + scale_y_continuous(expand=c(0,0))  # xmin=59800001 # xmax=170805979
         InitialPlot <- p1  ### <<<<<<<=========   This is to get an empty plot frame and draw CNV segment horizontal lines
  }

  # add copy number data
  # p1 <- p1 + geom_rect(data=x, mapping=aes_string(xmin='start',xmax='end', ymin='loss', ymax=0), fill=CN_low_colour)
  # p1 <- p1 + geom_rect(data=x, mapping=aes_string(xmin='start', xmax='end', ymin=0, ymax='gain'), fill=CN_high_colour)
  
  if (ChrNumb %in% c("Chr6p","Chr8p","Chr13p","Chr17p","Chr21p")  ) {   ### <<<<+++++ This is 'p' arms
          y <-  x %>% dplyr::filter(end <= ChrEndPosition); dim(x); dim(y) #### I added this line one 20151212 to make GenVisR plot for only Chr6q or Chr1q    start>1.21*10^8
  } else if(ChrNumb %in% c("Chr1q","Chr6q","Chr8q","Chr9q","Chr13q","Chr17q","Chr21q" ) ) {  ### <<<<+++++ This is 'q' arms
          y <- x %>% dplyr::filter(start > ChrStartPosition); dim(x); dim(y) #### I added this line one 20151212 to make GenVisR plot for only Chr6q or Chr1q
  } else if (ChrNumb %in% c("Chr1","Chr2","Chr3","Chr4", "Chr5",    "Chr6","Chr7","Chr8","Chr9","Chr10",     "Chr11","Chr12","Chr13","Chr14","Chr15",
                            "Chr16","Chr17","Chr18","Chr19","Chr20",    "Chr21","Chr22")) {
          y <- x
  }
  
  p1 <- p1 + geom_rect(data=y, mapping=aes(xmin=start,xmax=end, ymin=1*loss, ymax=0), fill=CN_low_colour)
  p1 <- p1 + geom_rect(data=y, mapping=aes(xmin=start,xmax=end, ymin=0, ymax=gain), fill=CN_high_colour)
  p1 <- p1 + geom_hline(aes(yintercept=0), linetype="dotted")

  # build the plot to "p1".   CNV barplot is done hre. 
  p1 <- p1 + ylabel + xlabel + facet + theme_bw() # + theme  #### CNV plot is done here.  Sometimes, plot is not displayed. I don't know why. Just run from the beginning again.  

  ############# ++++++++ Indicate deletion size on the top of the CNV plot. +++++++++++++ #########
  ############# ++++++++ Make horizontal lines to indicate Deltion region +++++++++++++ #########
   ##### 1) Chr1q, Chr6p, Chr6q, Chr8p, Chr8q, Chr9q Chr13p, Chr13q, Chr17p, Chr17q, Chr21p, Chr21q, Chr22p
  
  if(ChrNumb %in% c("Chr6q","Chr8p","Chr13q","Chr17p","Chr21q") ) {  ##### <<<<<=== This is for LOSS
          ########### +++++++++++++++ using the original CNV Chr6q Del data before cnFreq_qual(X)  to find 6q Del regions +++++++++++++= #################
          OriginalX_DelSize <- OriginalX %>% dplyr::filter(segmean < 2, start>ChrStartPosition, end<ChrEndPosition)  %>%
                                  dplyr::mutate(DelSize=(end-start));  dim(OriginalX_DelSize) # 38 18  # MBC-like 18 18 # PC-like: 13 18
          ######  Use "59*10^6"  instead of 61*10^6  Because some samples have 6p_Partial~6p_WholeArm loss. 
          OriginalX_DelSize_Sort <- OriginalX_DelSize[with(OriginalX_DelSize, order(DelSize, decreasing=FALSE)),]
        
          ### This number will be the same number of horizontal lines. 
          StartX <- OriginalX_DelSize_Sort$start; length(StartX); StartX[1:3] ## All: Chr6q: 38  # MBC-like:18 #  PC: Chr6q:13, 170716900 170613752 170606038 # NotAvail: Chr1q:22 Chr6q:18
          EndX <- OriginalX_DelSize_Sort$end; length(EndX); EndX[1:3] # All: Chr6q: 38. MBC-like:18.   PC:13  # NotAvail: 18
        
          # NumberLine <- length(StartX); print(NumberLine) # AllType: 68 # 18
          # SpacePerLine <- (0.5-0.25)/NumberLine; print(SpacePerLine) # 0.015625  <<- make horizontal lines in the y-axis 0.5~0.8 range. # PC: 0.01388889
          # HorizontalPosition <- seq(from=0.5, to=(0.25+SpacePerLine), by=-SpacePerLine); length(HorizontalPosition) # AllType, 68  # MBC, 32 # PC: 18
        
          length(OriginalX_DelSize_Sort$sample); table(OriginalX_DelSize_Sort$sample);length(unique(OriginalX_DelSize_Sort$sample)) ####### <<<<<<========= The number of samples. This will determine the number colors for the horizontal line.. 
          ## AllSmp: Chr6q:38 # MBC-ike:3s # PC-like Chr6q: 12 # NotAvail:. Chr6q:8

          ########### +++++++++++++++ using cnFreq_qual(X)  to find 6q Del regions +++++++++++++= #################
          ## left_join between x and OriginalX[, "start", "sample"] by "start"
          # x_SmpName <- dplyr::left_join(x, OriginalX[,c("start","sample")]) %>% dplyr::filter(!is.na(sample)) %>% dplyr::mutate(DelSize=(end-start)) %>%
          #             dplyr::filter(loss!=0); dim(x_SmpName)
          # x_SmpName_Sort <- x_SmpName[with(x_SmpName, order(DelSize, decreasing=FALSE)),]
          #
          # StartX <- x_SmpName_Sort$start; length(StartX); StartX[1:3] # All: 66  # MBC: 25 # PC: 18, 170716900 170613752 170606038
          # EndX <- x_SmpName_Sort$end; length(EndX); EndX[1:3 ] # PC:18, 170744200 170744318 170742076
          # length(x_SmpName_Sort$sample); table(x_SmpName_Sort$sample); length(unique(x_SmpName_Sort$sample)) ##  All: 24
          
  } else if(ChrNumb %in% c("Chr1q","Chr6p","Chr8q","Chr9q","Chr13p","Chr17q","Chr21p","Chr22p")  ) {    ## <<<<==== This is for GAIN
          ########### +++++++++++++++ using the original CNV Chr1q Gain data before cnFreq_qual(X)  to find 1q Gain regions +++++++++++++= #################
          table(OriginalX$sample, OriginalX$segmean)
          # OriginalX_DelSize <- OriginalX %>% dplyr::filter(!sample %in% c("01-076","01-115","01-131","01-163","01-190","04-006",  ## filter out samples that don't have Chr1q Gain. 
          #                                         "WM9","WM37","WM42","WM46","WM47","WM48","WM52","WM54","WM65","WM66","WM67","WM69","WM70",
          #                                         "WM73","WM74","WM80","WM82","WM83","WM84","WM85","WM89","WM95","WM96","WM97")) %>% dplyr::filter(segmean > 2) %>%
          OriginalX_DelSize <- OriginalX %>% dplyr::filter(segmean>2, start>ChrStartPosition, end<ChrEndPosition) %>% 
                                                    dplyr::mutate(DelSize=(end-start));  dim(OriginalX_DelSize)# View(OriginalX_DelSize) # AllSmp: 11 18
          OriginalX_DelSize_Sort <- OriginalX_DelSize[with(OriginalX_DelSize, order(DelSize, decreasing=FALSE)),]
          length(OriginalX_DelSize_Sort$sample); table(OriginalX_DelSize_Sort$sample); length(unique(OriginalX_DelSize_Sort$sample))
          ####### <<<<<<========== The number of samples.AllSubtype:11, PC:7 This number should match the number of colors
          
          StartX <- OriginalX_DelSize_Sort$start; length(StartX); StartX[1:3] ## All: Chr1q:11 # MBC Chr1q:15, Chr6q: 32 #  PC: Chr1q:17,  Chr6q:18, 170716900 170613752 170606038 # NotAvail: Chr1q:22 Chr6q:18
          EndX <- OriginalX_DelSize_Sort$end; length(EndX); EndX[1:3] # All: Chr1q:11 PC:17 # NotAvail: 3
  } else if(ChrNumb %in% c("Chr1","Chr2","Chr3","Chr4", "Chr5",    "Chr6","Chr7","Chr8","Chr9","Chr10",     "Chr11","Chr12","Chr13","Chr14","Chr15",
                           "Chr16","Chr17","Chr18","Chr19","Chr20",    "Chr21","Chr22")  ) {
          table(OriginalX$sample, OriginalX$segmean)
          OriginalX_DelSize <- OriginalX %>% dplyr::filter(segmean>2) %>% 
                          dplyr::mutate(DelSize=(end-start));  dim(OriginalX_DelSize)  # View(OriginalX_DelSize) # AllSmp: 9 19 
          
          OriginalX_DelSize_Sort <- OriginalX_DelSize[with(OriginalX_DelSize, order(DelSize, decreasing=FALSE)),]
          length(OriginalX_DelSize_Sort$sample); table(OriginalX_DelSize_Sort$sample); length(unique(OriginalX_DelSize_Sort$sample)) # 9  7 
          
          if(nrow(OriginalX_DelSize_Sort)==0) {  ########### <<<<<<=============== I added this line when OriginalX_DelSize has no sample
                    OriginalX_DelSize_Sort[1,]$sample <- "Fake"
                    StartX <- 0; length(StartX); StartX[1:3] ## All: Chr1q:11 # MBC Chr1q:15, Chr6q: 32 #  PC: Chr1q:17,  Chr6q:18, 170716900 170613752 170606038 # NotAvail: Chr1q:22 Chr6q:18
                    EndX <- 0; length(EndX); EndX[1:3] # All: Chr1q:11 PC:17 # NotAvail: 3
          }  else { 
                StartX <- OriginalX_DelSize_Sort$start; length(StartX); StartX[1:3] ## All: Chr1q:11 # MBC Chr1q:15, Chr6q: 32 #  PC: Chr1q:17,  Chr6q:18, 170716900 170613752 170606038 # NotAvail: Chr1q:22 Chr6q:18
                EndX <- OriginalX_DelSize_Sort$end; length(EndX); EndX[1:3] # All: Chr1q:11 PC:17 # NotAvail: 3
                ####### <<<<<<==========  This number should match the number of colors # Chr3: AllType: 7 
          }
  }
          
  ### Write "OriginalX_DelSize" which has Chr#, WMSubtype SegmentStart, SegmentStop
  fwrite(OriginalX_DelSize, paste0("/Users/lees130/Library/CloudStorage/OneDrive-NYULangoneHealth/N05b_WaldenstromMacroglobulinemia_Dylan/04b_WM_PROJECT_INFO_WMMultiOmics_20250518/WGS/code/OriginalX_DelSize_", ChrNumb,"_",WMSubtype,"_v1.7.txt"),
         sep="\t", quote=FALSE, row.names=FALSE   )
  
  if( WMSubtype=="AllType") {
          NumberLine <- length(StartX); print(NumberLine) # Chr1q, AllType: 11  # MBC, 25  # PC, 17 # Chr6q, 38 
          if (ChrNumb %in% c("Chr6q","Chr8p","Chr13q","Chr17p","Chr21q" ) ) {   ## This is for LOSS
                      SpacePerLine <- (0.75-0.2)/NumberLine; print(SpacePerLine) #AllType  0.009558824  <<- make horizontal lines in the y-axis 0.5~0.8 range. # PC: 0.018
                      HorizontalPosition <- seq(from=0.75, to=(0.2+SpacePerLine), by=-SpacePerLine); length(HorizontalPosition) # AllType, 68  # MBC, 25 # PC: 17
                      
                      length(unique(OriginalX_DelSize_Sort$sample)) # 20  This will determine the number of colors. 
                      # OriginalX_DelSize_Sort <- x_SmpName_Sort
                      set.seed=2
                      # DelRegionLineColor <- rand_color(length(unique(OriginalX_DelSize_Sort$sample))); print(DelRegionLineColor) #  24   #### By the number of samples. 
                        DelRegionLineColor <- c("green","black","purple","darkorange","grey",   "red","darkgoldenrod","darkgreen","deeppink","skyblue",     "bisque3","aquamarine1","pink","darkkhaki",'yellow',
                                              "coral","cyan","darkgoldenrod1","azure","darkolivegreen1",     "darkblue","deeppink4","darkorchid4",'brown1')
                                            # 01-190            WM42            WM82            WM96            WM43            WM22            WM30            WM84            WM70            WM33            WM86            WM69            WM66            WM48 
                                            # "green"         "black"        "purple"    "darkorange"          "grey"           "red" "darkgoldenrod"     "darkgreen"      "deeppink"       "skyblue"       "bisque3"   "aquamarine1"          "pink"     "darkkhaki" 
                                            # WM68            WM73            WM89 
                                            # "yellow"         "coral"          "cyan" 
                        
          }  else if (ChrNumb %in% c("Chr1q","Chr6p","Chr8q","Chr9q","Chr13p","Chr17q","Chr21p","Chr22p",
                                     "Chr1","Chr2","Chr3","Chr4", "Chr5",    "Chr6","Chr7","Chr8","Chr9","Chr10",     "Chr11","Chr12","Chr13","Chr14","Chr15",
                                     "Chr16","Chr17","Chr18","Chr19","Chr20",    "Chr21","Chr22") ) {   ## This is for GAIN
                      # SpacePerLine <- (-0.95+0.3)/NumberLine; print(SpacePerLine) #AllType  0.009558824  <<- make horizontal lines in the y-axis 0.5~0.8 range. # PC: 0.018
                      # HorizontalPosition <- seq(from=-0.95, to=(-0.3-SpacePerLine), by=-SpacePerLine); length(HorizontalPosition) # AllType, 68  # MBC, 25 # PC: 17
                      
                      SpacePerLine <- (0.95-0.2)/NumberLine; print(SpacePerLine) #AllType  0.009558824  <<- make horizontal lines in the y-axis 0.5~0.8 range. # PC: 0.018
                      HorizontalPosition <- seq(from=0.95, to=(-0.2+SpacePerLine), by=-SpacePerLine); length(HorizontalPosition) # AllType, 68  # MBC, 25 # PC: 17
                      
                      # OriginalX_DelSize_Sort <- x_SmpName_Sort
                      set.seed=2
                      # DelRegionLineColor <- rand_color(length(unique(OriginalX_DelSize_Sort$sample))); print(DelRegionLineColor) #  24
            
                       DelRegionLineColor <- c("green","black","purple","darkorange","grey",   "red","darkgoldenrod","darkgreen","deeppink","skyblue",     "bisque3", "aquamarine1","pink","darkkhaki",'yellow',
                                              "coral","cyan" ,"darkgoldenrod1","azure","darkolivegreen1",     "darkblue","deeppink4","darkorchid4","brown1", "cyan1" )    
                                   
          } 
          
          #### Keep only necessary colors according to the number of unique samples, length(unique(OriginalX_DelSize_Sort$sample))
          DelRegionLineColor_CutByNumbSmp <- DelRegionLineColor[1:length(unique(OriginalX_DelSize_Sort$sample))]   ## I don't need all 11 colors defined in line 416 
          names(DelRegionLineColor_CutByNumbSmp)  <- unique(OriginalX_DelSize_Sort$sample); print(DelRegionLineColor)
          DelRegionLineColor_DF <- data.frame(DelRegionLineColor_CutByNumbSmp) %>% tibble::rownames_to_column("sample")
          fwrite(DelRegionLineColor_DF, file=paste0("ColorCodeBySmp_", ChrNumb,"_",WMSubtype,"_v1.7.txt"), col.names=TRUE,row.names=FALSE, sep="\t", quote=FALSE)
          
          ### inner_join between OriginalX_DelSize_Sort and DelRegionLineColor_DF to map the line color per sample ID
          OriginalX_DelSize_LineColor <- dplyr::left_join(OriginalX_DelSize_Sort, DelRegionLineColor_DF); dim(OriginalX_DelSize_LineColor)

          Myggplot <- InitialPlot  ##. Myggplot <- p1  #### If I want to include horizontal line plot into CNV plot. 
          ## Myggplot <- p1
          
          # if (ChrNumb=="Chr6q") { 
                    Myggplot<- Myggplot+
                      geom_segment(aes(x=c(StartX[1]), y=c(HorizontalPosition[1]), xend=c(EndX[1]), yend=c(HorizontalPosition[1])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[1], size=1.5) +
                      geom_segment(aes(x=c(StartX[2]), y=c(HorizontalPosition[2]), xend=c(EndX[2]), yend=c(HorizontalPosition[2])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[2], size=1.5) +
                      geom_segment(aes(x=c(StartX[3]), y=c(HorizontalPosition[3]), xend=c(EndX[3]), yend=c(HorizontalPosition[3])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[3], size=1.5) +
                      geom_segment(aes(x=c(StartX[4]), y=c(HorizontalPosition[4]), xend=c(EndX[4]), yend=c(HorizontalPosition[4])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[4], size=1.5) +
                      geom_segment(aes(x=c(StartX[5]), y=c(HorizontalPosition[5]), xend=c(EndX[5]), yend=c(HorizontalPosition[5])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[5], size=1.5) +
                      geom_segment(aes(x=c(StartX[6]), y=c(HorizontalPosition[6]), xend=c(EndX[6]), yend=c(HorizontalPosition[6])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[6], size=1.5) +
                      geom_segment(aes(x=c(StartX[7]), y=c(HorizontalPosition[7]), xend=c(EndX[7]), yend=c(HorizontalPosition[7])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[7], size=1.5) +
                      geom_segment(aes(x=c(StartX[8]), y=c(HorizontalPosition[8]), xend=c(EndX[8]), yend=c(HorizontalPosition[8])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[8], size=1.5) +
                      geom_segment(aes(x=c(StartX[9]), y=c(HorizontalPosition[9]), xend=c(EndX[9]), yend=c(HorizontalPosition[9])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[9], size=1.5) +
                      geom_segment(aes(x=c(StartX[10]), y=c(HorizontalPosition[10]), xend=c(EndX[10]), yend=c(HorizontalPosition[10])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[10], size=1.5) +
          
                      geom_segment(aes(x=c(StartX[11]), y=c(HorizontalPosition[11]), xend=c(EndX[11]), yend=c(HorizontalPosition[11])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[11], size=1.5) +
                      geom_segment(aes(x=c(StartX[12]), y=c(HorizontalPosition[12]), xend=c(EndX[12]), yend=c(HorizontalPosition[12])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[12], size=1.5) +
                      geom_segment(aes(x=c(StartX[13]), y=c(HorizontalPosition[13]), xend=c(EndX[13]), yend=c(HorizontalPosition[13])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[13], size=1.5) +
                      geom_segment(aes(x=c(StartX[14]), y=c(HorizontalPosition[14]), xend=c(EndX[14]), yend=c(HorizontalPosition[14])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[14], size=1.5) +
                      geom_segment(aes(x=c(StartX[15]), y=c(HorizontalPosition[15]), xend=c(EndX[15]), yend=c(HorizontalPosition[15])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[15], size=1.5) +
                      geom_segment(aes(x=c(StartX[16]), y=c(HorizontalPosition[16]), xend=c(EndX[16]), yend=c(HorizontalPosition[16])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[16], size=1.5) +
                      geom_segment(aes(x=c(StartX[17]), y=c(HorizontalPosition[17]), xend=c(EndX[17]), yend=c(HorizontalPosition[17])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[17], size=1.5) +
                      geom_segment(aes(x=c(StartX[18]), y=c(HorizontalPosition[18]), xend=c(EndX[18]), yend=c(HorizontalPosition[18])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[18], size=1.5) +
                      geom_segment(aes(x=c(StartX[19]), y=c(HorizontalPosition[19]), xend=c(EndX[19]), yend=c(HorizontalPosition[19])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[19], size=1.5) +
                      geom_segment(aes(x=c(StartX[20]), y=c(HorizontalPosition[20]), xend=c(EndX[20]), yend=c(HorizontalPosition[20])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[20], size=1.5) +

                      geom_segment(aes(x=c(StartX[21]), y=c(HorizontalPosition[21]), xend=c(EndX[21]), yend=c(HorizontalPosition[21])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[21], size=1.5) +
                      geom_segment(aes(x=c(StartX[22]), y=c(HorizontalPosition[22]), xend=c(EndX[22]), yend=c(HorizontalPosition[22])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[22], size=1.5) +
                      geom_segment(aes(x=c(StartX[23]), y=c(HorizontalPosition[23]), xend=c(EndX[23]), yend=c(HorizontalPosition[23])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[23], size=1.5) +
                      geom_segment(aes(x=c(StartX[24]), y=c(HorizontalPosition[24]), xend=c(EndX[24]), yend=c(HorizontalPosition[24])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[24], size=1.5) +
                      geom_segment(aes(x=c(StartX[25]), y=c(HorizontalPosition[25]), xend=c(EndX[25]), yend=c(HorizontalPosition[25])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[25], size=1.5) +
                      geom_segment(aes(x=c(StartX[26]), y=c(HorizontalPosition[26]), xend=c(EndX[26]), yend=c(HorizontalPosition[26])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[26], size=1.5) +
                      geom_segment(aes(x=c(StartX[27]), y=c(HorizontalPosition[27]), xend=c(EndX[27]), yend=c(HorizontalPosition[27])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[27], size=1.5) +
                      geom_segment(aes(x=c(StartX[28]), y=c(HorizontalPosition[28]), xend=c(EndX[28]), yend=c(HorizontalPosition[28])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[28], size=1.5) +
                      geom_segment(aes(x=c(StartX[29]), y=c(HorizontalPosition[29]), xend=c(EndX[29]), yend=c(HorizontalPosition[29])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[29], size=1.5) +
                      geom_segment(aes(x=c(StartX[30]), y=c(HorizontalPosition[30]), xend=c(EndX[30]), yend=c(HorizontalPosition[30])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[30], size=1.5) +

                      geom_segment(aes(x=c(StartX[31]), y=c(HorizontalPosition[31]), xend=c(EndX[31]), yend=c(HorizontalPosition[31])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[31], size=1.5) +
                      geom_segment(aes(x=c(StartX[32]), y=c(HorizontalPosition[32]), xend=c(EndX[32]), yend=c(HorizontalPosition[32])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[32], size=1.5) +
                      geom_segment(aes(x=c(StartX[33]), y=c(HorizontalPosition[33]), xend=c(EndX[33]), yend=c(HorizontalPosition[33])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[33], size=1.5) +
                      geom_segment(aes(x=c(StartX[34]), y=c(HorizontalPosition[34]), xend=c(EndX[34]), yend=c(HorizontalPosition[34])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[34], size=1.5) +
                      geom_segment(aes(x=c(StartX[35]), y=c(HorizontalPosition[35]), xend=c(EndX[35]), yend=c(HorizontalPosition[35])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[35], size=1.5) +
                      geom_segment(aes(x=c(StartX[36]), y=c(HorizontalPosition[36]), xend=c(EndX[36]), yend=c(HorizontalPosition[36])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[36], size=1.5) +
                      geom_segment(aes(x=c(StartX[37]), y=c(HorizontalPosition[37]), xend=c(EndX[37]), yend=c(HorizontalPosition[37])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[37], size=1.5) +
                      geom_segment(aes(x=c(StartX[38]), y=c(HorizontalPosition[38]), xend=c(EndX[38]), yend=c(HorizontalPosition[38])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[38], size=1.5) 
          # } else if (ChrNumb %in% c("Chr1q", "Chr3","Chr4","Chr9","Chr12","Chr18")) {
          #     Myggplot <- Myggplot+
          #           geom_segment(aes(x=c(StartX[1]), y=c(HorizontalPosition[1]), xend=c(EndX[1]), yend=c(HorizontalPosition[1])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[1], size=1.5) +
          #           geom_segment(aes(x=c(StartX[2]), y=c(HorizontalPosition[2]), xend=c(EndX[2]), yend=c(HorizontalPosition[2])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[2], size=1.5) +
          #           geom_segment(aes(x=c(StartX[3]), y=c(HorizontalPosition[3]), xend=c(EndX[3]), yend=c(HorizontalPosition[3])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[3], size=1.5) +
          #           geom_segment(aes(x=c(StartX[4]), y=c(HorizontalPosition[4]), xend=c(EndX[4]), yend=c(HorizontalPosition[4])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[4], size=1.5) +
          #           geom_segment(aes(x=c(StartX[5]), y=c(HorizontalPosition[5]), xend=c(EndX[5]), yend=c(HorizontalPosition[5])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[5], size=1.5) +
          #           geom_segment(aes(x=c(StartX[6]), y=c(HorizontalPosition[6]), xend=c(EndX[6]), yend=c(HorizontalPosition[6])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[6], size=1.5) +
          #           geom_segment(aes(x=c(StartX[7]), y=c(HorizontalPosition[7]), xend=c(EndX[7]), yend=c(HorizontalPosition[7])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[7], size=1.5) +
          #           geom_segment(aes(x=c(StartX[8]), y=c(HorizontalPosition[8]), xend=c(EndX[8]), yend=c(HorizontalPosition[8])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[8], size=1.5) +
          #           geom_segment(aes(x=c(StartX[9]), y=c(HorizontalPosition[9]), xend=c(EndX[9]), yend=c(HorizontalPosition[9])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[9], size=1.5) +
          #           geom_segment(aes(x=c(StartX[10]), y=c(HorizontalPosition[10]), xend=c(EndX[10]), yend=c(HorizontalPosition[10])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[10], size=1.5) +
          #           geom_segment(aes(x=c(StartX[11]), y=c(HorizontalPosition[11]), xend=c(EndX[11]), yend=c(HorizontalPosition[11])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[11], size=1.5)  
          # } ## End of Chr1q or Chr6q

  } else if( WMSubtype=="MBCLike") {
            NumberLine <- length(StartX); print(NumberLine) # AllType: 66 # MBC, 1  # PC, 17         <<<<<<========. This number will be the number of horizontal lines. 
            if(NumberLine==0) NumberLine<-1
            
            if(ChrNumb %in% c("Chr6q","Chr8p","Chr13q","Chr17p","Chr21q" )) {  
                    SpacePerLine <- (0.75-0.15)/NumberLine; print(SpacePerLine) # 0.0078125  <<- make horizontal lines in the y-axis 0.5~0.8 range. # PC: 0.018
                    HorizontalPosition <- seq(from=0.75, to=(0.15+SpacePerLine), by=-SpacePerLine); length(HorizontalPosition) # AllType, 68  # MBC, 32 # PC: 17
        
                    ### Give different color to Del regions by different samples
                    # set.seed=2
                    # DelRegionLineColor <- rand_color(length(unique(OriginalX_DelSize_Sort$sample))); print(DelRegionLineColor) #  "#F674A7FF" "#9BFEB0FF" "#D4CF98FF" "#21FD76FF"
        
        
                    # OriginalX_DelSize_Sort <- x_SmpName_Sort
                    # print(paste0("Number of patients: ", length(unique(OriginalX_DelSize_Sort$sample)))) # 4
                    print(paste0("Number of patients: ", length(unique(OriginalX_DelSize_Sort$sample)))) # 4
                    DelRegionLineColor <- c("green","black","purple","darkorange","grey",   "red","darkgoldenrod","darkgreen","deeppink","skyblue",     "bisque3")
                                          # 01-190     WM70     WM33 
            } else if (ChrNumb %in% c("Chr1q","Chr6p","Chr8q","Chr9q","Chr13p","Chr17q","Chr21p","Chr22p",
                                      "Chr1","Chr2","Chr3","Chr4","Chr5",    "Chr6","Chr7","Chr8","Chr9","Chr10",     "Chr11","Chr12","Chr13","Chr14","Chr15",
                                      "Chr16","Chr17","Chr18","Chr19","Chr20",    "Chr21","Chr22")) {
                    # SpacePerLine <- (-0.9+0.20)/NumberLine; print(SpacePerLine) # 0.0078125  <<- make horizontal lines in the y-axis 0.5~0.8 range. # PC: 0.018
                    # HorizontalPosition <- seq(from=-0.9, to=(-0.20-SpacePerLine), by=-SpacePerLine); length(HorizontalPosition) # AllType, 68  # MBC, 32 # PC: 17
                    
                    SpacePerLine <- (0.95-0.20)/NumberLine; print(SpacePerLine) # 0.0078125  <<- make horizontal lines in the y-axis 0.5~0.8 range. # PC: 0.018
                    HorizontalPosition <- seq(from=0.95, to=(-0.20+SpacePerLine), by=-SpacePerLine); length(HorizontalPosition) # AllType, 68  # MBC, 32 # PC: 17
                    
                    print(paste0("Number of patients: ", length(unique(OriginalX_DelSize_Sort$sample)))) # 6
                    DelRegionLineColor <- c("green","black","purple","darkorange","grey",   "red","darkgoldenrod","darkgreen","deeppink","skyblue",     "bisque3")   # , "red","magenta1", "red4", "goldenrod3")
                                         #     WM33
            }

            DelRegionLineColor_CutByNumbSmp <- DelRegionLineColor[1:length(unique(OriginalX_DelSize_Sort$sample))]   ## I don't need all 11 colors defined in line 416 
            names(DelRegionLineColor_CutByNumbSmp) <- unique(OriginalX_DelSize_Sort$sample); print(DelRegionLineColor_CutByNumbSmp)
            if(is.na(names(DelRegionLineColor_CutByNumbSmp)[1])) {
                    names(DelRegionLineColor_CutByNumbSmp) <- "NoSample"
            }
            DelRegionLineColor_DF <- data.frame(DelRegionLineColor_CutByNumbSmp) %>% tibble::rownames_to_column("sample")
            fwrite(DelRegionLineColor_DF, file=paste0("ColorCodeBySmp_", ChrNumb,"_",WMSubtype,"_v1.7.txt"), col.names=TRUE,row.names=FALSE, sep="\t", quote=FALSE)
            
            ### inner_join between OriginalX_DelSize_Sort and DelRegionLineColor_DF to map the line color per sample ID
            # OriginalX_DelSize_LineColor <- dplyr::inner_join(OriginalX_DelSize_Sort, DelRegionLineColor_DF); dim(OriginalX_DelSize_LineColor)
            OriginalX_DelSize_LineColor <- dplyr::left_join(OriginalX_DelSize_Sort, DelRegionLineColor_DF); dim(OriginalX_DelSize_LineColor)  ## MBC-like Chr6q: 18 19  Chr1q: 15 19

            ### Add Del Region lines with different color by sample IDs. <<<<<< ===== Not working
            # p1 <- p1 + geom_segment(aes(x=c(7*10^7,7*10^7), y=c(0.5, 0.7), xend=c(10*10^7, 11*10^7), yend=c(0.5,0.7) ), color=c("green", "red"), size=1)
            # p1 <- p1 + geom_segment(aes(x=c(7*10^7,7*10^7), y=c(0.5, 0.7), xend=c(10*10^7, 11*10^7), yend=c(0.5,0.7) , color=c("green", "red")), size=1)

            # Horizontal line segment data

            ## geom_segment(aes(x=c(7*10^7), y=c(0.5), xend=c(10*10^7), yend=c(0.5)), color=c("#9BFEB0FF") , size = 1) ### <<<= This works
            ## p1 +geom_segment(aes(x=c(99637055), y=c(0.5), xend=c(109854926), yend=c(0.5)), color=c("#9BFEB0FF") , size = 1) ### <<<= This works
            ## if x and xend are too close, the line is not visible.

            Myggplot <- InitialPlot  ##. Myggplot <- p1  #### If I want to include horizontal line plot into CNV plot. 
                # p1<-p1 + geom_segment(aes(x=c(StartX[1]), y=c(HorizontalPosition[1]), xend=c(EndX[1]+8000000), yend=c(HorizontalPosition[1])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[1], size=1.0)
                # p2<-p2+geom_segment(aes(x=c(StartX[i]), y=c(HorizontalPosition[i]), xend=c(EndX[i]), yend=c(HorizontalPosition[i])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[i], size=6) # this works but only last line.
                # Myggplot<-Myggplot+geom_segment(aes(x=c(99637055, 79967200), y=c(0.5, 0.7), xend=c(10*10^7, 11*10^7), yend=c(0.5, 0.7)), color=rep("green", 2), size=6) # This doesn't work.
           
            # if (ChrNumb=="Chr6q") { 
                     Myggplot<- Myggplot+
                          geom_segment(aes(x=c(StartX[1]), y=c(HorizontalPosition[1]), xend=c(EndX[1]), yend=c(HorizontalPosition[1])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[1], size=1.0) +
                          geom_segment(aes(x=c(StartX[2]), y=c(HorizontalPosition[2]), xend=c(EndX[2]), yend=c(HorizontalPosition[2])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[2], size=1.0) +
                          geom_segment(aes(x=c(StartX[3]), y=c(HorizontalPosition[3]), xend=c(EndX[3]), yend=c(HorizontalPosition[3])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[3], size=1.0) +
                          geom_segment(aes(x=c(StartX[4]), y=c(HorizontalPosition[4]), xend=c(EndX[4]), yend=c(HorizontalPosition[4])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[4], size=1.0) +
                          geom_segment(aes(x=c(StartX[5]), y=c(HorizontalPosition[5]), xend=c(EndX[5]), yend=c(HorizontalPosition[5])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[5], size=1.0) +
                          geom_segment(aes(x=c(StartX[6]), y=c(HorizontalPosition[6]), xend=c(EndX[6]), yend=c(HorizontalPosition[6])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[6], size=1.0) +
                          geom_segment(aes(x=c(StartX[7]), y=c(HorizontalPosition[7]), xend=c(EndX[7]), yend=c(HorizontalPosition[7])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[7], size=1.0) +
                          geom_segment(aes(x=c(StartX[8]), y=c(HorizontalPosition[8]), xend=c(EndX[8]), yend=c(HorizontalPosition[8])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[8], size=1.0) +
                          geom_segment(aes(x=c(StartX[9]), y=c(HorizontalPosition[9]), xend=c(EndX[9]), yend=c(HorizontalPosition[9])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[9], size=1.0) +
                          geom_segment(aes(x=c(StartX[10]), y=c(HorizontalPosition[10]), xend=c(EndX[10]), yend=c(HorizontalPosition[10])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[10], size=1.0) +
          
                          geom_segment(aes(x=c(StartX[11]), y=c(HorizontalPosition[11]), xend=c(EndX[11]), yend=c(HorizontalPosition[11])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[11], size=1.0) +
                          geom_segment(aes(x=c(StartX[12]), y=c(HorizontalPosition[12]), xend=c(EndX[12]), yend=c(HorizontalPosition[12])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[12], size=1.0) +
                          geom_segment(aes(x=c(StartX[13]), y=c(HorizontalPosition[13]), xend=c(EndX[13]), yend=c(HorizontalPosition[13])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[13], size=1.0) +
                          geom_segment(aes(x=c(StartX[14]), y=c(HorizontalPosition[14]), xend=c(EndX[14]), yend=c(HorizontalPosition[14])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[14], size=1.0) +
                          geom_segment(aes(x=c(StartX[15]), y=c(HorizontalPosition[15]), xend=c(EndX[15]), yend=c(HorizontalPosition[15])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[15], size=1.0) +
                          geom_segment(aes(x=c(StartX[16]), y=c(HorizontalPosition[16]), xend=c(EndX[16]), yend=c(HorizontalPosition[16])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[16], size=1.0) +
                          geom_segment(aes(x=c(StartX[17]), y=c(HorizontalPosition[17]), xend=c(EndX[17]), yend=c(HorizontalPosition[17])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[17], size=1.0) +
                          geom_segment(aes(x=c(StartX[18]), y=c(HorizontalPosition[18]), xend=c(EndX[18]), yend=c(HorizontalPosition[18])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[18], size=1.0) 
            # } else if (ChrNumb=="Chr1q") { 
                   # Myggplot<- Myggplot+
                   #       geom_segment(aes(x=c(StartX[1]), y=c(HorizontalPosition[1]), xend=c(EndX[1]), yend=c(HorizontalPosition[1])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[1], size=1.0) +
                   #       geom_segment(aes(x=c(StartX[2]), y=c(HorizontalPosition[2]), xend=c(EndX[2]), yend=c(HorizontalPosition[2])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[2], size=1.0) +
                   #       geom_segment(aes(x=c(StartX[3]), y=c(HorizontalPosition[3]), xend=c(EndX[3]), yend=c(HorizontalPosition[3])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[3], size=1.0) +
                   #       geom_segment(aes(x=c(StartX[4]), y=c(HorizontalPosition[4]), xend=c(EndX[4]), yend=c(HorizontalPosition[4])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[4], size=1.0) +
                   #       geom_segment(aes(x=c(StartX[5]), y=c(HorizontalPosition[5]), xend=c(EndX[5]), yend=c(HorizontalPosition[5])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[5], size=1.0) +
                   #       geom_segment(aes(x=c(StartX[6]), y=c(HorizontalPosition[6]), xend=c(EndX[6]), yend=c(HorizontalPosition[6])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[6], size=1.0) +
                   #       geom_segment(aes(x=c(StartX[7]), y=c(HorizontalPosition[7]), xend=c(EndX[7]), yend=c(HorizontalPosition[7])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[7], size=1.0) +
                   #       geom_segment(aes(x=c(StartX[8]), y=c(HorizontalPosition[8]), xend=c(EndX[8]), yend=c(HorizontalPosition[8])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[8], size=1.0) +
                   #       geom_segment(aes(x=c(StartX[9]), y=c(HorizontalPosition[9]), xend=c(EndX[9]), yend=c(HorizontalPosition[9])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[9], size=1.0) +
                   #       geom_segment(aes(x=c(StartX[10]), y=c(HorizontalPosition[10]), xend=c(EndX[10]), yend=c(HorizontalPosition[10])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[10], size=1.0) +
                   #       
                   #       geom_segment(aes(x=c(StartX[11]), y=c(HorizontalPosition[11]), xend=c(EndX[11]), yend=c(HorizontalPosition[11])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[11], size=1.0) +
                   #       geom_segment(aes(x=c(StartX[12]), y=c(HorizontalPosition[12]), xend=c(EndX[12]), yend=c(HorizontalPosition[12])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[12], size=1.0) +
                   #       geom_segment(aes(x=c(StartX[13]), y=c(HorizontalPosition[13]), xend=c(EndX[13]), yend=c(HorizontalPosition[13])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[13], size=1.0) +
                   #       geom_segment(aes(x=c(StartX[14]), y=c(HorizontalPosition[14]), xend=c(EndX[14]), yend=c(HorizontalPosition[14])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[14], size=1.0) +
                   #       geom_segment(aes(x=c(StartX[15]), y=c(HorizontalPosition[15]), xend=c(EndX[15]), yend=c(HorizontalPosition[15])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[15], size=1.0)
           # } 
  } else if (WMSubtype=="PCLike") {
          NumberLine <- length(StartX); print(NumberLine) # AllType: 66 # MBC, 25  # PC,  13  <<<<<<<=============== Decide the number of horizontal lines. 
          if(NumberLine==0) NumberLine<-1;
          if (ChrNumb %in% c("Chr6q","Chr8p","Chr13q","Chr17p","Chr21q" ) ) {           
                  SpacePerLine <- (0.75-0.2)/NumberLine; print(SpacePerLine) # 0.0078125  <<- make horizontal lines in the y-axis 0.5~0.8 range. # PC: 0.018
                  HorizontalPosition <- seq(from=0.75, to=(0.2+SpacePerLine), by=-SpacePerLine); length(HorizontalPosition) # AllType, 68  # MBC, 25 # PC: 17
        
                  # OriginalX_DelSize_Sort <- x_SmpName_Sort
                  print(paste0("Number of patients: ", length(unique(OriginalX_DelSize_Sort$sample)))) # 22
                  DelRegionLineColor <- c("green","black","purple","darkorange","grey",   "red","darkgoldenrod","darkgreen","deeppink","skyblue",     "bisque3", "aquamarine1","pink","darkkhaki",'yellow',
                                      "coral","cyan" ,"darkgoldenrod1","azure","darkolivegreen1",     "darkblue","deeppink4","darkorchid4","brown1", "cyan1" )
                                    #   WM42            WM43            WM22            WM30            WM69            WM66            WM48            WM68            WM73             WM8            WM32            WM37 
                                    # "green"         "black"        "purple"    "darkorange"          "grey"           "red" "darkgoldenrod"     "darkgreen"      "deeppink"       "skyblue"       "bisque3"   "aquamarine1" 
          } else if (ChrNumb %in% c("Chr1q","Chr6p","Chr8q","Chr9q","Chr13p","Chr17q","Chr21p","Chr22p",
                                    "Chr1","Chr2","Chr3","Chr4","Chr5",    "Chr6","Chr7","Chr8","Chr9","Chr10",     "Chr11","Chr12","Chr13","Chr14","Chr15",
                                    "Chr16","Chr17","Chr18","Chr19","Chr20",    "Chr21","Chr22")) { 
                  # SpacePerLine <- (-0.9+0.3)/NumberLine; print(SpacePerLine) # 0.0078125  <<- make horizontal lines in the y-axis 0.5~0.8 range. # PC: 0.018
                  # HorizontalPosition <- seq(from=-0.9, to=(-0.3-SpacePerLine), by=-SpacePerLine); length(HorizontalPosition) # AllType, 68  # MBC, 25 # PC: 17
                  
                  SpacePerLine <- (0.95-0.4)/NumberLine; print(SpacePerLine) # 0.0078125  <<- make horizontal lines in the y-axis 0.5~0.8 range. # PC: 0.018
                  HorizontalPosition <- seq(from=0.95, to=(-0.4+SpacePerLine), by=-SpacePerLine); length(HorizontalPosition) # AllType, 68  # MBC, 25 # PC: 17
                  
                  # OriginalX_DelSize_Sort <- x_SmpName_Sort
                  # print(paste0("Number of patients: ", length(unique(OriginalX_DelSize_Sort$sample)))) # 12
                  print(paste0("Number of patients: ", length(unique(OriginalX_DelSize_Sort$sample)))) # 12
                  DelRegionLineColor <- c("green","black","purple","darkorange","grey",   "red","darkgoldenrod","darkgreen","deeppink","skyblue",     "bisque3","aquamarine1","pink","darkkhaki",'yellow',          
                                          "coral","cyan" ,"darkgoldenrod1","azure","darkolivegreen1",     "darkblue","deeppink4","darkorchid4","brown1", "cyan1")
                                          # WM49,   WM32,   WM72,   WM43,         WM26,     WM27, WM8 
          }
          
          names(DelRegionLineColor) <- unique(OriginalX_DelSize_Sort$sample); print(DelRegionLineColor)
            # WM42            WM22            WM29            WM30            WM69            WM43          01-131          01-115            WM66            WM48            WM68            WM73
          # "green"         "black"        "purple"    "darkorange"        "grey"           "red" "darkgoldenrod"     "darkgreen"      "deeppink"       "skyblue"       "bisque3"   "aquamarine1"

          DelRegionLineColor_CutByNumbSmp <- DelRegionLineColor[1:length(unique(OriginalX_DelSize_Sort$sample))]   ## I don't need all 11 colors defined in line 416 
          names(DelRegionLineColor_CutByNumbSmp) <- unique(OriginalX_DelSize_Sort$sample); print(DelRegionLineColor)
          if(is.na(names(DelRegionLineColor_CutByNumbSmp)[1])) {
                  names(DelRegionLineColor_CutByNumbSmp) <- "NoSample"
          }
          DelRegionLineColor_DF <- data.frame(DelRegionLineColor_CutByNumbSmp) %>% tibble::rownames_to_column("sample")
          fwrite(DelRegionLineColor_DF, file=paste0("ColorCodeBySmp_", ChrNumb,"_",WMSubtype,"_v1.7.txt"), col.names=TRUE,row.names=FALSE, sep="\t", quote=FALSE)
          
          ### inner_join between OriginalX_DelSize_Sort and DelRegionLineColor_DF to map the line color per sample ID
          # OriginalX_DelSize_LineColor <- dplyr::inner_join(OriginalX_DelSize_Sort, DelRegionLineColor_DF); dim(OriginalX_DelSize_LineColor)
          OriginalX_DelSize_LineColor <- dplyr::inner_join(OriginalX_DelSize_Sort, DelRegionLineColor_DF); dim(OriginalX_DelSize_LineColor) ## PC-like Chr1q 7 19  ## Chr6q, 18 19

          ## Chr1q: WM49 green, WM32 black, WM72 purple, WM43 darkorange, WM26 grey, WM27 red, WM8 darkgoldenrod
          
          ### Add Del Region lines with different color by sample IDs. <<<<<< ===== Not working
          # p1 <- p1 + geom_segment(aes(x=c(7*10^7,7*10^7), y=c(0.5, 0.7), xend=c(10*10^7, 11*10^7), yend=c(0.5,0.7) ), color=c("green", "red"), size=1)

          # Horizontal line segment data
          ## geom_segment(aes(x=c(7*10^7), y=c(0.5), xend=c(10*10^7), yend=c(0.5)), color=c("#9BFEB0FF") , size = 1) ### <<<= This works
          ## p1 +geom_segment(aes(x=c(99637055), y=c(0.5), xend=c(109854926), yend=c(0.5)), color=c("#9BFEB0FF") , size = 1) ### <<<= This works
          ## if x and xend are too close, the line is not visible.

          Myggplot <- InitialPlot  ##. Myggplot <- p1  #### If I want to include horizontal line plot into CNV plot. 
          # if (ChrNumb=="Chr6q") {    
                  # p1<-p1 + geom_segment(aes(x=c(StartX[1]), y=c(HorizontalPosition[1]), xend=c(EndX[1]+8000000), yend=c(HorizontalPosition[1])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[1], size=1.0)
                  # p2<-p2+geom_segment(aes(x=c(StartX[i]), y=c(HorizontalPosition[i]), xend=c(EndX[i]), yend=c(HorizontalPosition[i])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[i], size=6) # this works but only last line.
                  # Myggplot<-Myggplot+geom_segment(aes(x=c(99637055, 79967200), y=c(0.5, 0.7), xend=c(10*10^7, 11*10^7), yend=c(0.5, 0.7)), color=rep("green", 2), size=6) # This doesn't work.
                  Myggplot<- Myggplot+
                    geom_segment(aes(x=c(StartX[1]), y=c(HorizontalPosition[1]), xend=c(EndX[1]), yend=c(HorizontalPosition[1])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[1], size=1.0) +
                    geom_segment(aes(x=c(StartX[2]), y=c(HorizontalPosition[2]), xend=c(EndX[2]), yend=c(HorizontalPosition[2])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[2], size=1.0) +
                    geom_segment(aes(x=c(StartX[3]), y=c(HorizontalPosition[3]), xend=c(EndX[3]), yend=c(HorizontalPosition[3])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[3], size=1.0) +
                    geom_segment(aes(x=c(StartX[4]), y=c(HorizontalPosition[4]), xend=c(EndX[4]), yend=c(HorizontalPosition[4])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[4], size=1.0) +
                    geom_segment(aes(x=c(StartX[5]), y=c(HorizontalPosition[5]), xend=c(EndX[5]), yend=c(HorizontalPosition[5])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[5], size=1.0) +
                    geom_segment(aes(x=c(StartX[6]), y=c(HorizontalPosition[6]), xend=c(EndX[6]), yend=c(HorizontalPosition[6])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[6], size=1.0) +
                    geom_segment(aes(x=c(StartX[7]), y=c(HorizontalPosition[7]), xend=c(EndX[7]), yend=c(HorizontalPosition[7])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[7], size=1.0) +
                    geom_segment(aes(x=c(StartX[8]), y=c(HorizontalPosition[8]), xend=c(EndX[8]), yend=c(HorizontalPosition[8])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[8], size=1.0) +
                    geom_segment(aes(x=c(StartX[9]), y=c(HorizontalPosition[9]), xend=c(EndX[9]), yend=c(HorizontalPosition[9])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[9], size=1.0) +
                    geom_segment(aes(x=c(StartX[10]), y=c(HorizontalPosition[10]), xend=c(EndX[10]), yend=c(HorizontalPosition[10])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[10], size=1.0) +
        
                    geom_segment(aes(x=c(StartX[11]), y=c(HorizontalPosition[11]), xend=c(EndX[11]), yend=c(HorizontalPosition[11])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[11], size=1.0) +
                    geom_segment(aes(x=c(StartX[12]), y=c(HorizontalPosition[12]), xend=c(EndX[12]), yend=c(HorizontalPosition[12])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[12], size=1.0) +
                    geom_segment(aes(x=c(StartX[13]), y=c(HorizontalPosition[13]), xend=c(EndX[13]), yend=c(HorizontalPosition[13])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[13], size=1.0) +
                    geom_segment(aes(x=c(StartX[14]), y=c(HorizontalPosition[14]), xend=c(EndX[14]), yend=c(HorizontalPosition[14])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[14], size=1.0) +
                    geom_segment(aes(x=c(StartX[15]), y=c(HorizontalPosition[15]), xend=c(EndX[15]), yend=c(HorizontalPosition[15])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[15], size=1.0) +
                    geom_segment(aes(x=c(StartX[16]), y=c(HorizontalPosition[16]), xend=c(EndX[16]), yend=c(HorizontalPosition[16])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[16], size=1.0) +
                    geom_segment(aes(x=c(StartX[17]), y=c(HorizontalPosition[17]), xend=c(EndX[17]), yend=c(HorizontalPosition[17])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[17], size=1.0) +
                    geom_segment(aes(x=c(StartX[18]), y=c(HorizontalPosition[18]), xend=c(EndX[18]), yend=c(HorizontalPosition[18])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[18], size=1.0) +
                    geom_segment(aes(x=c(StartX[19]), y=c(HorizontalPosition[19]), xend=c(EndX[19]), yend=c(HorizontalPosition[19])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[19], size=1.0) +
                    geom_segment(aes(x=c(StartX[20]), y=c(HorizontalPosition[20]), xend=c(EndX[20]), yend=c(HorizontalPosition[20])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[20], size=1.0) +
                    
                    geom_segment(aes(x=c(StartX[21]), y=c(HorizontalPosition[21]), xend=c(EndX[21]), yend=c(HorizontalPosition[21])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[21], size=1.0) +
                    geom_segment(aes(x=c(StartX[22]), y=c(HorizontalPosition[22]), xend=c(EndX[22]), yend=c(HorizontalPosition[22])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[22], size=1.0) 
            # } else if (ChrNumb=="Chr1q") {    
            #       Myggplot<- Myggplot+
            #         geom_segment(aes(x=c(StartX[1]), y=c(HorizontalPosition[1]), xend=c(EndX[1]), yend=c(HorizontalPosition[1])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[1], size=1.0) +
            #         geom_segment(aes(x=c(StartX[2]), y=c(HorizontalPosition[2]), xend=c(EndX[2]), yend=c(HorizontalPosition[2])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[2], size=1.0) +
            #         geom_segment(aes(x=c(StartX[3]), y=c(HorizontalPosition[3]), xend=c(EndX[3]), yend=c(HorizontalPosition[3])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[3], size=1.0) +
            #         geom_segment(aes(x=c(StartX[4]), y=c(HorizontalPosition[4]), xend=c(EndX[4]), yend=c(HorizontalPosition[4])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[4], size=1.0) +
            #         geom_segment(aes(x=c(StartX[5]), y=c(HorizontalPosition[5]), xend=c(EndX[5]), yend=c(HorizontalPosition[5])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[5], size=1.0) +
            #         geom_segment(aes(x=c(StartX[6]), y=c(HorizontalPosition[6]), xend=c(EndX[6]), yend=c(HorizontalPosition[6])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[6], size=1.0) +
            #         geom_segment(aes(x=c(StartX[7]), y=c(HorizontalPosition[7]), xend=c(EndX[7]), yend=c(HorizontalPosition[7])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[7], size=1.0) 
            #         
            #         geom_segment(aes(x=c(StartX[11]), y=c(HorizontalPosition[11]), xend=c(EndX[11]), yend=c(HorizontalPosition[11])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[11], size=1.0) 
            # }
  } else if (WMSubtype=="NotAvail") {
          NumberLine <- length(StartX); print(NumberLine) # Chr3:  # AllType: 66 # MBC, 25  # PC, 17 # NotAvail: 4
          if(NumberLine==0) NumberLine<-1;
          if (ChrNumb %in% c("Chr6q","Chr8p","Chr13q","Chr17p","Chr21q" ) ) {    
                SpacePerLine <- (0.75-0.2)/NumberLine; print(SpacePerLine) # 0.0078125  <<- make horizontal lines in the y-axis 0.5~0.8 range. # PC: 0.018
                HorizontalPosition <- seq(from=0.75, to=(0.2+SpacePerLine), by=-SpacePerLine); length(HorizontalPosition) # AllType, 68  # MBC, 25 # PC: 17
      
                # OriginalX_DelSize_Sort <- x_SmpName_Sort
                # print(paste0("Number of patients: ", length(unique(OriginalX_DelSize_Sort$sample))))
                print(paste0("Number of patients: ", length(unique(OriginalX_DelSize_Sort$sample)))) # 8
                DelRegionLineColor <- c("green","black","purple","darkorange","grey")     #  ,"red","darkgoldenrod","darkgreen")
                                    #   # WM82   WM96      WM86      WM84     WM89 
          } else if (ChrNumb %in% c("Chr1q","Chr6p","Chr8q","Chr9q","Chr13p","Chr17q","Chr21p","Chr22p",
                                    "Chr1","Chr2","Chr3","Chr4","Chr5",    "Chr6","Chr7","Chr8","Chr9","Chr10",     "Chr11","Chr12","Chr13","Chr14","Chr15",
                                    "Chr16","Chr17","Chr18","Chr19","Chr20",    "Chr21","Chr22")) {
                # SpacePerLine <- (-0.95+0.3)/NumberLine; print(SpacePerLine) # 0.0078125  <<- make horizontal lines in the y-axis 0.5~0.8 range. # PC: 0.018
                # HorizontalPosition <- seq(from=-0.9, to=(-0.3-SpacePerLine), by=-SpacePerLine); length(HorizontalPosition) # AllType, 68  # MBC, 25 # PC: 17
                
                SpacePerLine <- (0.95-0.4)/NumberLine; print(SpacePerLine) # 0.0078125  <<- make horizontal lines in the y-axis 0.5~0.8 range. # PC: 0.018
                HorizontalPosition <- seq(from=0.95, to=(-0.4+SpacePerLine), by=-SpacePerLine); length(HorizontalPosition) # AllType, 68  # MBC, 25 # PC: 17
                
                # OriginalX_DelSize_Sort <- x_SmpName_Sort
                # print(paste0("Number of patients: ", length(unique(OriginalX_DelSize_Sort$sample))))
                print(paste0("Number of patients: ", length(unique(OriginalX_DelSize_Sort$sample)))) # 8
                DelRegionLineColor <- c("green","black","purple","darkorange","grey",   "red","darkgoldenrod","darkgreen","deeppink","skyblue",     "bisque3","aquamarine1","pink","darkkhaki",'yellow',          
                                        "coral","cyan" ,"darkgoldenrod1","azure","darkolivegreen1",     "darkblue","deeppink4","darkorchid4","brown1", "cyan1")
                
          }
                
          names(DelRegionLineColor) <- unique(OriginalX_DelSize_Sort$sample); print(DelRegionLineColor)
          # WM82         WM96         WM86         WM84         WM89 
          # "green"      "black"     "purple" "darkorange"       "grey" 

          DelRegionLineColor_CutByNumbSmp <- DelRegionLineColor[1:length(unique(OriginalX_DelSize_Sort$sample))]   ## I don't need all 11 colors defined in line 416 
          names(DelRegionLineColor_CutByNumbSmp) <- unique(OriginalX_DelSize_Sort$sample); print(DelRegionLineColor)
          if(is.na(names(DelRegionLineColor_CutByNumbSmp)[1])) {
                  names(DelRegionLineColor_CutByNumbSmp) <- "NoSample"
          }
          DelRegionLineColor_DF <- data.frame(DelRegionLineColor_CutByNumbSmp) %>% tibble::rownames_to_column("sample")
          fwrite(DelRegionLineColor_DF, file=paste0("ColorCodeBySmp_", ChrNumb,"_",WMSubtype,"_v1.7.txt"), col.names=TRUE,row.names=FALSE, sep="\t", quote=FALSE)
          
          ### inner_join between OriginalX_DelSize_Sort and DelRegionLineColor_DF to map the line color per sample ID
          # OriginalX_DelSize_LineColor <- dplyr::inner_join(OriginalX_DelSize_Sort, DelRegionLineColor_DF); dim(OriginalX_DelSize_LineColor)
          OriginalX_DelSize_LineColor <- dplyr::inner_join(OriginalX_DelSize_Sort, DelRegionLineColor_DF); dim(OriginalX_DelSize_LineColor) # Chr1q: 22 19  #  Chr6q 11 19  <<<<<===== 
          #### This is the number of horizontal lines. 

          ### Add Del Region lines with different color by sample IDs. <<<<<< ===== Not working
          # p1 <- p1 + geom_segment(aes(x=c(7*10^7,7*10^7), y=c(0.5, 0.7), xend=c(10*10^7, 11*10^7), yend=c(0.5,0.7) ), color=c("green", "red"), size=1)

          # Horizontal line segment data
          ## geom_segment(aes(x=c(7*10^7), y=c(0.5), xend=c(10*10^7), yend=c(0.5)), color=c("#9BFEB0FF") , size = 1) ### <<<= This works
          ## p1 +geom_segment(aes(x=c(99637055), y=c(0.5), xend=c(109854926), yend=c(0.5)), color=c("#9BFEB0FF") , size = 1) ### <<<= This works
          ## if x and xend are too close, the line is not visible.

          Myggplot <- InitialPlot  ##. Myggplot <- p1  #### If I want to include horizontal line plot into CNV plot. 
          
          # if (ChrNumb=="Chr6q") {    
                # p1<-p1 + geom_segment(aes(x=c(StartX[1]), y=c(HorizontalPosition[1]), xend=c(EndX[1]+8000000), yend=c(HorizontalPosition[1])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[1], size=1.0)
                # p2<-p2+geom_segment(aes(x=c(StartX[i]), y=c(HorizontalPosition[i]), xend=c(EndX[i]), yend=c(HorizontalPosition[i])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[i], size=6) # this works but only last line.
                # Myggplot<-Myggplot+geom_segment(aes(x=c(99637055, 79967200), y=c(0.5, 0.7), xend=c(10*10^7, 11*10^7), yend=c(0.5, 0.7)), color=rep("green", 2), size=6) # This doesn't work.
                Myggplot<- Myggplot+
                  geom_segment(aes(x=c(StartX[1]), y=c(HorizontalPosition[1]), xend=c(EndX[1]), yend=c(HorizontalPosition[1])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[1], size=1.0) +
                  geom_segment(aes(x=c(StartX[2]), y=c(HorizontalPosition[2]), xend=c(EndX[2]), yend=c(HorizontalPosition[2])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[2], size=1.0) +
                  geom_segment(aes(x=c(StartX[3]), y=c(HorizontalPosition[3]), xend=c(EndX[3]), yend=c(HorizontalPosition[3])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[3], size=1.0) +
                  geom_segment(aes(x=c(StartX[4]), y=c(HorizontalPosition[4]), xend=c(EndX[4]), yend=c(HorizontalPosition[4])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[4], size=1.0) +
                  geom_segment(aes(x=c(StartX[5]), y=c(HorizontalPosition[5]), xend=c(EndX[5]), yend=c(HorizontalPosition[5])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[5], size=1.0) +
                  geom_segment(aes(x=c(StartX[6]), y=c(HorizontalPosition[6]), xend=c(EndX[6]), yend=c(HorizontalPosition[6])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[6], size=1.0) +
                  geom_segment(aes(x=c(StartX[7]), y=c(HorizontalPosition[7]), xend=c(EndX[7]), yend=c(HorizontalPosition[7])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[7], size=1.0) +
                  geom_segment(aes(x=c(StartX[8]), y=c(HorizontalPosition[8]), xend=c(EndX[8]), yend=c(HorizontalPosition[8])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[8], size=1.0) +
                  geom_segment(aes(x=c(StartX[9]), y=c(HorizontalPosition[9]), xend=c(EndX[9]), yend=c(HorizontalPosition[9])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[9], size=1.0) +
                  geom_segment(aes(x=c(StartX[10]), y=c(HorizontalPosition[10]), xend=c(EndX[10]), yend=c(HorizontalPosition[10])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[10], size=1.0) +
                  
                  geom_segment(aes(x=c(StartX[11]), y=c(HorizontalPosition[11]), xend=c(EndX[11]), yend=c(HorizontalPosition[11])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[11], size=1.0) +
                  geom_segment(aes(x=c(StartX[12]), y=c(HorizontalPosition[12]), xend=c(EndX[12]), yend=c(HorizontalPosition[12])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[12], size=1.0) +
                  geom_segment(aes(x=c(StartX[13]), y=c(HorizontalPosition[13]), xend=c(EndX[13]), yend=c(HorizontalPosition[13])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[13], size=1.0) +
                  geom_segment(aes(x=c(StartX[14]), y=c(HorizontalPosition[14]), xend=c(EndX[14]), yend=c(HorizontalPosition[14])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[14], size=1.0) +
                  geom_segment(aes(x=c(StartX[15]), y=c(HorizontalPosition[15]), xend=c(EndX[15]), yend=c(HorizontalPosition[15])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[15], size=1.0) +
                  geom_segment(aes(x=c(StartX[16]), y=c(HorizontalPosition[16]), xend=c(EndX[16]), yend=c(HorizontalPosition[16])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[16], size=1.0) 
          # } else if (ChrNumb=="Chr1q") {
          #       Myggplot<- Myggplot+
          #         geom_segment(aes(x=c(StartX[1]), y=c(HorizontalPosition[1]), xend=c(EndX[1]), yend=c(HorizontalPosition[1])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[1], size=1.0) +
          #         geom_segment(aes(x=c(StartX[2]), y=c(HorizontalPosition[2]), xend=c(EndX[2]), yend=c(HorizontalPosition[2])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[2], size=1.0) +
          #         geom_segment(aes(x=c(StartX[3]), y=c(HorizontalPosition[3]), xend=c(EndX[3]), yend=c(HorizontalPosition[3])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[3], size=1.0) +
          #         geom_segment(aes(x=c(StartX[4]), y=c(HorizontalPosition[4]), xend=c(EndX[4]), yend=c(HorizontalPosition[4])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[4], size=1.0) +
          #         geom_segment(aes(x=c(StartX[5]), y=c(HorizontalPosition[5]), xend=c(EndX[5]), yend=c(HorizontalPosition[5])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[5], size=1.0) +
          #         geom_segment(aes(x=c(StartX[6]), y=c(HorizontalPosition[6]), xend=c(EndX[6]), yend=c(HorizontalPosition[6])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[6], size=1.0) +
          #         geom_segment(aes(x=c(StartX[7]), y=c(HorizontalPosition[7]), xend=c(EndX[7]), yend=c(HorizontalPosition[7])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[7], size=1.0) +
          #         geom_segment(aes(x=c(StartX[8]), y=c(HorizontalPosition[8]), xend=c(EndX[8]), yend=c(HorizontalPosition[8])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[8], size=1.0) +
          #         geom_segment(aes(x=c(StartX[9]), y=c(HorizontalPosition[9]), xend=c(EndX[9]), yend=c(HorizontalPosition[9])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[9], size=1.0) +
          #         geom_segment(aes(x=c(StartX[10]), y=c(HorizontalPosition[10]), xend=c(EndX[10]), yend=c(HorizontalPosition[10])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[10], size=1.0) +
          #         
          #         geom_segment(aes(x=c(StartX[11]), y=c(HorizontalPosition[11]), xend=c(EndX[11]), yend=c(HorizontalPosition[11])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[11], size=1.0) +
          #         geom_segment(aes(x=c(StartX[12]), y=c(HorizontalPosition[12]), xend=c(EndX[12]), yend=c(HorizontalPosition[12])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[12], size=1.0) +
          #         geom_segment(aes(x=c(StartX[13]), y=c(HorizontalPosition[13]), xend=c(EndX[13]), yend=c(HorizontalPosition[13])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[13], size=1.0) +
          #         geom_segment(aes(x=c(StartX[14]), y=c(HorizontalPosition[14]), xend=c(EndX[14]), yend=c(HorizontalPosition[14])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[14], size=1.0) +
          #         geom_segment(aes(x=c(StartX[15]), y=c(HorizontalPosition[15]), xend=c(EndX[15]), yend=c(HorizontalPosition[15])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[15], size=1.0) +
          #         geom_segment(aes(x=c(StartX[16]), y=c(HorizontalPosition[16]), xend=c(EndX[16]), yend=c(HorizontalPosition[16])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[16], size=1.0) +
          #         geom_segment(aes(x=c(StartX[17]), y=c(HorizontalPosition[17]), xend=c(EndX[17]), yend=c(HorizontalPosition[17])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[17], size=1.0) +
          #         geom_segment(aes(x=c(StartX[18]), y=c(HorizontalPosition[18]), xend=c(EndX[18]), yend=c(HorizontalPosition[18])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[18], size=1.0) +
          #         geom_segment(aes(x=c(StartX[19]), y=c(HorizontalPosition[19]), xend=c(EndX[19]), yend=c(HorizontalPosition[19])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[19], size=1.0) +
          #         geom_segment(aes(x=c(StartX[20]), y=c(HorizontalPosition[20]), xend=c(EndX[20]), yend=c(HorizontalPosition[20])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[20], size=1.0) +
          #         
          #         geom_segment(aes(x=c(StartX[21]), y=c(HorizontalPosition[21]), xend=c(EndX[21]), yend=c(HorizontalPosition[21])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[21], size=1.0) +
          #         geom_segment(aes(x=c(StartX[22]), y=c(HorizontalPosition[22]), xend=c(EndX[22]), yend=c(HorizontalPosition[22])), color=OriginalX_DelSize_LineColor$DelRegionLineColor[22], size=1.0) 
          # }     
  }
  
  if(ChrNumb %in% c("Chr6q","Chr8p","Chr13q","Chr17p","Chr21q" ) ) {  #### <<<< ==== For LOSS
        ggplot2::ggsave(filename=paste0("GenVisR_HorizontalLinePlot_", ChrNumb, "Only_",WMSubtype, "_Loss_v1.7.pdf"), Myggplot, width=8, height=4)
  } else if(ChrNumb %in% c("Chr1q","Chr6p","Chr8q","Chr9q","Chr13p","Chr17q","Chr21p","Chr22p") ) { 
        ggplot2::ggsave(filename=paste0("GenVisR_HorizontalLinePlot_", ChrNumb, "Only_",WMSubtype, "_Gain_v1.7.pdf"), Myggplot, width=8, height=4)
  } else if(ChrNumb %in% c("Chr1","Chr2","Chr3","Chr4","Chr5",    "Chr6","Chr7","Chr8","Chr9","Chr10",     "Chr11","Chr12","Chr13","Chr14","Chr15",
                          "Chr16","Chr17","Chr18","Chr19","Chr20",    "Chr21","Chr22")) { 
        ggplot2::ggsave(filename=paste0("GenVisR_HorizontalLinePlot_", ChrNumb, "_",WMSubtype, "_WholeRegionGain_v1.7.pdf"), Myggplot, width=8, height=4)
  } 
############ +++++++++++++++++ ################# ++++++++++++++++++++++ ############## ++++++++

  # if there are other layers, add them
  if(!is.null(plotLayer))   {
    ## p1 <- Myggplot + plotLayer  ####.  <<===== If I want to include Myggplot (horizontal lines) and p1(CNV plot) into one plot. 
    p1 <- p1 + plotLayer ####  <<===== If I want to make only p1(CNV plot) one plot. 
  }

  # if title is supplied plot it
  if(!is.null(plot_title))   {
    p1 <- Myggplot + ggtitle(plot_title)
  }

  return(p1)
}



######### ++++++++++++ ######### ++++++++++++ ######### ++++++++++++ ######### ++++++++++++ ######### ++++++++++++ ######### ++++++++++++ ######### ++++++++++++ 
######### ++++++++++++ ######### ++++++++++++ ######### ++++++++++++ ######### ++++++++++++ ######### ++++++++++++ ######### ++++++++++++ ######### ++++++++++++ 



#############################################################
## original cnFreq()

cnData_NoWM91WM93WM94 <- readRDS("cnData_NoWM91WM93WM94.rds")
x=cnData_NoWM91WM93WM94

cnFreq <- function(x, CN_low_cutoff=1.5, CN_high_cutoff=2.5, plot_title=NULL,
                   CN_Loss_colour='#002EB8', CN_Gain_colour='#A30000',
                   x_title_size=12, y_title_size=12, facet_lab_size=10,
                   plotLayer=NULL, plotType="proportion", genome="hg19",
                   plotChr=NULL, out="plot") {
  # Perform quality check on input data
  x <- cnFreq_qual(x)
  samples <- unique(x$sample)

  # Calculate a columns of Observed CN gains/losses/and obs samples in the
  # cohort for each segment
  gainFreq <- function(x){length(x[x >= CN_high_cutoff])}
  gainFrequency <- aggregate(segmean ~ chromosome + start + end, data=x, gainFreq)$segmean

  lossFreq <- function(x){length(x[x <= CN_low_cutoff])}
  lossFrequency <- aggregate(segmean ~ chromosome + start + end, data=x, lossFreq)$segmean

  x <- aggregate(segmean ~ chromosome + start + end, data=x, length)
  colnames(x)[which(colnames(x) %in% "segmean")] <- "sampleFrequency"
  x$gainFrequency <- gainFrequency
  x$lossFrequency <- lossFrequency

  # check for coordinate space, if any widths are 1 it might indicate a problem
  if(max(x$sampleFrequency) > length(samples)){
    memo <- paste0("Detected additional sample rows after disjoin operation",
                   " typically this indicates coordinates are 0-based, please convert",
                   " coordinates to 1-base for accurate results")
    warning(memo)
  }

  # Calculate the proportion
  x$gainProportion <- x$gainFrequency/length(samples)
  x$lossProportion <- x$lossFrequency/length(samples)

  # get the dummy data for plot boundaries
  preloaded <- c("hg38", "hg19", "mm10", "mm9", "rn5")
  if(any(genome == preloaded)){
    message("genome specified is preloaded, retrieving data...")
    UCSC_Chr_pos <- GenVisR::cytoGeno[GenVisR::cytoGeno$genome == genome,]
    UCSC_Chr_pos <- multi_chrBound(UCSC_Chr_pos)
  } else {
    # Obtain data for UCSC genome and extract relevant columns
    memo <- paste0("attempting to query UCSC mySQL database for chromosome",
                   " positions")
    message(memo)
    cyto_data <- suppressWarnings(multi_cytobandRet(genome))
    UCSC_Chr_pos <- multi_chrBound(cyto_data)
  }

  # check that the dummy data has a size
  if(nrow(UCSC_Chr_pos) < 1)
  {
    memo <- paste0("did not recognize genome ", genome,
                   ", plotting provided data and ignoring chromosome ",
                   "boundaries! Output could be decieving!")
    warning(memo)
  }

  dummy_data <- lapply(unique(x$sample),
                       function(sample, chr_pos) cbind(chr_pos, sample),
                       UCSC_Chr_pos)
  dummy_data <- do.call("rbind", dummy_data)
  chr_order <- gtools::mixedsort(unique(dummy_data$chromosome))
  dummy_data$chromosome <- factor(dummy_data$chromosome, levels=chr_order)

  # select chromosomes to plot
  if(!is.null(plotChr)){
    if(any(!plotChr %in% dummy_data$chromosome)) {
      missingChr <- plotChr[!plotChr %in% dummy_data$chromosome]
      plotChr <- plotChr[!plotChr %in% missingChr]
      memo <- paste0("The following chromosomes: ", toString(missingChr),
                     ", could not be found! Valid chromosomes are: ",
                     toString(unique(dummy_data$chromosome)))
      warning(memo)
    }

    dummy_data <- dummy_data[dummy_data$chromosome %in% plotChr,]
    dummy_data$chromosome <- factor(dummy_data$chromosome, levels=plotChr)
    x <- x[x$chromosome %in% plotChr,]
    x$chromosome <- factor(x$chromosome, levels=plotChr)
  }

  # build the plot
  p1 <- cnFreq_buildMain(x, plotType, dummy_data = dummy_data, plot_title=plot_title,
                         CN_low_colour=CN_Loss_colour,
                         CN_high_colour=CN_Gain_colour,
                         x_lab_size=x_title_size,
                         y_lab_size=y_title_size,
                         facet_lab_size=facet_lab_size,
                         plotLayer=plotLayer)

  # Decide what to output
  output <- multi_selectOut(data=list("data"=x), plot=p1, out=out)
  return(output)
}

