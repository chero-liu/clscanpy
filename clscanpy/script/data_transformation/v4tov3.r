## module load OESingleCell/3.0.d &&  Rscript scripts/R/v4tov3.r -i {input.cluster_rds} -f h5seurat -o result/rds &&
## module unload OESingleCell/3.0.d && module load OESingleCell/2.0.0 && Rscript scripts/R/v4tov3.r -o result/rds

getNumCores= function () 
{
    n.cores <- as.integer(Sys.getenv("LSB_DJOB_NUMPROC"))
    if (is.na(n.cores)) {
        n.cores <- parallel::detectCores()
    }
    if (n.cores > parallel::detectCores() | n.cores < 1) {
        stop("illegal number of cores (", n.cores, ") specified.")
    }
    return(n.cores)
}

saveRDSMC=function (object, file, threads = getNumCores()) 
{
    message("using ", threads, " threads for compression.")
    con <- pipe(paste0("pigz -p ", threads, " -9 -f > ", file),         "wb")
    saveRDS(object, file = con)
    on.exit(if (exists("con")) close(con))
}
readRDSMC= function (file, cores)
{
    tictoc::tic()
    con <- pipe(paste0("cat ", file, " | pigz -dcp ", cores),
        "rb")
    object <- tryCatch({
        readRDS(file = con)
    }, error = function(err) {
        stop("could not read file\n", file, ":\n", err)
    }, finally = {
        if (exists("con"))
            close(con)
    })
    tictoc::toc()
    return(object)
}

suppressPackageStartupMessages( library("argparse") )
suppressPackageStartupMessages( library("Seurat") )

parser = ArgumentParser(description = "single cell sequencing data manipulating toolsets.",
                        usage = "%(prog)s [global options]" )
parser$add_argument("-i", "--input", type = "character",
             help = "The input seurat object in several possible format.")
parser$add_argument("-f", "--informat", type = "character", default = "h5seurat",
             help = "The format of data object, the possible choices can be:h5seurat,(seurat)rds,(sce)rds, loom.[default: %(default)s]")
parser$add_argument("-o", "--output", type = "character", default = "./",
             help = "the output directory of results."  )
opt = parser$parse_args()

if ( is.null(opt$output) ){
    print("NO output directory specified,the current directory will be used!")
    output_dir = getwd()
}else{
    if ( file.exists(opt$output) ){
        output_dir = opt$output
    }else{
        output_dir = opt$output
        dir.create(output_dir,recursive = T)
    }
}
## ==========================================
if (substring(packageVersion("Seurat"),1,1) == 4) {
    print("saving v4 results...")
    output_dir = normalizePath(output_dir )
    data_ob = OESingleCell::ReadX(input = opt$input, informat = opt$informat,verbose = F)
    data_list=list()
    data_list$counts =  data_ob@assays$RNA@counts
    data_list$data =  data_ob@assays$RNA@data
    data_list$scale.data =  data_ob@assays$RNA@scale.data
    data_list$meta.data =  data_ob@meta.data
    data_list$var.features <- data_ob@assays$RNA@var.features
    data_list$meta.features <- data_ob@assays$RNA@meta.features
    #commands
    for (i in names(data_ob@commands)){
    #    class(data_ob@commands[[i]]) = "SeuratCommand"
        attributes(class(data_ob@commands[[i]]))$package = "Seurat"
    } 
    data_list$commands =  data_ob@commands
    data_list$reduction=list()
    for ( i in names(data_ob@reductions)) {
        data_list$reduction[[i]]=data_ob@reductions[[i]]@cell.embeddings
    }
    saveRDSMC(data_list,file.path(output_dir,"data_list.rds"))
    # data_ob@assays # counts data scale.data # var.features # key
    # data_ob@active.assay
    # data_ob@reductions # cell.embeddings assay.used key
    # data_ob@meta.data
    # data_ob@active.ident
} else {
    print("building v3 object...")
    if ( file.exists(file.path(output_dir,"data_list.rds")) ){
        data_list= readRDSMC(file.path(output_dir,"data_list.rds"),20)
    } else {
        stop(file.path(output_dir,"data_list.rds"), "does not exist. Please run in seurat v4 first.")
    }
    data_ob=CreateSeuratObject(counts=data_list$counts,meta.data=data_list$meta.data)
    data_ob@assays$RNA@data = data_list$data
    data_ob@assays$RNA@scale.data = data_list$scale.data
    data_ob@assays$RNA@var.features <- data_list$var.features
    data_ob@assays$RNA@meta.features <- data_list$meta.features
    data_ob@commands = data_list$commands
    # reduction
    reductions = list()
    for ( i in names(data_list$reduction )) {
        reductions[[i]] = CreateDimReducObject(embeddings=data_list$reduction[[i]] ,key= paste0(toupper(i),"_"),global=T,assay="RNA")
    }
    data_ob@reductions = reductions
    # rename
    if(grepl("-",Cells(data_ob)[1])){
        newnames= unlist(lapply( strsplit(Cells(data_ob),split="-"),
            FUN = function (x) {paste(x[1],x[2],sep="_")}))
        data_ob = RenameCells(data_ob,new.names= newnames)
    }
    data_ob=SetIdent(data_ob,value="seurat_clusters")
    # rm intermediate
    flag =  try(saveRDSMC(data_ob,file.path(output_dir,"data_ob_v3.rds")))
    if (class(flag)=="try-error" ) { 
        print("An error occur while saving v3rds.")
    } else {
        print("Removing tmp file...")
        file.remove(file.path(output_dir,"data_list.rds"))
    }    
}