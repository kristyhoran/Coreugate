require(readr)
require(tidyverse)
require(readxl)

# install.packages("Rcpp")
library(Rcpp)

# A C function to calculate pairwise allelic distance quickly!!
cppFunction('NumericMatrix distance_kh(CharacterMatrix df, bool remove_na, bool prop, String na){
            // rows and cols are the lenght of isolates, from L
            //a matrix for distance matrix, intialised to number of rows (isolates)
            int L = df.nrow();
            // cols for the df
            int C = df.ncol();
            
            //Rcout << L << std::endl << C << std::endl;
            NumericMatrix mat(L,L);
            //rows up to L
            for (int i=0; i <L; ++i){
            // rows from i to L so as not to repeat
            for (int j = i; j<L; ++j){
            // d is allelic distance = number of mismatches
            double d = 0;
            int pos = 0;  
            // each loci up to C
            for (int k=0;k<C;++k){
            // Rcout << df(i,k) << std::endl << df(j,k)  <<std::endl;
            
            if (remove_na){
            if (df(i,k) != na && df(j,k) != na){
            pos = pos + 1;
            //Rcout << pos << std::endl;
            if (df(i,k) != df(j,k)){
            d = d +1;
            //Rcout << d << std::endl;
            } // end of if df(i,k) != df(j,k)
            } // end if for na string
            
            } // end if remove na
            else {
            pos = pos + 1;
            if (df(i,k) != df(j,k)){
            d = d +1;
            //Rcout << d << std::endl;
            } // end of if df(i,k) != df(j,k)
            } // end else
            } //end of k loop
            //Rcout << pos << std::endl;
            //Rcout << d << std::endl;
            //Rcout << prop << std::endl;
            if (prop){
            d = d/pos;
            
            }
            
            mat(i,j) =d;
            mat(j,i) = d;
            
            
            } //end of j loop
            
            } //end of i loop
            return mat;
            
            }
            ')

# calculate distance 
distance_kh_C <- function(alleles, prop = F, remove.NA = T){
  #  test for tibble??
  # prop = T ==> takes into account the number of loci available for comparison
  # prop = F ==> is the simple allelic difference. I would suggest this should only be used if comparing all loci ie remove.NA = F
  # ensure first column is Isolate
  if (colnames(alleles)[1] != "Isolate"){
    
    colnames(alleles)[1] <- "Isolate"
    
  }
  # size of matrix
  L <- length(alleles$Isolate)
  
  # set cols and the length of allele profile
  cols <- length(colnames(alleles))
  A <- c(2:cols)
  lenA <- length(A)
  # set all questionable alleled to LNF
  alleles[ alleles == "PLOT3"| alleles == "PLOT5" | alleles == "ASM"| alleles == "ALM" | alleles == "NIPH" | alleles == "NIPHEM" ] <- "LNF"
  # make alleles a matrix
  input <- as.matrix(alleles[,A])
  # find distance
  mat <- distance_kh(input,remove_na = remove.NA,prop = prop,"LNF")
  # # set matrix row and col names
  colnames(mat) <- alleles$Isolate
  rownames(mat) <- alleles$Isolate
  # return matrix
  mat
}



# source(file.path("{{script_dir}}",'pad.R'))
# get allele profile
alleles <- read_delim(file.path("overall_alleles.tsv"), "\t", escape_double =  F, trim_ws = T)
alleles$FILE <- gsub(".fa", "", alleles$FILE)
# calculate pairwise allelic distance
pad <- distance_kh_C(alleles = alleles, prop = F, remove.NA = T )
write.table(as.data.frame(pad) , file = file.path("distances.tab"), sep = "\t", row.names = T, quote = F)