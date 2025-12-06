#!/bin/bash


 pway=${1}
ipway=${2}


cd /scratch/mparis15/Goodarz/ZZZ_OLD_ASTHMA/
R --no-save --quiet << EOF
options( width=121 )
library( "RhpcBLASctl" )
# get_num_cores()
# get_num_procs()
blas_set_num_threads( 1 )

my <- read.table(
  "list_hugo_CAUSE.txt",
  sep="\t", header=TRUE,
  as.is=TRUE, stringsAsFactors=FALSE )


library("fgsea")
  strGO <- "/scratch/mparis15/PATHWAY_2024jun/Human_${pway}_fgsea.gmt"
   goBP <- gmtPathways( strGO )
gmtPathwaysDesc <- function( gmt.file ) 
{
 pathwayLines <- strsplit(readLines(gmt.file), "\t")
 H1 <- rep( "N/A", length( pathwayLines ) )
 for( i in 1:length( pathwayLines ) ) {
   H1[ i ] <- pathwayLines[ i ][[1]][1] }
 H2 <- rep( "N/A", length( H1 ) )
 for( i in 1:length( pathwayLines ) ) {
   H2[ i ] <- pathwayLines[ i ][[1]][2] }
 data.frame( pathway=H1, desc=H2, stringsAsFactors=FALSE )
}
desc <- gmtPathwaysDesc( strGO )
all_genes <- unique( do.call( c, goBP ) )


M <- merge( 
  data.frame( x=all_genes, stringsAsFactors=FALSE ), my, by=1 )

w_neg <- which( M[,2] < 0 )
w_pos <- which( M[,2] > 0 )


set.seed( ${ipway} )
 maxsim <- 1e8
lDF <- data.frame( stringsAsFactors=FALSE )
for( i in ${ipway} )
{
  genes <- goBP[[ i ]]
  # remove empty gene names (KEGG)
  genes <- genes[ nzchar( genes ) ]
  l_pway <- length( genes )
  if( l_pway > 999999 ) next


  #----- test causative (negative) direction
  w_in  <- intersect( w_neg, which( M[,1] %in% genes ) )
  l_in  <- length( w_in )
  w_out <- intersect( w_neg, setdiff( 1:nrow(M), w_in ) )
  l_out <- length( w_out )
  if( ( l_in >= 2 ) & ( l_out >= ( 2 * l_in ) ) )
  {
    s_in   <- sum( M[ w_in, 2 ] )
    i_sim  <- 1
    nb_low <- 0

    while( ( i_sim <= maxsim ) & 
           ( nb_low < 100    ) ) {
    w_rnd  <- sample( w_out, l_in )
    s_rnd  <- sum( M[ w_rnd, 2 ] )
    # count the number of times random simulation
    # was lower than that observed
    nb_low <- nb_low + ifelse( s_rnd < s_in, 1, 0 )
    i_sim  <- i_sim + 1
    }

    vec_in <- M[ w_in, 2 ]
   #names( vec_in ) <- M[ w_in, 1 ]
    names( vec_in ) <- paste0( M[ w_in, 1 ], "(", M[ w_in, 6 ], ")" )
    vec_in <- sort( vec_in, decreasing=( s_in > 0 ) )
    pval   <- ( 1 + nb_low ) / ( 1 + i_sim )

    df <- data.frame( 
        pway=desc[i,1],
        onto="${pway}",
        desc=desc[i,2],
         dir="cause",
         sum=s_in,
        pval=pval,
        padj=1,
        l_in=l_in,
       l_out=l_out,
      l_pway=l_pway,
         idx=i,
leading_edge=paste( 
             head( names( vec_in ), min( l_in, 25 ) ),
             collapse = ',' ),
     stringsAsFactors=FALSE )
     # print( df )
  if( 0 == nrow(lDF) ){ lDF <- df } else { lDF <- rbind( lDF, df ) }
  }


  #----- test pleiotropic (positive) direction
  w_in  <- intersect( w_pos, which( M[,1] %in% genes ) )
  l_in  <- length( w_in )
  w_out <- intersect( w_pos, setdiff( 1:nrow(M), w_in ) )
  l_out <- length( w_out )
  if( ( l_in >= 2 ) & ( l_out >= ( 2 * l_in ) ) )
  {
    s_in   <- sum( M[ w_in, 2 ] )
    i_sim  <- 1
    nb_low <- 0

    while( ( i_sim <= maxsim ) & 
           ( nb_low < 100    ) ) {
    w_rnd  <- sample( w_out, l_in )
    s_rnd  <- sum( M[ w_rnd, 2 ] )
    # count the number of times random simulation
    # was higher than that observed
    nb_low <- nb_low + ifelse( s_rnd > s_in, 1, 0 )
    i_sim  <- i_sim + 1
    # if( 0 == ( i_sim %% 1000 ) ){ 
    #   cat( "pleio", i_sim, l_in, l_out, "\n", sep=" " ) }
    }

    vec_in <- M[ w_in, 2 ]
   #names( vec_in ) <- M[ w_in, 1 ]
    names( vec_in ) <- paste0( M[ w_in, 1 ], "(", M[ w_in, 6 ], ")" )
    vec_in <- sort( vec_in, decreasing=( s_in > 0 ) )
    pval   <- ( 1 + nb_low ) / ( 1 + i_sim )

    df <- data.frame( 
        pway=desc[i,1],
        onto="${pway}",
        desc=desc[i,2],
         dir="pleio",
         sum=s_in,
        pval=pval,
        padj=1,
        l_in=l_in,
       l_out=l_out,
      l_pway=l_pway,
         idx=i,
leading_edge=paste( 
             head( names( vec_in ), min( l_in, 25 ) ),
             collapse = ',' ),
     stringsAsFactors=FALSE )
    # print( df )
  if( 0 == nrow(lDF) ){ lDF <- df } else { lDF <- rbind( lDF, df ) }
  }
  }


write.table( lDF, file="pway_2024jun_${pway}_${ipway}.txt",
  sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE )
quit()
EOF
