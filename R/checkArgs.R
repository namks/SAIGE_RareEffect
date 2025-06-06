checkArgBool = function(arg0, arg0name){
  if(!arg0 %in% c(TRUE, FALSE))
    stop(arg0name, " should be TRUE or FALSE.")
}	

checkArgNumeric = function(arg0, arg0name, minVal = NULL, maxVal = NULL, incMin = TRUE, incMax = TRUE){
 if(!is.numeric(arg0)){
     stop(arg0name, " should be a numeric value.")
  }	 

  if(!is.null(minVal)){
   if(incMin){  
    if(arg0 < minVal){
      stop(arg0name, " should be a numeric value greater than or equal to ", minVal,"\n")
    }
   }else{ 
    if(arg0 <= minVal){
      stop(arg0name, " should be a numeric value greater than", minVal,"\n")
    }

   }	   
  }

  if(!is.null(maxVal)){
   if(incMax){	  
    if(arg0 > maxVal){
      stop(arg0name, " should be a numeric value less than or equal to ", maxVal,"\n")
    }
   }else{
     if(arg0 >= maxVal){
	   stop(arg0name, " should be a numeric value less than", maxVal,"\n")  	   
     } 	     
   }
  }
}




checkArgsListNumeric =function(start,
end,
max_missing,
min_MAC,
min_MAF,
min_Info,
SPAcutoff,
dosage_zerod_cutoff,
dosage_zerod_MAC_cutoff,
markers_per_chunk, 
groups_per_chunk,
minGroupMAC_in_BurdenTest,
max_MAC_use_ER){

	checkArgNumeric(start, deparse(substitute(start)), 1, 250000000)
        checkArgNumeric(end, deparse(substitute(end)), 1, 250000000)
	checkArgNumeric(max_missing, deparse(substitute(max_missing,)), 0, 1)
	checkArgNumeric(min_MAC, deparse(substitute(min_MAC,)), minVal=0.5, incMin=T)
	checkArgNumeric(min_MAF, deparse(substitute(min_MAF,)), 0, 0.5)
	checkArgNumeric(min_Info, deparse(substitute(min_Info,)), 0, 1)
	checkArgNumeric(SPAcutoff, deparse(substitute(SPAcutoff)), 0.5, 4)
	checkArgNumeric(dosage_zerod_cutoff, deparse(substitute(dosage_zerod_cutoff)), 0, 0.5)
	checkArgNumeric(dosage_zerod_MAC_cutoff, deparse(substitute(dosage_zerod_MAC_cutoff)), dosage_zerod_cutoff, 100)
        cat("Any dosages <= ", dosage_zerod_cutoff, " for genetic variants with MAC <= ", dosage_zerod_MAC_cutoff, " are set to be 0.\n")
	checkArgNumeric(markers_per_chunk, deparse(substitute(markers_per_chunk)), minVal=1000)
	checkArgNumeric(groups_per_chunk, deparse(substitute(groups_per_chunk)), minVal=1)
	checkArgNumeric(max_MAC_use_ER, deparse(substitute(max_MAC_use_ER)), minVal=0)
	#cat("min_MAC ", min_MAC, "\n")
	#if(length(which(minGroupMAC_in_BurdenTest > 0)) > 0){
	#	cat("minGroupMAC_in_BurdenTest need to be >= min_MAC\n")
	#}
	
}




checkArgsListBool = function(is_imputed_data,
                     LOCO,
                     is_output_moreDetails,
		     is_overwrite_output){
                     #is_rewrite_XnonPAR_forMales){

	checkArgBool(is_imputed_data, deparse(substitute(is_imputed_data)))
	checkArgBool(LOCO, deparse(substitute(LOCO)))
	checkArgBool(is_output_moreDetails, deparse(substitute(is_output_moreDetails)))
	checkArgBool(is_overwrite_output, deparse(substitute(is_overwrite_output)))
	#checkArgBool(is_rewrite_XnonPAR_forMales, deparse(substitute(is_rewrite_XnonPAR_forMales)))		
}



#checkArgsList_for_Region = function(method_to_CollapseUltraRare, 
				    #DosageCutoff_for_UltraRarePresence,
checkArgsList_for_Region = function( 
				    MACCutoff_to_CollapseUltraRare,
				    maxMAF_in_groupTest,
				    maxMAC_in_groupTest,
				    markers_per_chunk_in_groupTest){


    #if (method_to_CollapseUltraRare != "") {
        if (MACCutoff_to_CollapseUltraRare < 0) {
            stop("MACCutoff_to_CollapseUltraRare must be >= 0\n")
        }
        #if (DosageCutoff_for_UltraRarePresence <= 0 | DosageCutoff_for_UltraRarePresence > 2) {
        #    stop("DosageCutoff_for_UltraRarePresenc needs be to larger than 0 and less or equal to 2\n")
        #}

        #if (method_to_CollapseUltraRare == "absence_or_presence") {
                  cat("Ultra rare variants with MAC <= ", MACCutoff_to_CollapseUltraRare,
                    " will be collpased for set-based tests in the 'absence or presence' way.\n")
                    
		    
		    #"For the resulted collpased marker, any individual having ",
                    #DosageCutoff_for_UltraRarePresence, "<= dosage < ",
                    #(1 + DosageCutoff_for_UltraRarePresence),
                    #" for any ultra rare variant has 1 in the genotype vector, having dosage >= ",
                    #(1 + DosageCutoff_for_UltraRarePresence),
                    #" for any ultra rare variant has 2 in the genotype vector, otherwise 0. \n")
        #}
	#else if (method_to_CollapseUltraRare == "sum_geno") {
        #          cat("Ultra rare variants with MAC <= ", MACCutoff_to_CollapseUltraRare,
        #            " will be collpased for set-based tests in the 'sum_geno' way. ",
        #            "The resulted collpased marker equals weighted sum of the genotypes of all ultra rare variantsi. NOTE: this option currently is not active\n")
        #}
    #}
    #else {
    #    cat("Ultra rare variants won't be collpased for set-based association tests\n")
    #}

    if(length(maxMAF_in_groupTest) < 1){
      stop("maxMAF_in_groupTest should contain at least one numeric value between 0 and 0.5\n")
    }else{	    
      for(i in 1:length(maxMAF_in_groupTest)){
        checkArgNumeric(maxMAF_in_groupTest[i], deparse(substitute(maxMAF_in_groupTest[i])), 0, 0.5, FALSE, TRUE)
      } 
    }

    
    if(length(maxMAC_in_groupTest) < 1){
      stop("maxMAC_in_groupTest should contain at least one numeric value >= 0\n")
    }else{	    
      for(i in 1:length(maxMAC_in_groupTest)){
        checkArgNumeric(maxMAC_in_groupTest[i], deparse(substitute(maxMAC_in_groupTest[i])),  minVal=0)
      } 
    }

    checkArgNumeric(markers_per_chunk_in_groupTest, deparse(substitute(markers_per_chunk_in_groupTest)), 1, 10000, TRUE, TRUE)
}

