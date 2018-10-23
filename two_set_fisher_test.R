# two_set_fisher_test.R 

# commonly used stats functions
two_set_fisher_test = function(universe, setA, setB, alternative = "greater"){
  p = NA; OR = NA
  
  fisher_elements = c(sum(!universe %in% c(setA,setB)),sum((universe %in% setA) & !(universe %in% setB)),
                      sum(!(universe %in% setA) & (universe %in% setB)),sum((universe %in% setA) & (universe %in% setB)))
  if (fisher_elements[1] > 0 && fisher_elements[3] > 0 && fisher_elements[2] >= 0 && fisher_elements[4] >= 0){
    test.table = matrix(as.numeric(fisher_elements), nrow=2)
    f.test = fisher.test(test.table, alternative = alternative) # only test for enrichment for now
    OR = f.test$estimate
    p = f.test$p.value
    
    count00 = test.table[1,1]
    count10 = test.table[2,1]
    count01 = test.table[1,2]
    count11 = test.table[2,2]
    
    # return(list("p"=p, "OR"=OR, "notAnotB" = count00
    #             , "AnotB" = count10, "BnotA" = count01, "AB" = count11))
    
    result_row = c(fisher_elements,p,OR)
    return(result_row)
  }
}