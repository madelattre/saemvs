#' @export
setGeneric(
  "summary",
  function(res) {
    standardGeneric("summary")
  }
)

#' @exportMethod summary
setMethod(
  "summary",
  signature(
    res = "resSAEMVS"
  ),
  function(res) {
    # Vérifier que res n'est pas vide

    selected_support <- res@unique_support[[which.min(res@crit_values)]]
    
    dim_s <- dim(selected_support)

    print("---- Selected variables : ----")

    for (j in seq(dim_s[2])) {
        select_var <- which(selected_support[,j] == TRUE)
        print(paste0("- Parameter ", j, " :"))
        print(select_var)
    }
  }
)