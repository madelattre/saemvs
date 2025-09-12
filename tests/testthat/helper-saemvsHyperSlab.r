# ------------------------------------------------------------
# Helper for saemvsHyperSlab / saemvsHyperSpikeAndSlab tests
# ------------------------------------------------------------

# Valid object for saemvsHyperSlab
valid_hyper_slab <- list(
  slab_parameter = 12000,
  cov_re_prior_scale = diag(2),
  cov_re_prior_df = 2
)

# Invalid cases for saemvsHyperSlab
invalid_cases_slab <- list(
  "negative slab_parameter" = list(
    args = modifyList(valid_hyper_slab, list(slab_parameter = -1)),
    msg = "slab_parameter must be NULL or strictly positive."
  ),
  "zero cov_re_prior_df" = list(
    args = modifyList(valid_hyper_slab, list(cov_re_prior_df = 0)),
    msg = "cov_re_prior_df must be NULL or strictly positive."
  ),
  "cov_re_prior_scale non-square" = list(
    args = modifyList(valid_hyper_slab, list(
      cov_re_prior_scale = matrix(1:6, nrow = 2, ncol = 3)
    )),
    msg = "cov_re_prior_scale"
  ),
  "cov_re_prior_scale non-numeric" = list(
    args = modifyList(valid_hyper_slab, list(
      cov_re_prior_scale = matrix(c("a", "b", "c", "d"), nrow = 2)
    )),
    msg = "cov_re_prior_scale"
  )
)

# Valid object for saemvsHyperSpikeAndSlab
valid_hyper_spike <- list(
  spike_parameter = 0.1,
  hyper_slab = do.call(saemvsHyperSlab, valid_hyper_slab)
)

# Invalid cases for saemvsHyperSpikeAndSlab
invalid_cases_spike <- list(
  "negative spike_parameter" = list(
    args = modifyList(valid_hyper_spike, list(spike_parameter = -0.1)),
    msg = "spike_parameter must be NULL or strictly positive."
  ),
  "zero spike_parameter" = list(
    args = modifyList(valid_hyper_spike, list(spike_parameter = 0)),
    msg = "spike_parameter must be NULL or strictly positive."
  )
)
