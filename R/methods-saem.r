setGeneric(
  "run_saem",
  function(data, model, init, tuning_algo, hyperparam) {
    standardGeneric("run_saem")
  }
)

setMethod(
  "run_saem",
  signature(
    data = "saemvsData", model = "modelC", init = "initC",
    tuning_algo = "tuningC", hyperparam = "fullHyperC"
  ),
  function(data, model, init, tuning_algo, hyperparam) {
    # Mise en forme et vérification de la cohérence des objets
    # Attention, l'ordre de ces étapes est important

    check_data(data, model)
    data_processed <- prepare_data(data, model)
    # if not default_init check_init(init,data,model)
    # else create init object with the default method
    check_init(init, data_processed, model)
    init_alg <- prepare_init(init, model)
    check_hyper(hyperparam, model, tuning_algo)
    hyperparam <- prepare_hyper(hyperparam, data_processed, model)
    config <- make_config(
      data_processed, model, tuning_algo, init_alg, hyperparam
    )
    state <- init_state(config)
    case <- get_case(config)

    # Utiliser les bonnes implémentations des SA et M steps selon les cas

    switch(case,
      map_full_select = {
        state <- run_generic_saem(
          config, state, sa_step_hdim, m_step_full_map, update_prop_mh_hdim
        )
      },
      map_part_select_nofixed = {
        state <- run_generic_saem(
          config, state, sa_step_split, m_step_mix, update_prop_mh_mix
        )
      },
      map_part_select_fixed = {
        state <- run_fixed_saem(
          config, state, sa_step_split, m_step_mix, update_prop_mh_mix
        )
      },
      mle_nofixed = {
        state <- run_generic_saem(
          config, state, sa_step_ldim, m_step_mle, update_prop_mh_ldim
        )
      },
      mle_fixed = {
        state <- run_fixed_saem(
          config, state, sa_step_ldim, m_step_mle, update_prop_mh_ldim
        )
      }
    )

    return(state)
  }
)
