setGeneric(
  "run_saem",
  function(data, model, init, tuning_algo, hyperparam) {
    standardGeneric("run_saem")
  }
)

setMethod(
  "run_saem",
  signature(
    data = "saemvsData", model = "saemvsModel", init = "saemvsInit",
    tuning_algo = "saemvsTuning", hyperparam = "saemvsHyperSpikeAndSlab"
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
          config, state, sa_step_hdim, m_step_full_map, update_proposal_mh_to_select
        )
      },
      map_part_select_nofixed = {
        state <- run_generic_saem(
          config, state, sa_step_split, m_step_mix, update_proposal_mh_all
        )
      },
      map_part_select_fixed = {
        state <- run_fixed_saem(
          config, state, sa_step_split, m_step_mix, update_proposal_mh_all
        )
      },
      mle_nofixed = {
        state <- run_generic_saem(
          config, state, sa_step_ldim, m_step_mle, update_proposal_mh_not_to_select
        )
      },
      mle_fixed = {
        state <- run_fixed_saem(
          config, state, sa_step_ldim, m_step_mle, update_proposal_mh_not_to_select
        )
      }
    )

    return(state)
  }
)
