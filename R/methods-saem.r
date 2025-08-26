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
        state <- run_saem_full(
          config, state, sa_step_to_select, m_step_map_to_select, update_proposal_mh_to_select
        )
      },
      map_part_select_nofixed = {
        state <- run_saem_full(
          config, state, sa_step_all, m_step_map_all, update_proposal_mh_all
        )
      },
      map_part_select_fixed = {
        state <- run_saem_with_fixed_effects(
          config, state, sa_step_all, m_step_map_all, update_proposal_mh_all
        )
      },
      mle_nofixed = {
        state <- run_saem_full(
          config, state, sa_step_not_to_select, m_step_mle, update_proposal_mh_not_to_select
        )
      },
      mle_fixed = {
        state <- run_saem_with_fixed_effects(
          config, state, sa_step_not_to_select, m_step_mle, update_proposal_mh_not_to_select
        )
      }
    )

    return(state)
  }
)
