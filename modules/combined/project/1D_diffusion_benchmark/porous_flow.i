#########################################
#                                       #
# File written by create_input_files.py #
#                                       #
#########################################
# PorousFlow simulation of injection and production in a simplified GeoTES aquifer
# Much of this file is standard porous-flow stuff.  The unusual aspects are:
# - transfer of the rates of changes of each species (kg.s) to the aquifer_geochemistry.i simulation.  This is achieved by saving these changes from the PorousFlowMassTimeDerivative residuals
# - transfer of the temperature field to the aquifer_geochemistry.i simulation
# Interesting behaviour can be simulated by this file without its 'parent' simulation, exchanger.i.  exchanger.i provides mass-fractions injected via the injection_rate_massfrac_* variables, but since these are more-or-less constant throughout the duration of the exchanger.i simulation, the initial_conditions specified below may be used.  Similar, exchanger.i provides injection_temperature, but that is also constant.

[Mesh]
  [cmg]
    type = CartesianMeshGenerator
    dim = 1
    dx = '0.5 1 1.5 1.5 1 5 10 20 2.5'
    ix = '1 10 30 30 10 10 10 10 1'
    subdomain_id = '1 2 3 4 5 6 7 8 9'
  []

  [concrete]
    type = RenameBlockGenerator
    input = cmg
    old_block = '1 2 3'
    new_block = '1 1 1'
  []

  [rename_concrete]
    type = RenameBlockGenerator
    input = concrete
    old_block = '1'
    new_block = 'concrete'
  []

  [clayey]
    type = RenameBlockGenerator
    input = rename_concrete
    old_block = '4 5 6 7 8 9'
    new_block = '2 2 2 2 2 2'
  []

  [rename_clayey]
    type = RenameBlockGenerator
    input = clayey
    old_block = '2'
    new_block = 'clayey'
  []
[]

[Problem]
  type = FEProblem
  coord_type = 'RZ'
  rz_coord_axis = Y
[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 0'
[]

# add BC conditions for the rest species
[BCs]
  [constant_injection_porepressure]
    type = DirichletBC
    variable = porepressure
    value = 10
    boundary = left
  []
  [constant_outer_porepressure]
    type = DirichletBC
    variable = porepressure
    value = 10
    boundary = right
  []

  [f_H_in]
    type = DirichletBC
    variable = f_H
    boundary = left
    value = 0.1
  []

  [f_H]
    type = PorousFlowOutflowBC
    variable = f_H
    boundary = right
    include_relperm = false
    mass_fraction_component = 0
  []

  # [f_H20_in]
  #   type = DirichletBC
  #   variable = porepressure
  #   boundary = left
  #   value = 0.9
  # []

  # [f_H20]
  #   type = PorousFlowOutflowBC
  #   variable = porepressure
  #   boundary = right
  #   include_relperm = false
  #   mass_fraction_component = 1
  # []
[]

[Modules]
  [./FluidProperties]
    [./the_simple_fluid]
      type = SimpleFluidProperties
      thermal_expansion = 0
      bulk_modulus = 2E9
      viscosity = 1E-3
      density0 = 1000
      cv = 4000.0
      cp = 4000.0
    [../]
  [../]
[]

[PorousFlowFullySaturated]
  coupling_type = Hydro
  porepressure = porepressure
  mass_fraction_vars = 'f_H'
  save_component_rate_in = 'rate_H rate_H2O'
  fp = the_simple_fluid
  temperature_unit = Celsius
  pressure_unit = Pa
  #multiply_by_density = true
  add_darcy_aux = true
  stabilization = Full
[]

[Kernels]
  [f_H]
    type = PorousFlowDispersiveFlux
  []
[]

[Materials]
  [./porosity_clayey]
    type = PorousFlowPorosityConst # this simulation has no porosity changes from dissolution
    block = clayey
    porosity = 0.01
  [../]
  [./porosity_concrete]
    type = PorousFlowPorosityConst # this simulation has no porosity changes from dissolution
    block = concrete
    porosity = 0.063
  [../]
  [./permeability_clayey]
    type = PorousFlowPermeabilityConst
    block = clayey
    permeability = '1E-18 0 0   0 1E-18 0   0 0 1E-18'
  [../]
  [./permeability_concrete]
    type = PorousFlowPermeabilityConst
    block = concrete
    permeability = '1.7E-15 0 0   0 1.7E-15 0   0 0 4.1E-16'
  [../]
[]

[Preconditioning]
  active = typically_efficient
  [./typically_efficient]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_hypre_type'
    petsc_options_value = ' hypre    boomeramg'
  [../]
  [./strong]
    type = SMP
    full = true
    petsc_options = '-ksp_diagonal_scale -ksp_diagonal_scale_fix'
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap'
    petsc_options_value = ' asm      ilu           NONZERO                   2'
  [../]
  [./probably_too_strong]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
    petsc_options_value = ' lu       mumps'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = Newton
  end_time = 7.76E6 # 90 days
  [./TimeStepper]
    type = FunctionDT
    function = 'min(3E4, max(1E4, 0.2 * t))'
  [../]
[]
[Outputs]
  exodus = true
[]

[ICs]
  [f_H]
    type = ConstantIC
    block = concrete
    value = 0.1
    variable = f_H
  []
[]

[Variables]
  [./f_H]
  #  initial_condition = 0.1
  [../]
  # [./f_Cl]
  #   initial_condition = 0.04870664551708
  # [../]
  # [./f_SO4]
  #   initial_condition = 0.0060359986852517
  # [../]
  # [./f_HCO3]
  #   initial_condition = 5.0897287594019e-05
  # [../]
  # [./f_SiO2aq]
  #   initial_condition = 3.0246609868421e-05
  # [../]
  # [./f_Al]
  #   initial_condition = 3.268028901929e-08
  # [../]
  # [./f_Ca]
  #   initial_condition = 0.00082159428184586
  # [../]
  # [./f_Mg]
  #   initial_condition = 1.8546347062146e-05
  # [../]
  # [./f_Fe]
  #   initial_condition = 4.3291908204093e-05
  # [../]
  # [./f_K]
  #   initial_condition = 6.8434768308898e-05
  # [../]
  # [./f_Na]
  #   initial_condition = 0.033298053919671
  # [../]
  # [./f_Sr]
  #   initial_condition = 1.2771866652177e-05
  # [../]
  # [./f_F]
  #   initial_condition = 5.5648860174073e-06
  # [../]
  # [./f_BOH]
  #   initial_condition = 0.0003758574621917
  # [../]
  # [./f_Br]
  #   initial_condition = 9.0315286107068e-05
  # [../]
  # [./f_Ba]
  #   initial_condition = 1.5637460875161e-07
  # [../]
  # [./f_Li]
  #   initial_condition = 8.3017067912701e-05
  # [../]
  # [./f_NO3]
  #   initial_condition = 0.00010958455036169
  # [../]
  # [./f_O2aq]
  #   initial_condition = -7.0806852373351e-05
  # [../]
  [./porepressure]
    initial_condition = 30E4
  [../]
[]

[AuxVariables]
  [./rate_H]
  [../]
  # [./rate_Cl]
  # [../]
  # [./rate_SO4]
  # [../]
  # [./rate_HCO3]
  # [../]
  # [./rate_SiO2aq]
  # [../]
  # [./rate_Al]
  # [../]
  # [./rate_Ca]
  # [../]
  # [./rate_Mg]
  # [../]
  # [./rate_Fe]
  # [../]
  # [./rate_K]
  # [../]
  # [./rate_Na]
  # [../]
  # [./rate_Sr]
  # [../]
  # [./rate_F]
  # [../]
  # [./rate_BOH]
  # [../]
  # [./rate_Br]
  # [../]
  # [./rate_Ba]
  # [../]
  # [./rate_Li]
  # [../]
  # [./rate_NO3]
  # [../]
  # [./rate_O2aq]
  # [../]
  [./rate_H2O]
  [../]
[]

# [MultiApps]
#   [./react]
#     type = TransientMultiApp
#     input_files = aquifer_geochemistry.i
#     clone_master_mesh = true
#     execute_on = 'timestep_end'
#   [../]
# []
# [Transfers]
#   [./changes_due_to_flow]
#     type = MultiAppCopyTransfer
#     direction = to_multiapp
#     source_variable = 'rate_H rate_Cl rate_SO4 rate_HCO3 rate_SiO2aq rate_Al rate_Ca rate_Mg rate_Fe rate_K rate_Na rate_Sr rate_F rate_BOH rate_Br rate_Ba rate_Li rate_NO3 rate_O2aq rate_H2O temperature'
#     variable = 'pf_rate_H pf_rate_Cl pf_rate_SO4 pf_rate_HCO3 pf_rate_SiO2aq pf_rate_Al pf_rate_Ca pf_rate_Mg pf_rate_Fe pf_rate_K pf_rate_Na pf_rate_Sr pf_rate_F pf_rate_BOH pf_rate_Br pf_rate_Ba pf_rate_Li pf_rate_NO3 pf_rate_O2aq pf_rate_H2O temperature'
#     multi_app = react
#   [../]
#   [./massfrac_from_geochem]
#     type = MultiAppCopyTransfer
#     direction = from_multiapp
#     source_variable = 'massfrac_H massfrac_Cl massfrac_SO4 massfrac_HCO3 massfrac_SiO2aq massfrac_Al massfrac_Ca massfrac_Mg massfrac_Fe massfrac_K massfrac_Na massfrac_Sr massfrac_F massfrac_BOH massfrac_Br massfrac_Ba massfrac_Li massfrac_NO3 massfrac_O2aq '
#     variable = 'f_H f_Cl f_SO4 f_HCO3 f_SiO2aq f_Al f_Ca f_Mg f_Fe f_K f_Na f_Sr f_F f_BOH f_Br f_Ba f_Li f_NO3 f_O2aq '
#     multi_app = react
#   [../]
# []
