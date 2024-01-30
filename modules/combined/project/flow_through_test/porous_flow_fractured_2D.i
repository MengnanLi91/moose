# geometry
# core sample
# 8in(0.2032) x 1.5in(0.0381) Cylinder
radius = 0.01905 #m
length = 0.2032 #m
# fracture
# 0.24 mm along the axis of the core
frac_len = 2.4e-4 # m
# initial condition
inlet_pressure = 1.78e6
outlet_pressure = 1.63e6
domain_porosity = 0.127
domain_permeability = 2.10E-19 # 3D 7.96E-20
frac_porosity = ${fparse  domain_porosity*frac_len}
frac_permeability = ${fparse frac_len*frac_len*frac_len/12}
Al_mass = 6.147940e-07
Ca_mass = 2.272036e-05
Cl_mass = 4.537986e-06
F_mass =  1.875223e-07
Fe_mass =  2.356492e-12
H_mass = 3.392909e-04
K_mass =  3.127768e-06
Mg_mass =  1.549918e-05
NO3_mass = 3.947198e-07
Na_mass = 8.083672e-06
O2_mass = 6.412920e-06
SO4_mass = 1.513277e-05
SiO2_mass = 8.907870e-06
HCO3_mass = 2.068976e-02
H2O_mass = 9.788853E-01
porepressure = ${inlet_pressure}
end_time = 259200
density = 1010 # kg/m^3

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 0'
[]


[Mesh]
  [cmg]
    type = CartesianMeshGenerator
    dim = 2
    dx = '${radius} ${frac_len} ${radius}'
    dy = '${length}'
    ix = '10 1 10'
    iy = '80'
    subdomain_id = '1 2 3'
  []

  [matrix_subdomain]
    type = RenameBlockGenerator
    input = cmg
    old_block = '1 3'
    new_block = 'matrix matrix'
  []

  [fracture_subdomain]
    type = RenameBlockGenerator
    input = matrix_subdomain
    old_block = 2
    new_block = fracture
  []

  [rename]
    type = RenameBoundaryGenerator
    input = fracture_subdomain
    old_boundary = '0 1 2 3'
    new_boundary = 'inlet right_wall outlet left_wall'
  []

  [break_boundary]
    input = rename
    type = BreakBoundaryOnSubdomainGenerator
    boundaries = 'inlet outlet'
  []

  [rename_fracture]
    type = RenameBoundaryGenerator
    input = break_boundary
    old_boundary = '4 6'
    new_boundary = 'frac_inlet frac_outlet'
  []

[]

[ICs]
  [f_Al]
    type = ConstantIC
    variable = f_Al
    value = ${Al_mass}
    block = fracture
  []

  [f_Al_matrix]
    type = ConstantIC
    variable = f_Al
    value = 0
    block = matrix
  []

  [f_HCO3]
    type = ConstantIC
    variable = f_HCO3
    value = ${HCO3_mass}
    block = fracture
  []

  [f_HCO3_matrix]
    type = ConstantIC
    variable = f_HCO3
    value = 0
    block = matrix
  []
[]

[Variables]
  [f_Al]
#    initial_condition = ${Al_mass}
  []
  [f_HCO3]
#    initial_condition = ${HCO3_mass}
  []
  [porepressure]
    initial_condition = ${porepressure}
  []
[]

[AuxVariables]
  [rate_Al]
  []
  [rate_HCO3]
  []
  [rate_H2O]
  []
  [chem_porosity]
    initial_condition = ${domain_porosity}
  []
  [velocity_x]
    family = MONOMIAL
    order = CONSTANT
    block = fracture
  []
  [velocity_y]
    family = MONOMIAL
    order = CONSTANT
    block = fracture
  []
[]


[AuxKernels]
  [velocity_x]
    type = PorousFlowDarcyVelocityComponent
    variable = velocity_x
    component = x
  []
  [velocity_y]
    type = PorousFlowDarcyVelocityComponent
    variable = velocity_y
    component = y
  []
[]


[BCs]
  [constant_injection_porepressure]
    type = DirichletBC
    variable = porepressure
    value = ${inlet_pressure}
    boundary = inlet
  []
  [constant_outer_porepressure]
    type = DirichletBC
    variable = porepressure
    value = ${outlet_pressure}
    boundary = outlet
  []

  [f_Al_in]
    type = DirichletBC
    variable = f_Al
    boundary = frac_inlet
    value = ${Al_mass}
  []
  [f_HCO3_in]
    type = DirichletBC
    variable = f_HCO3
    boundary = frac_inlet
    value = ${HCO3_mass}
  []
  # [f_H2O_in]
  #   type = DirichletBC
  #   variable = porepressure
  #   boundary = inlet_node
  #   value = ${H2O_mass}
  # []

  [f_Al]
    type = PorousFlowOutflowBC
    variable = f_Al
    boundary = frac_outlet
    include_relperm = false
    mass_fraction_component = 0
  []
  [f_HCO3]
    type = PorousFlowOutflowBC
    variable = f_HCO3
    boundary = frac_outlet
    include_relperm = false
    mass_fraction_component = 1
  []
  # [f_Al]
  #   type = PorousFlowOutflowBC
  #   variable = f_Al
  #   boundary = outlet
  #   include_relperm = false
  #   mass_fraction_component = 0
  # []
  # [f_HCO3]
  #   type = PorousFlowOutflowBC
  #   variable = f_HCO3
  #   boundary = outlet
  #   include_relperm = false
  #   mass_fraction_component = 1
  # []
  # [f_H2O]
  #   type = PorousFlowOutflowBC
  #   variable = porepressure
  #   boundary = outlet
  #   include_relperm = false
  # []

[]

# [DiracKernels]
#   [outflow_Al]
#     type = PorousFlowPointSourceOutflow
#     point = ' 0 ${length} 0'
#     variable = f_Al
#     multiply_by_density = false
#     include_relperm = false
#   []

#   [outflow_HCO3]
#     type = PorousFlowPointSourceOutflow
#     point = ' 0 ${length} 0'
#     variable = f_HCO3
#     multiply_by_density = false
#     include_relperm = false
#   []
# []

[Postprocessors]
  # [nnn]
  #   type = NearestNodeNumber
  #   nearest_node_number_uo = nnn_uo
  #   execute_on = 'initial timestep_begin'
  # []
  [mass_extracted_Al_in]
    type = SideAverageValue
    variable = f_Al
    boundary = frac_inlet
    execute_on = 'initial timestep_end'
  []
  [mass_extracted_HCO3_in]
    type = SideAverageValue
    variable = f_HCO3
    boundary = frac_inlet
    execute_on = 'initial timestep_end'
  []
  [mass_extracted_Al_out]
    type = SideAverageValue
    variable = f_Al
    boundary = frac_outlet
    execute_on = 'initial timestep_end'
  []
  [mass_extracted_HCO3_out]
    type = SideAverageValue
    variable = f_HCO3
    boundary = frac_outlet
    execute_on = 'initial timestep_end'
  []
  [darcy_velocity]
    type = SideAverageValue
    variable = darcy_vel_y
    boundary = inlet
    execute_on = 'initial timestep_end'
  []
  [mass_extracted]
    type = LinearCombinationPostprocessor
    pp_names = 'mass_extracted_Al_out mass_extracted_HCO3_out'
    pp_coefs = '1 1'
    execute_on = 'initial timestep_end'
  []
  [avg_porosity]
    type = ElementAverageValue
    variable = chem_porosity
  []
[]

[FluidProperties]
  [the_simple_fluid]
    type = SimpleFluidProperties
    bulk_modulus = 2E9
    viscosity = 1.03E-03
    density0 = 1010
  []
[]

[PorousFlowFullySaturated]
  coupling_type = Hydro
  porepressure = porepressure
  mass_fraction_vars = 'f_Al f_HCO3'
  save_component_rate_in = 'rate_Al rate_HCO3 rate_H2O' # change in kg at every node / dt
  fp = the_simple_fluid
  temperature_unit = Celsius
  pressure_unit = Pa
  #multiply_by_density = true
  add_darcy_aux = true
  stabilization = Full
[]

[Materials]
  [poro_fracture]
    type = PorousFlowPorosityConst
    porosity = ${frac_porosity}
    block = 'fracture'
  []
  [poro_matrix]
    type = PorousFlowPorosityChemCoupled
    porosity = chem_porosity
    block = 'matrix'
  []
  [permeability_matrix]
    type = PorousFlowPermeabilityConst
    permeability = '${domain_permeability} 0 0   0 ${domain_permeability} 0   0 0 ${domain_permeability}'
    block = 'matrix'
  []
  [permeability_fracture]
    type = PorousFlowPermeabilityConst
    permeability = '${frac_permeability} 0 0  0 ${frac_permeability} 0  0 0 ${frac_permeability}'
    block = 'fracture'
  []
[]

[Preconditioning]
  active = typically_efficient
  [typically_efficient]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_hypre_type'
    petsc_options_value = ' hypre    boomeramg'
  []
  [strong]
    type = SMP
    full = true
    petsc_options = '-ksp_diagonal_scale -ksp_diagonal_scale_fix'
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap'
    petsc_options_value = ' asm      ilu           NONZERO                   2'
  []
  [probably_too_strong]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
    petsc_options_value = ' lu       mumps'
  []
[]

[Executioner]
  type = Transient
  solve_type = Newton
  end_time = ${end_time} #126092.85
  nl_rel_tol = 1E-8
  nl_abs_tol = 1e-12
  [TimeSteppers]
    [funcDT]
      type = FunctionDT
      function = 'if(t>500, 100, 10)'
    []
  []
[]


[Outputs]
  exodus = true
  csv = true
  # [my_checkpoint]
  #   type = Checkpoint
  #   num_files = 4
  #   interval = 5
  # []
[]

# [MultiApps]
#   [react]
#     type = TransientMultiApp
#     input_files = basalt_mineral_test.i
#     clone_master_mesh = true
#     execute_on = 'timestep_end'
#   []
# []
# [Transfers]
#   [changes_due_to_flow]
#     type = MultiAppCopyTransfer
#     source_variable = 'rate_Al rate_Ca rate_Cl rate_F rate_Fe rate_H rate_HCO3 rate_K rate_Mg rate_Na rate_NO3 rate_O2 rate_SiO2 rate_SO4 rate_H2O porepressure'
#     variable = 'pf_rate_Al pf_rate_Ca pf_rate_Cl pf_rate_F pf_rate_Fe pf_rate_H pf_rate_HCO3 pf_rate_K pf_rate_Mg pf_rate_Na pf_rate_NO3 pf_rate_O2 pf_rate_SiO2 pf_rate_SO4 pf_rate_H2O pressure'
#     to_multi_app = react
#   []
#   [massfrac_from_geochem]
#     type = MultiAppCopyTransfer
#     source_variable = 'massfrac_Al massfrac_Ca massfrac_Cl massfrac_F massfrac_Fe massfrac_H massfrac_HCO3 massfrac_K massfrac_Mg massfrac_Na massfrac_NO3 massfrac_O2 massfrac_SiO2 massfrac_SO4 porosity'
#     variable = 'f_Al f_Ca f_Cl f_F f_Fe f_H f_HCO3 f_K f_Mg f_Na f_NO3 f_O2 f_SiO2 f_SO4 chem_porosity'
#     from_multi_app = react
#   []
# []
