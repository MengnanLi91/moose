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
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 20
    ny = 80
    xmin = ${fparse -radius}
    xmax = ${radius} #m
    ymin = 0
    ymax = ${fparse length} #m
  []
  [rename]
    type = RenameBoundaryGenerator
    input = gen
    old_boundary = '0 1 2 3'
    new_boundary = 'inlet right_wall outlet left_wall'
  []
  [matrix_subdomain]
    type = RenameBlockGenerator
    input = rename
    old_block = 0
    new_block = matrix
  []
  [fracture_sideset]
    type = ParsedGenerateSideset
    input = matrix_subdomain
    combinatorial_geometry = 'x>${fparse -frac_len/2} & x<${fparse frac_len/2}'
    normal = '1 0 0'
    new_sideset_name = fracture_sideset
  []
  [fracture_subdomain]
    type = LowerDBlockFromSidesetGenerator
    input = fracture_sideset
    new_block_id = 1
    new_block_name = fracture
    sidesets = fracture_sideset
  []

  [inlet_node]
    type = ExtraNodesetGenerator
    new_boundary = 'inlet_node'
    coord = '0.0 0.0 0.0'
    input = fracture_subdomain
  []

  # [outlet_node]
  #   type = ExtraNodesetGenerator
  #   new_boundary = 'outlet_node'
  #   coord = '0.0 ${length} 0.0'
  #   input = inlet_node
  # []
  # construct_side_list_from_node_list=true
[]

# [ICs]
# [f_Al]
#   type = ConstantIC
#   variable = f_Al
#   value = ${Al_mass}
#   block = fracture
# []

# [f_HCO3]
#   type = ConstantIC
#   variable = f_HCO3
#   value = ${HCO3_mass}
#   block = fracture
# []
# []

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
    type = PorousFlowDarcyVelocityComponentLowerDimensional
    variable = velocity_x
    component = x
    aperture = ${frac_len}
  []
  [velocity_y]
    type = PorousFlowDarcyVelocityComponentLowerDimensional
    variable = velocity_y
    component = y
    aperture = ${frac_len}
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
    boundary = inlet_node
    value = ${Al_mass}
  []
  [f_HCO3_in]
    type = DirichletBC
    variable = f_HCO3
    boundary = inlet_node
    value = ${HCO3_mass}
  []
  # [f_H2O_in]
  #   type = DirichletBC
  #   variable = porepressure
  #   boundary = inlet_node
  #   value = ${H2O_mass}
  # []

  # [f_Al_node]
  #   type = PorousFlowOutflowBC
  #   variable = f_Al
  #   boundary = outlet_node
  #   include_relperm = false
  # []
  # [f_HCO3_node]
  #   type = PorousFlowOutflowBC
  #   variable = f_HCO3
  #   boundary = outlet_node
  #   include_relperm = false
  # []
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

# [UserObjects]
#   [nnn_uo]
#     type = NearestNodeNumberUO
#     point = '0 ${length} 0'
#     execute_on = 'initial timestep_begin'
#   []
# []

[Postprocessors]
  # [nnn]
  #   type = NearestNodeNumber
  #   nearest_node_number_uo = nnn_uo
  #   execute_on = 'initial timestep_begin'
  # []
  [mass_extracted_Al_in]
    type = NodalExtremeValue
    variable = f_Al
    boundary = inlet
    execute_on = 'initial timestep_end'
  []
  [mass_extracted_HCO3_in]
    type = NodalExtremeValue
    variable = f_HCO3
    boundary = inlet
    execute_on = 'initial timestep_end'
  []
  [mass_extracted_Al_out]
    type = NodalExtremeValue
    variable = f_Al
    boundary = outlet
    execute_on = 'initial timestep_end'
  []
  [mass_extracted_HCO3_out]
    type = NodalExtremeValue
    variable = f_HCO3
    boundary = outlet
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
    permeability = '${domain_permeability} 0 0  0 ${frac_permeability} 0  0 0 ${domain_permeability}'
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
  end_time = 2 #${end_time} #126092.85
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
