# Benchmarks for multicomponent reactive transport across a cement/clay interface
# Geometry


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

  [clayey]
    type = RenameBlockGenerator
    input = concrete
    old_block = '4 5 6 7 8 9'
    new_block = '2 2 2 2 2 2'
  []
[]

[Problem]
  type = FEProblem
  coord_type = 'RZ'
  rz_coord_axis = Y
[]

[GlobalParams]
  PorousFlowDictator = 'dictator'
  gravity = '0 0 0'
[]

[AuxVariables]

[]

[AuxKernels]


[]

[Variables]
  [H20]
    initial_condition = 12e6
  []
[]

[Kernels]
  [mass0]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = H20
  []
  [flux0]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    variable = H20
  []
[]

[UserObjects]

[]

[FluidProperties]


[]

[Materials]
  [temperature]
    type = PorousFlowTemperature
    temperature = '45'
  []
  [brineco2]
    type = PorousFlowFluidState
    gas_porepressure = 'pgas'
    z = 'zi'
    temperature_unit = Celsius
    xnacl = 'xnacl'
    capillary_pressure = pc
    fluid_state = fs
  []
  [porosity]
    type = PorousFlowPorosityConst
    porosity = '0.12'
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1e-13 0 0 0 1e-13 0 0 0 1e-13'
  []
  [relperm_water]
    type = PorousFlowRelativePermeabilityVG
    m = 0.457
    phase = 0
    s_res = 0.3
    sum_s_res = 0.35
  []
[]

[BCs]
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type'
    petsc_options_value = 'gmres bjacobi lu NONZERO'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  end_time = 8.64e8
  nl_max_its = 25
  l_max_its = 100
  dtmax = 5e6
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 100
  []
[]

[Postprocessors]
[]

[Outputs]
[]
