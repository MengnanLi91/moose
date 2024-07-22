# Hibiki (2001)
# flow regime
# void fraction          0.045    0.085     0.15
# mixture velocity(m/s)  0.61     0.63      0.77
# particle diameter(m)   0.0025   0.0026    0.003

mu = 1.002e-3
rho = 998.19
mu_d = 1.825e-5
rho_d = 1.204
pipe_diameter = 0.0508 # m
length = 2.7178 # m
U = 0.77 # m/s
dp = 0.003
inlet_phase_2 =  0.15
g = 9.81
advected_interp_method = 'upwind'
velocity_interp_method = 'rc'

[GlobalParams]
  rhie_chow_user_object = 'rc'
  density_interp_method = 'average'
  mu_interp_method = 'average'
[]

[UserObjects]
  [rc]
    type = INSFVRhieChowInterpolator
    u = vel_x
    v = vel_y
    pressure = pressure
  []
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = '${length}'
    ymin = '${fparse -pipe_diameter / 2}'
    ymax = '${fparse pipe_diameter / 2}'
    nx = 400
    ny = 30
  []
 # file = validation_2010_cp/LATEST
[]

# [Problem]
#   restart_file_base = validation_2010_cp/LATEST
# []

[Variables]
  [vel_x]
    type = INSFVVelocityVariable
  []
  [vel_y]
    type = INSFVVelocityVariable
  []
  [pressure]
    type = INSFVPressureVariable
  []
  [phase_2]
    type = INSFVScalarFieldVariable
    initial_condition = ${inlet_phase_2}
  []
[]

[FVKernels]
 # inactive = 'v_drift u_drift'
  [mass]
    type = INSFVMassAdvection
    variable = pressure
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = 'rho_mixture'
  []

  [u_time]
    type = INSFVMomentumTimeDerivative
    variable = vel_x
    rho = 'rho_mixture'
    momentum_component = 'x'
  []
  [u_advection]
    type = INSFVMomentumAdvection
    variable = vel_x
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = 'rho_mixture'
    momentum_component = 'x'
  []
  [u_drift]
    type = WCNSFV2PMomentumDriftFlux
    variable = vel_x
    rho_d = ${rho_d}
    rho = ${rho}
    fd = 'phase_2'
    u_slip = 'vel_slip_x'
    v_slip = 'vel_slip_y'
    momentum_component = 'x'
  []
  [u_viscosity]
    type = INSFVMomentumDiffusion
    variable = vel_x
    mu = 'mu_mixture'
    limit_interpolation = true
    momentum_component = 'x'
  []
  [u_pressure]
    type = INSFVMomentumPressure
    variable = vel_x
    momentum_component = 'x'
    pressure = pressure
  []
  [u_friction]
    type = PINSFVMomentumFriction
    Darcy_name = Darcy_coefficient_vec
    is_porous_medium = false
    momentum_component = x
    mu = mu_mixture
    rho = rho_mixture
    variable = vel_x
  []

  [u_buoyant]
    type = INSFVMomentumGravity
    variable = vel_x
    rho = 'rho_mixture'
    momentum_component = 'x'
    gravity = '${fparse -g} 0 0'
  []

  [v_time]
    type = INSFVMomentumTimeDerivative
    variable = vel_y
    rho = 'rho_mixture'
    momentum_component = 'y'
  []
  [v_advection]
    type = INSFVMomentumAdvection
    variable = vel_y
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = 'rho_mixture'
    momentum_component = 'y'
  []
  [v_drift]
    type = WCNSFV2PMomentumDriftFlux
    variable = vel_y
    rho_d = ${rho_d}
    rho = ${rho}
    fd = 'phase_2'
    u_slip = 'vel_slip_x'
    v_slip = 'vel_slip_y'
    momentum_component = 'y'
  []
  [v_viscosity]
    type = INSFVMomentumDiffusion
    variable = vel_y
    mu = 'mu_mixture'
    limit_interpolation = true
    momentum_component = 'y'
  []
  [v_pressure]
    type = INSFVMomentumPressure
    variable = vel_y
    momentum_component = 'y'
    pressure = pressure
  []
  [v_friction]
    type = PINSFVMomentumFriction
    Darcy_name = Darcy_coefficient_vec
    is_porous_medium = false
    momentum_component = y
    mu = mu_mixture
    rho = rho_mixture
    variable = vel_y
  []
  # [v_buoyant]
  #   type = INSFVMomentumGravity
  #   variable = vel_y
  #   rho = 'rho_mixture'
  #   momentum_component = 'y'
  #   gravity = '${fparse -g} 0 0'
  # []
  [phase_2_time]
    type = FVFunctorTimeKernel
    variable = phase_2
    functor = phase_2
  []
  [phase_2_advection]
    type = INSFVScalarFieldAdvection
    variable = phase_2
    u_slip = 'vel_x'
    v_slip = 'vel_y'
    velocity_interp_method = ${velocity_interp_method}
    advected_interp_method = 'upwind'
  []
  [phase_2_src]
    type = NSFVMixturePhaseInterface
    variable = phase_2
    phase_coupled = phase_1
    alpha = 0 # volumetric exchange coefficient
  []
[]

[FVBCs]
  [inlet-u]
    type = INSFVInletVelocityBC
    boundary = 'left'
    variable = vel_x
    functor = '${U}'
  []
  [inlet-v]
    type = INSFVInletVelocityBC
    boundary = 'left'
    variable = vel_y
    functor = '0'
  []
  [walls-u]
    type = INSFVNoSlipWallBC
    boundary = 'top bottom'
    variable = vel_x
    function = 0
  []
  [walls-v]
    type = INSFVNoSlipWallBC
    boundary = 'top bottom'
    variable = vel_y
    function = 0
  []
  [outlet_p]
    type = INSFVOutletPressureBC
    boundary = 'right'
    variable = pressure
    function = '0'
  []
  [inlet_phase_2]
    type = FVDirichletBC
    boundary = 'left'
    variable = phase_2
    value = ${inlet_phase_2}
  []
[]

[AuxVariables]
  [drag_coefficient]
    type = MooseVariableFVReal
  []
  [rho_mixture_var]
    type = MooseVariableFVReal
  []
  [mu_mixture_var]
    type = MooseVariableFVReal
  []
  [vel_slip_x_var]
    type = MooseVariableFVReal
  []
  [vel_slip_y_var]
    type = MooseVariableFVReal
  []
  [slip_pressure_term]
    type = MooseVariableFVReal
  []
  [vg_x]
    type = MooseVariableFVReal
  []
  [vg_y]
    type = MooseVariableFVReal
  []
[]

[AuxKernels]
  [populate_cd]
    type = FunctorAux
    variable = drag_coefficient
    functor = 'Darcy_coefficient'
  []
  [populate_rho_mixture_var]
    type = FunctorAux
    variable = rho_mixture_var
    functor = 'rho_mixture'
  []
  [populate_mu_mixture_var]
    type = FunctorAux
    variable = mu_mixture_var
    functor = 'mu_mixture'
  []
  [populate_vx_slip_var]
    type = FunctorAux
    variable = vel_slip_x_var
    functor = 'vel_slip_x'
  []
  [populate_vy_slip_var]
    type = FunctorAux
    variable = vel_slip_y_var
    functor = 'vel_slip_y'
  []
  [vg_x]
    type = ParsedAux
    variable = vg_x
    coupled_variables = 'vel_x vel_slip_x_var rho_mixture_var'
    expression = 'vel_x + ${rho}*vel_slip_x_var/rho_mixture_var'
  []
  [vg_y]
    type = ParsedAux
    variable = vg_y
    coupled_variables = 'vel_y vel_slip_y_var rho_mixture_var'
    expression = 'vel_y + ${rho}*vel_slip_y_var/rho_mixture_var'
  []
[]

[FunctorMaterials]
  [populate_u_slip]
    type = WCNSFV2PSlipVelocityFunctorMaterial
    slip_velocity_name = 'vel_slip_x'
    momentum_component = 'x'
    u = 'vel_x'
    v = 'vel_y'
    rho = ${rho}
    mu = 'mu_mixture'
    rho_d = ${rho_d}
    particle_diameter = ${dp}
    linear_coef_name = 'Darcy_coefficient'
    gravity = '${fparse -g} 0 0'
  []
  [populate_v_slip]
    type = WCNSFV2PSlipVelocityFunctorMaterial
    slip_velocity_name = 'vel_slip_y'
    momentum_component = 'y'
    u = 'vel_x'
    v = 'vel_y'
    rho = ${rho}
    mu = 'mu_mixture'
    rho_d = ${rho_d}
    particle_diameter = ${dp}
    linear_coef_name = 'Darcy_coefficient'
    gravity = '${fparse -g} 0 0'
  []
  [compute_phase_1]
    type = ADParsedFunctorMaterial
    property_name = phase_1
    functor_names = 'phase_2'
    expression = '1 - phase_2'
  []
  [CD]
    type = NSFVDispersePhaseDragFunctorMaterial
    rho = 'rho_mixture'
    mu = mu_mixture
    u = 'vel_x'
    v = 'vel_y'
    particle_diameter = ${dp}
  []
  [mixing_material]
    type = NSFVMixtureFunctorMaterial
    phase_2_names = '${rho} ${mu}'
    phase_1_names = '${rho_d} ${mu_d}'
    prop_names = 'rho_mixture mu_mixture'
    phase_1_fraction = 'phase_2'
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  nl_rel_tol = 1e-6
  end_time = 100
  automatic_scaling = true
  #line_search = 'basic'
  [TimeSteppers]
    [iter_dt]
      type = IterationAdaptiveDT
      optimal_iterations = 10
      iteration_window = 5
      growth_factor = 2.0
      cutback_factor = 0.25
      dt = 0.0001
    []
    [const_dt]
      type = ConstantDT
      dt = 5
    []
  []
  steady_state_detection = true
  steady_state_tolerance = 2e-4
  #steady_state_start_time = 9
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_shift_type'
    petsc_options_value = 'lu       NONZERO'
  []
[]

[Outputs]
  exodus = true
  file_base = 'Hibiki_2001'
  # [my_checkpoint]
  #   type = Checkpoint
  #   num_files = 2
  #   interval = 5
  # []
  [CSV]
    type = CSV
    execute_on = 'TIMESTEP_END'
  []
[]

[Postprocessors]

  [rho_outlet]
    type = SideAverageValue
    boundary = 'right'
    variable = 'rho_mixture_var'
  []
  [vslip_x]
    type = SideAverageValue
    boundary = 'right'
    variable = 'vel_slip_x_var'
  []
  [vslip_y]
    type = SideAverageValue
    boundary = 'right'
    variable = 'vel_slip_y_var'
  []
  [void_fraction]
    type = SideAverageValue
    boundary = 'left'
    variable = 'phase_2'
  []
  [vslip_value]
    type = ParsedPostprocessor
    expression = 'sqrt(vslip_x*vslip_x + vslip_y*vslip_y)*vslip_x/abs(vslip_x)'
    pp_names = 'vslip_x vslip_y'
  []
  [vg_x]
    type = SideAverageValue
    boundary = 'right'
    variable = 'vg_x'
  []
  [vg_y]
    type = SideAverageValue
    boundary = 'right'
    variable = 'vg_y'
  []
  [Re_g]
    type = ParsedPostprocessor
    expression = '${rho_d} * ${fparse pipe_diameter} * vg_x* void_fraction/${mu_d}'
    pp_names = 'vg_x void_fraction'
  []

[]
