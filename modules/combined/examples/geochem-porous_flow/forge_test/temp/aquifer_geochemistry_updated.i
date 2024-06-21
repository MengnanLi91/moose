
[UserObjects]
  [rate_Albite]
    type = GeochemistryKineticRate
    kinetic_species_name = Albite
    intrinsic_rate_constant = 1E-17
    multiply_by_mass = true
    area_quantity = 10
    activation_energy = 69.8E3
    one_over_T0 = 0.003354
  []
  [rate_Anhydrite]
    type = GeochemistryKineticRate
    kinetic_species_name = Anhydrite
    intrinsic_rate_constant = 1.0E-7
    multiply_by_mass = true
    area_quantity = 10
    activation_energy = 14.3E3
    one_over_T0 = 0.003354
  []
  [rate_Anorthite]
    type = GeochemistryKineticRate
    kinetic_species_name = Anorthite
    intrinsic_rate_constant = 1.0E-13
    multiply_by_mass = true
    area_quantity = 10
    activation_energy = 17.8E3
    one_over_T0 = 0.003354
  []
  [rate_Calcite]
    type = GeochemistryKineticRate
    kinetic_species_name = Calcite
    intrinsic_rate_constant = 1.0E-10
    multiply_by_mass = true
    area_quantity = 10
    activation_energy = 23.5E3
    one_over_T0 = 0.003354
  []
  [rate_Clinochl-7A]
    type = GeochemistryKineticRate
    kinetic_species_name = Clinochl-7A
    intrinsic_rate_constant = 1.0E-17
    multiply_by_mass = true
    area_quantity = 10
    activation_energy = 88.0E3
    one_over_T0 = 0.003354
  []
  [rate_Illite]
    type = GeochemistryKineticRate
    kinetic_species_name = Illite
    intrinsic_rate_constant = 1E-17
    multiply_by_mass = true
    area_quantity = 10
    activation_energy = 29E3
    one_over_T0 = 0.003354
  []
  [rate_K-feldspar]
    type = GeochemistryKineticRate
    kinetic_species_name = K-feldspar
    intrinsic_rate_constant = 1E-17
    multiply_by_mass = true
    area_quantity = 10
    activation_energy = 38E3
    one_over_T0 = 0.003354
  []
  [rate_Quartz]
    type = GeochemistryKineticRate
    kinetic_species_name = Quartz
    intrinsic_rate_constant = 1E-18
    multiply_by_mass = true
    area_quantity = 10
    activation_energy = 90.1E3
    one_over_T0 = 0.003354
  []
  [rate_Phlogopite]
    type = GeochemistryKineticRate
    kinetic_species_name = Phlogopite
    intrinsic_rate_constant = 1E-17
    multiply_by_mass = true
    area_quantity = 10
    activation_energy = 22E3
    one_over_T0 = 0.003354
  []
  [rate_Dolomite]
    type = GeochemistryKineticRate
    kinetic_species_name = Dolomite
    intrinsic_rate_constant = 1E-12
    multiply_by_mass = true
    area_quantity = 10
    activation_energy = 52.2E3
    one_over_T0 = 0.003354
  []
  [rate_Tremolite]
    type = GeochemistryKineticRate
    kinetic_species_name = Tremolite
    intrinsic_rate_constant = 1E-15
    multiply_by_mass = true
    area_quantity = 10
    activation_energy = 94.4E3
    one_over_T0 = 0.003354
  []
  [definition]
    type = GeochemicalModelDefinition
    database_file = '../../../../geochemistry/database/moose_geochemdb.json'
    basis_species = 'H2O H+ Na+ K+ Ca++ Mg++ Al+++ SiO2(aq) Cl- SO4-- HCO3- B(OH)3'
    remove_all_extrapolated_secondary_species = true
    kinetic_minerals = 'Albite Anorthite Calcite Illite K-feldspar Quartz  Phlogopite Anhydrite Clinochl-7A Dolomite Tremolite'
    kinetic_rate_descriptions = 'rate_Albite rate_Anhydrite rate_Anorthite rate_Calcite rate_Clinochl-7A rate_Illite rate_K-feldspar rate_Quartz rate_Phlogopite rate_Dolomite rate_Tremolite'
  []
  [nodal_void_volume_uo]
    type = NodalVoidVolume
    porosity = porosity
    execute_on = 'initial timestep_end' # "initial" means this is evaluated properly for the first timestep
  []
[]

[SpatialReactionSolver]
  model_definition = definition
  geochemistry_reactor_name = reactor
  charge_balance_species = 'Cl-'
  constraint_species = 'H2O                 H+               Na+              K+                Ca++              Mg++           Al+++            SiO2(aq)              Cl-                SO4--               HCO3-             B(OH)3'
  constraint_value = '  1.000              1.280e-07           0.08937         0.005307         6.70E-05          3.45E-08        1.776e-10         0.004823              0.07743             0.006454            0.001328           0.002301         '
  constraint_meaning = 'kg_solvent_water free_concentration free_concentration  free_concentration  free_concentration free_concentration      free_concentration    free_concentration    bulk_composition   free_concentration   free_concentration      free_concentration'
  constraint_unit = '   kg                 molal     molal            molal              molal             molal               molal           molal               moles                molal              molal               molal            '
  initial_temperature = 225
  temperature = temperature
  kinetic_species_name = '        Albite             Anorthite           Calcite       Illite    K-feldspar       Quartz       Phlogopite   Anhydrite Clinochl-7A Dolomite  Tremolite   '
  kinetic_species_initial_value = '996               104.3               0.001           0.2607   467.4       7828         167.9    220.5  36.01     0.5427   110.9   '
  kinetic_species_unit = '         moles              moles              moles         moles       moles          moles          moles        moles     moles     moles       moles  '
  evaluate_kinetic_rates_always = true # otherwise will easily "run out" of dissolving species
  source_species_names = 'H2O H+ Na+ K+ Ca++ Mg++ Al+++ SiO2(aq) Cl- SO4-- HCO3- B(OH)3'
  source_species_rates = 'rate_H2O_per_1l rate_H_per_1l rate_Na_per_1l rate_K_per_1l rate_Ca_per_1l rate_Mg_per_1l rate_Al_per_1l rate_SiO2_per_1l rate_Cl_per_1l rate_SO4_per_1l rate_HCO3_per_1l rate_BOH3_per_1l'
  ramp_max_ionic_strength_initial = 0 # max_ionic_strength in such a simple problem does not need ramping
  execute_console_output_on = ''
  add_aux_molal = false # save some memory and reduce variables in output exodus
  add_aux_mg_per_kg = false # save some memory and reduce variables in output exodus
  add_aux_free_mg = false # save some memory and reduce variables in output exodus
  add_aux_activity = false # save some memory and reduce variables in output exodus
  add_aux_bulk_moles = false # save some memory and reduce variables in output exodus
  adaptive_timestepping = true
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 15
    ny = 10
    xmin = -100
    xmax = 200
    ymin = -100
    ymax = 100
  []
  [injection_node]
    input = gen
    type = ExtraNodesetGenerator
    new_boundary = injection_node
    coord = '0 0 0'
  []
[]

[Executioner]
  type = Transient
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 50
    dt = 500
  []
  end_time = 31536000
[]

[AuxVariables]
  [temperature]
    initial_condition = 220.0
  []
  [porosity]
    initial_condition = 1.8e-4
  []
  [nodal_void_volume]
  []
  [free_cm3_Kfeldspar] # necessary because of the minus sign in K-feldspar which does not parse correctly in the porosity AuxKernel
  []
  [free_cm3_Clinochl7A] # necessary because of the minus sign in Clinochl-7A which does not parse correctly in the porosity AuxKernel
  []
  [pf_rate_H] # change in H mass (kg/s) at each node provided by the porous-flow simulation
  []
  [pf_rate_Na]
  []
  [pf_rate_K]
  []
  [pf_rate_Ca]
  []
  [pf_rate_Mg]
  []
  [pf_rate_SiO2]
  []
  [pf_rate_Al]
  []
  [pf_rate_Fe]
  []
  [pf_rate_Cl]
  []
  [pf_rate_SO4]
  []
  [pf_rate_HCO3]
  []
  [pf_rate_BOH3]
  []
  [pf_rate_H2O] # change in H2O mass (kg/s) at each node provided by the porous-flow simulation
  []
  [rate_H_per_1l]
  []
  [rate_Na_per_1l]
  []
  [rate_K_per_1l]
  []
  [rate_Ca_per_1l]
  []
  [rate_Mg_per_1l]
  []
  [rate_SiO2_per_1l]
  []
  [rate_Al_per_1l]
  []
  [rate_Fe_per_1l]
  []
  [rate_Cl_per_1l]
  []
  [rate_SO4_per_1l]
  []
  [rate_HCO3_per_1l]
  []
  [rate_BOH3_per_1l]
  []
  [rate_H2O_per_1l]
  []
  [transported_H]
  []
  [transported_Na]
  []
  [transported_K]
  []
  [transported_Ca]
  []
  [transported_Mg]
  []
  [transported_SiO2]
  []
  [transported_Al]
  []
  [transported_Fe]
  []
  [transported_Cl]
  []
  [transported_SO4]
  []
  [transported_HCO3]
  []
  [transported_BOH3]
  []
  [transported_H2O]
  []
  [transported_mass]
  []
  [massfrac_H]
  []
  [massfrac_Na]
  []
  [massfrac_K]
  []
  [massfrac_Ca]
  []
  [massfrac_Mg]
  []
  [massfrac_SiO2]
  []
  [massfrac_Al]
  []
  [massfrac_Fe]
  []
  [massfrac_Cl]
  []
  [massfrac_SO4]
  []
  [massfrac_HCO3]
  []
  [massfrac_BOH3]
  []
  [massfrac_H2O]
  []
[]

[AuxKernels]
  [free_cm3_Kfeldspar]
    type = GeochemistryQuantityAux
    variable = free_cm3_Kfeldspar
    species = 'K-feldspar'
    quantity = free_cm3
    execute_on = 'timestep_begin timestep_end'
  []
  [free_cm3_Clinochl7A]
    type = GeochemistryQuantityAux
    variable = free_cm3_Clinochl7A
    species = 'Clinochl-7A'
    quantity = free_cm3
    execute_on = 'timestep_begin timestep_end'
  []
  [porosity_auxk]
    type = ParsedAux
    coupled_variables = 'free_cm3_Albite free_cm3_Anhydrite free_cm3_Anorthite free_cm3_Calcite  free_cm3_Clinochl7A free_cm3_Illite free_cm3_Kfeldspar  free_cm3_Quartz  free_cm3_Phlogopite free_cm3_Dolomite free_cm3_Tremolite  '
    expression = '1000.0 / (1000.0 + free_cm3_Albite + free_cm3_Anhydrite + free_cm3_Anorthite + free_cm3_Calcite  + free_cm3_Clinochl7A + free_cm3_Illite + free_cm3_Kfeldspar  + free_cm3_Quartz  + free_cm3_Phlogopite + free_cm3_Dolomite+ free_cm3_Tremolite)'
    variable = porosity
    execute_on = 'timestep_end'
  []
  [nodal_void_volume_auxk]
    type = NodalVoidVolumeAux
    variable = nodal_void_volume
    nodal_void_volume_uo = nodal_void_volume_uo
    execute_on = 'initial timestep_end' # "initial" to ensure it is properly evaluated for the first timestep
  []
  [rate_H_per_1l_auxk]
    type = ParsedAux
    coupled_variables = 'pf_rate_H nodal_void_volume'
    variable = rate_H_per_1l
    expression = 'pf_rate_H / 1.0079 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [rate_Na_per_1l_auxk]
    type = ParsedAux
    coupled_variables = 'pf_rate_Na nodal_void_volume'
    variable = rate_Na_per_1l
    expression = 'pf_rate_Na / 22.9898 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [rate_K_per_1l_auxk]
    type = ParsedAux
    coupled_variables = 'pf_rate_K nodal_void_volume'
    variable = rate_K_per_1l
    expression = 'pf_rate_K / 39.0983 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [rate_Ca_per_1l_auxk]
    type = ParsedAux
    coupled_variables = 'pf_rate_Ca nodal_void_volume'
    variable = rate_Ca_per_1l
    expression = 'pf_rate_Ca / 40.08 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [rate_Mg_per_1l_auxk]
    type = ParsedAux
    coupled_variables = 'pf_rate_Mg nodal_void_volume'
    variable = rate_Mg_per_1l
    expression = 'pf_rate_Mg / 24.305 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [rate_SiO2_per_1l_auxk]
    type = ParsedAux
    coupled_variables = 'pf_rate_SiO2 nodal_void_volume'
    variable = rate_SiO2_per_1l
    expression = 'pf_rate_SiO2 / 60.0843 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [rate_Al_per_1l_auxk]
    type = ParsedAux
    coupled_variables = 'pf_rate_Al nodal_void_volume'
    variable = rate_Al_per_1l
    expression = 'pf_rate_Al / 26.9815 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [rate_Fe_per_1l_auxk]
    type = ParsedAux
    coupled_variables = 'pf_rate_Fe nodal_void_volume'
    variable = rate_Fe_per_1l
    expression = 'pf_rate_Fe / 55.847 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [rate_Cl_per_1l_auxk]
    type = ParsedAux
    coupled_variables = 'pf_rate_Cl nodal_void_volume'
    variable = rate_Cl_per_1l
    expression = 'pf_rate_Cl / 35.453 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [rate_SO4_per_1l_auxk]
    type = ParsedAux
    coupled_variables = 'pf_rate_SO4 nodal_void_volume'
    variable = rate_SO4_per_1l
    expression = 'pf_rate_SO4 / 96.0576 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [rate_HCO3_per_1l_auxk]
    type = ParsedAux
    coupled_variables = 'pf_rate_HCO3 nodal_void_volume'
    variable = rate_HCO3_per_1l
    expression = 'pf_rate_HCO3 / 61.0171 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [rate_BOH3_per_1l_auxk]
    type = ParsedAux
    coupled_variables = 'pf_rate_BOH3 nodal_void_volume'
    variable = rate_BOH3_per_1l
    expression = 'pf_rate_BOH3 / 61.8329 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [rate_H2O_per_1l_auxk]
    type = ParsedAux
    coupled_variables = 'pf_rate_H2O nodal_void_volume'
    variable = rate_H2O_per_1l
    expression = 'pf_rate_H2O / 18.01801802 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [transported_H_auxk]
    type = GeochemistryQuantityAux
    variable = transported_H
    species = 'H+'
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_begin'
  []
  [transported_Na_auxk]
    type = GeochemistryQuantityAux
    variable = transported_Na
    species = 'Na+'
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_begin'
  []
  [transported_K_auxk]
    type = GeochemistryQuantityAux
    variable = transported_K
    species = 'K+'
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_begin'
  []
  [transported_Ca_auxk]
    type = GeochemistryQuantityAux
    variable = transported_Ca
    species = 'Ca++'
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_begin'
  []
  [transported_Mg_auxk]
    type = GeochemistryQuantityAux
    variable = transported_Mg
    species = 'Mg++'
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_begin'
  []
  [transported_SiO2_auxk]
    type = GeochemistryQuantityAux
    variable = transported_SiO2
    species = 'SiO2(aq)'
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_begin'
  []
  [transported_Al_auxk]
    type = GeochemistryQuantityAux
    variable = transported_Al
    species = 'Al+++'
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_begin'
  []
  [transported_Cl_auxk]
    type = GeochemistryQuantityAux
    variable = transported_Cl
    species = 'Cl-'
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_begin'
  []
  [transported_SO4_auxk]
    type = GeochemistryQuantityAux
    variable = transported_SO4
    species = 'SO4--'
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_begin'
  []
  [transported_HCO3_auxk]
    type = GeochemistryQuantityAux
    variable = transported_HCO3
    species = 'HCO3-'
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_begin'
  []
  [transported_BOH3_auxk]
    type = GeochemistryQuantityAux
    variable = transported_HCO3
    species = 'B(OH)3'
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_begin'
  []
  [transported_H2O_auxk]
    type = GeochemistryQuantityAux
    variable = transported_H2O
    species = 'H2O'
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_begin'
  []
  [transported_mass_auxk]
    type = ParsedAux
    coupled_variables = ' transported_H transported_Na transported_K transported_Ca transported_Mg transported_SiO2 transported_Al transported_Fe transported_Cl transported_SO4 transported_HCO3 transported_BOH3 transported_H2O'
    variable = transported_mass
    expression = 'transported_H * 1.0079 + transported_Cl * 35.453 + transported_SO4 * 96.0576 + transported_HCO3 * 61.0171 + transported_SiO2 * 60.0843 + transported_Al * 26.9815 + transported_Fe * 55.847 + transported_Ca * 40.08 + transported_Mg * 24.305 + transported_K * 39.0983 + transported_Na * 22.9898 + transported_BOH3*61.8329+ transported_H2O * 18.01801802'
    execute_on = 'timestep_end'
  []
  [massfrac_H_auxk]
    type = ParsedAux
    coupled_variables = 'transported_H transported_mass'
    variable = massfrac_H
    expression = 'transported_H * 1.0079 / transported_mass'
    execute_on = 'timestep_end'
  []
  [massfrac_Na_auxk]
    type = ParsedAux
    coupled_variables = 'transported_Na transported_mass'
    variable = massfrac_Na
    expression = 'transported_Na * 22.9898 / transported_mass'
    execute_on = 'timestep_end'
  []
  [massfrac_K_auxk]
    type = ParsedAux
    coupled_variables = 'transported_K transported_mass'
    variable = massfrac_K
    expression = 'transported_K * 39.0983 / transported_mass'
    execute_on = 'timestep_end'
  []
  [massfrac_Ca_auxk]
    type = ParsedAux
    coupled_variables = 'transported_Ca transported_mass'
    variable = massfrac_Ca
    expression = 'transported_Ca * 40.08 / transported_mass'
    execute_on = 'timestep_end'
  []
  [massfrac_Mg_auxk]
    type = ParsedAux
    coupled_variables = 'transported_Mg transported_mass'
    variable = massfrac_Mg
    expression = 'transported_Mg * 24.305 / transported_mass'
    execute_on = 'timestep_end'
  []
  [massfrac_SiO2_auxk]
    type = ParsedAux
    coupled_variables = 'transported_SiO2 transported_mass'
    variable = massfrac_SiO2
    expression = 'transported_SiO2 * 60.0843 / transported_mass'
    execute_on = 'timestep_end'
  []
  [massfrac_Al_auxk]
    type = ParsedAux
    coupled_variables = 'transported_Al transported_mass'
    variable = massfrac_Al
    expression = 'transported_Al * 26.9815 / transported_mass'
    execute_on = 'timestep_end'
  []
  [massfrac_Fe_auxk]
    type = ParsedAux
    coupled_variables = 'transported_Fe transported_mass'
    variable = massfrac_Fe
    expression = 'transported_Fe * 55.847 / transported_mass'
    execute_on = 'timestep_end'
  []
  [massfrac_Cl_auxk]
    type = ParsedAux
    coupled_variables = 'transported_Cl transported_mass'
    variable = massfrac_Cl
    expression = 'transported_Cl * 35.453 / transported_mass'
    execute_on = 'timestep_end'
  []
  [massfrac_SO4_auxk]
    type = ParsedAux
    coupled_variables = 'transported_SO4 transported_mass'
    variable = massfrac_SO4
    expression = 'transported_SO4 * 96.0576 / transported_mass'
    execute_on = 'timestep_end'
  []
  [massfrac_HCO3_auxk]
    type = ParsedAux
    coupled_variables = 'transported_HCO3 transported_mass'
    variable = massfrac_HCO3
    expression = 'transported_HCO3 * 61.0171 / transported_mass'
    execute_on = 'timestep_end'
  []
  [massfrac_BOH3_auxk]
    type = ParsedAux
    coupled_variables = 'transported_BOH3 transported_mass'
    variable = massfrac_BOH3
    expression = 'transported_BOH3 * 61.8329 / transported_mass'
    execute_on = 'timestep_end'
  []
  [massfrac_H2O_auxk]
    type = ParsedAux
    coupled_variables = 'transported_H2O transported_mass'
    variable = massfrac_H2O
    expression = 'transported_H2O * 18.01801802 / transported_mass'
    execute_on = 'timestep_end'
  []
[]

[GlobalParams]
  point = '0 0 0'
  reactor = reactor
[]

[Postprocessors]
  [temperature]
    type = PointValue
    variable = 'solution_temperature'
  []
  [porosity]
    type = PointValue
    variable = porosity
  []
  [solution_temperature]
    type = PointValue
    variable = solution_temperature
  []
  [massfrac_H]
    type = PointValue
    variable = massfrac_H
  []
  [massfrac_Na]
    type = PointValue
    variable = massfrac_Na
  []
  [massfrac_K]
    type = PointValue
    variable = massfrac_K
  []
  [massfrac_Ca]
    type = PointValue
    variable = massfrac_Ca
  []
  [massfrac_Mg]
    type = PointValue
    variable = massfrac_Mg
  []
  [massfrac_SiO2]
    type = PointValue
    variable = massfrac_SiO2
  []
  [massfrac_Al]
    type = PointValue
    variable = massfrac_Al
  []
  [massfrac_Fe]
    type = PointValue
    variable = massfrac_Fe
  []
  [massfrac_Cl]
    type = PointValue
    variable = massfrac_Cl
  []
  [massfrac_SO4]
    type = PointValue
    variable = massfrac_SO4
  []
  [massfrac_HCO3]
    type = PointValue
    variable = massfrac_HCO3
  []
  [massfrac_BOH3]
    type = PointValue
    variable = massfrac_BOH3
  []
  [massfrac_H2O]
    type = PointValue
    variable = massfrac_H2O
  []
  [cm3_Albite]
    type = PointValue
    variable = 'free_cm3_Albite'
  []
  [cm3_Anhydrite]
    type = PointValue
    variable = 'free_cm3_Anhydrite'
  []
  [cm3_Anorthite]
    type = PointValue
    variable = 'free_cm3_Anorthite'
  []
  [cm3_Calcite]
    type = PointValue
    variable = 'free_cm3_Calcite'
  []
  [cm3_Clinochl-7A]
    type = PointValue
    variable = 'free_cm3_Clinochl-7A'
  []
  [cm3_Illite]
    type = PointValue
    variable = 'free_cm3_Illite'
  []
  [cm3_K-feldspar]
    type = PointValue
    variable = 'free_cm3_K-feldspar'
  []
  [cm3_Quartz]
    type = PointValue
    variable = 'free_cm3_Quartz'
  []
  [cm3_Phlogopite]
    type = PointValue
    variable = 'free_cm3_Phlogopite'
  []
  [cm3_Dolomite]
    type = PointValue
    variable = 'free_cm3_Dolomite'
  []
  [cm3_Tremolite]
    type = PointValue
    variable = 'free_cm3_Tremolite'
  []
  [cm3_mineral]
    type = LinearCombinationPostprocessor
    pp_names = 'cm3_Albite cm3_Anhydrite cm3_Anorthite cm3_Calcite cm3_Clinochl-7A cm3_Illite cm3_K-feldspar cm3_Quartz cm3_Phlogopite cm3_Dolomite cm3_Tremolite'
    pp_coefs = '1 1 1 1 1 1 1 1 1 1 1 '
  []
  [pH]
    type = PointValue
    variable = 'pH'
  []
[]

[Outputs]
  [exo]
    type = Exodus
    execute_on = final
  []
  csv = true
[]

