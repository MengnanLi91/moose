[Mesh]
  [cmg]
    type = CartesianMeshGenerator
    dx = 1
    dim = 1
  []
[]

[Times]
  [external_input]
    type = ControllableInputTimes
    next_time = "0.1 0.2"
    execute_on = 'initial timestep_begin FINAL'
  []
[]

[Problem]
  solve = false
[]

[Functions]
  [dts]
    type = ParsedFunction
    expression = 0.1*exp(t*0.5)
  []
[]


[Executioner]
  type = Transient
  num_steps = 20
 [TimeSteppers]
    [external_time]
      type = TimeSequenceFromTimes
      times = external_input
    []

    [ConstDT1]
      type = FunctionDT
      function = dts
      min_dt = 0.1
    []
 []

[]

[Controls]
  [web_server]
    type = WebServerControl
    execute_on = 'INITIAL TIMESTEP_BEGIN FINAL'
  []
[]


[Outputs]
  [out]
    type = JSON
    execute_system_information_on = none
  []
[]
