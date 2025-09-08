[Mesh]
  [generated_mesh]
    type = AnnularMeshGenerator
    nt = 45
    nr = 10
    rmax = 10.0
    rmin = 0.0
  []
[]

[Variables]
  [T]
    initial_condition = 300.0
    family = LAGRANGE
    order = FIRST
  []
[]

[Kernels]
  [heat_conduction]
    type = HeatConduction
    variable = T
  []
  [heat_source]
    type = HeatSource
    variable = T
    value = 100.0
  []
  [time_derivative]
    type = HeatConductionTimeDerivative
    variable = T
  []
[]

[Materials]
  [heat_conductivity]
    type = HeatConductionMaterial
    thermal_conductivity = 45.0
    specific_heat = 0.5
  []
  [density]
    type = GenericConstantMaterial
    prop_names = 'density'
    prop_values = 8000.0
  []
[]

[BCs]
  [outer_boundary]
    type = DirichletBC
    variable = T
    boundary = rmax
    value = 300.0
  []
[]

[Executioner]
  type = Transient
  end_time = 1e30
  dt = 100.0
[]

[Outputs]
  exodus = true
[]
