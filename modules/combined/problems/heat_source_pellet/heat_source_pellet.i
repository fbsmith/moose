[Mesh]
  [fuel_pellet]
    type = AnnularMeshGenerator
    nt = 45
    nr = 10
    rmin = 0.0
    rmax = 0.005 # 5 mm                                                              # this one is solvable
    # rmax = 0.004096 # Dooley = Design of an optimized nuclear fuel pellet          # this one is not solvable, it always crashes for me
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
    value = 5.97e8 # W / m3    taken from AI overview, want to find better source
  []
  [time_derivative]
    type = HeatConductionTimeDerivative
    variable = T
  []
[]

[Materials]
  [fuel_heat_conductivity]
    type = HeatConductionMaterial
    thermal_conductivity = 2.06 # W / m-K    taken from IAEA report at 2000 K for 95% dense UO2
    specific_heat = 373.0 # J / kg-K    taken from IAEA report at 2000 K
  []
  [fuel_density]
    type = GenericConstantMaterial
    prop_names = 'density'
    prop_values = 10320.0 # kg / m3    taken from IAEA report at 2000 K
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
  end_time = 100.0
  dt = 1.0

  # solve_type = 'PJFNK'

  # petsc_options_iname = '-pc_type -pc_hypre_type'
  # petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  exodus = true
[]
