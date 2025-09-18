[GlobalParams]
  displacements = 'disp_r disp_theta'
[]

[Mesh]
  type = AnnularMesh
  nt = 45
  nr = 10
  rmax = 0.005
  rmin = 0.0
[]

[Physics/SolidMechanics/QuasiStatic]
  [all]
    add_variables = true
    displacements = 'disp_r disp_theta'
    generate_output = 'vonmises_stress'
    material_output_order = SECOND
    out_of_plane_direction = z
  []
[]

[Materials]
  [elasticity]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1e9
    poissons_ratio = 0.3
  []
  [stress]
    type = ComputeLinearElasticStress
  []
[]

[BCs]
  [r_pin]
    type = DirichletBC
    boundary = rmin
    value = 0.0
    variable = disp_r
  []
  [theta_pin]
    type = DirichletBC
    boundary = rmin
    value = 0.0
    variable = disp_theta
  []
  [Pressure]
    [gap]
        boundary = rmax
        function = 1e7*t
    []
  []
[]

[Executioner]
  type = Transient
  end_time = 5
  dt = 1
[]

[Outputs]
  exodus = true
[]
