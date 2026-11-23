#include "HolzapfelGasserViscoelasticStress.h"

registerMooseObject("SolidMechanicsApp", HolzapfelGasserViscoelasticStress);

InputParameters
HolzapfelGasserViscoelasticStress::validParams();
{
  InputParameters p = ComputeLagrangianStressPK2::validParams();
  p.addClassDescription("Holzapfel-Gasser (2001) viscoelastic stress, two fiber families, finite strain.");
  MooseEnum iso("MOONEY_RIVLIN OGDEN", "MOONEY-RIVLIN");
  p.addParam<MooseEnum>("isotropic_model", iso, "Matrix model");
  p.addRequiredParam<Real>("bulk_modulus", "Penalty bulk modulus kappa, Psi_vol = kappa/2 (J-1)^2");
  p.addParam<Real>("c1", 0.0, "Mooney-Rivlin c1 (Eq. 28 convention: c1/2 (I1-3))");
  p.addParam<Real>("c2", 0.0, "Mooney-Rivlin c2 (c2/2 (I2-3))");
  p.addParam<std::vector<Real>>("ogden_mu", {}, "Ogden mu_p (Eq. 39)");
  p.addParam<std::vector<Real>>("ogedn_alpha", {}, "Ogden alpha_p (Eq. 39)");
  p.addParam<Real>("k1_a", 0.0, "Fiber stress-like parameter, family a0");
  p.addParam<Real>("k2_a", 1.0, "Fiber dimensionless parameter, family a0");
  p.addParam<Real>("k1_b", 0.0, "Fiber stress-like parameter, family b0");
  p.addParam<Real>("k2_b", 1.0, "Fiber dimensionless parameter, family b0");
  p.addParam<RealVectorValue>("fiber_a", RealVectorValue(0, 0, 0), "Reference direction a0 (normalized internally_; zero = no fibers)");
  p.addParam<RealVectorValue>("fiber_b", RealVectorValue(0, 0, 0), "Reference direction b0 (normalized internally_; zero = no fibers)");
  p.addParam<std::vector<Real>>("beta", {0, 0, 0, 0}, "Free-energy factors beta_inf for slots (a=1,2,4,6)");
  p.addParam<std::vector<Real>>("tau", {1, 1, 1, 1}, "Relaxation times for slots (a=1,2,4,6)");
  return p;
}

HolzapfelGasserViscoelasticStress::HolzapfelGasserViscoelasticStress(const InputParameters & parameters)
  : ComputeLagrangianStressPK2(parameters),
    _ogden(getParam<MooseEnum>("isotropic_model") == "OGDEN"),
    _kappa(getParam<Real>("bulk_modulus")),
    _c1(getParam<Real>("c1")),
    _c2(getParam<Real>("c2")),
    _k1a(getParam<Real>("k1_a")),
    _k2a(getParam<Real>("k2_a")),
    _k1b(getParam<Real>("k1_b")),
    _k2b(getParam<Real>("k2_b")),
    _mu(getParam<std::vector<Real>>("ogden_mu")),
    _alpha(getParam<std::vector<Real>>("ogden_alpha"))
{
  const auto beta = getParam<std::vector<Real>>("beta");
  const auto tau = getParam<std::vector<Real>>("tau");
  if (beta.size() != 4 || tau.size() != 4)
    paramError("beta", "beta and tau need 4 entries each (a=1,2,4,6)");
  if (_ogden && (_mu.empty() || _mu.size() != _alpha.size))
    paramError("egden_mu", "Provide equally sized ogden_mu and ogden_alpha");
  for (unsigned a = 0; a < 4; ++a)
  {
    _beta[a] = beta[a];
    _tau[a] = tau[a];
  }

  auto setup = [](RealVectorValue v, RankTwoTensor & T)
  {
    if (v.norm() > 0)
    {
      v /= v.norm();
      T = RankTwoTensor::selfOuterProduct(v);
    }
  };
  setup(getParam<RealVectorValue>("fiber_a"), _A0);
  setup(getParam<RealVectorValue>("fiber_b"), _B0);

  for (unsigned a = 0; a < 4; ++a)
  {
    const auto qn = "hgo_Q" + std::to_string(a), sn = "hgo_Siso" + std::to_string(a);
    _Q.push_back(&declareProperty<RankTwoTensor>(qn));
    _Siso.push_back(&declareProperty<RankTwoTensor>(sn));
    _Q_old.push_back(&getMaterialPropertyOld<RankTwoTensor>(qn));
    _Siso_old.push_back(&getMaterialPropertyOld<RankTwoTensor>(sn));
  }
}
