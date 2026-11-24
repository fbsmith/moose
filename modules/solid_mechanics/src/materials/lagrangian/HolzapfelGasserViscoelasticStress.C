#include "HolzapfelGasserViscoelasticStress.h"

registerMooseObject("SolidMechanicsApp", HolzapfelGasserViscoelasticStress);

InputParameters
HolzapfelGasserViscoelasticStress::validParams()
{
  InputParameters p = ComputeLagrangianStressPK2::validParams();
  p.addClassDescription("Holzapfel-Gasser (2001) viscoelastic stress, two fiber families, finite strain.");
  MooseEnum iso("MOONEY_RIVLIN OGDEN", "MOONEY_RIVLIN");
  p.addParam<MooseEnum>("isotropic_model", iso, "Matrix model");
  p.addRequiredParam<Real>("bulk_modulus", "Penalty bulk modulus kappa, Psi_vol = kappa/2 (J-1)^2");
  p.addParam<Real>("c1", 0.0, "Mooney-Rivlin c1 (Eq. 28 convention: c1/2 (I1-3))");
  p.addParam<Real>("c2", 0.0, "Mooney-Rivlin c2 (c2/2 (I2-3))");
  p.addParam<std::vector<Real>>("ogden_mu", {}, "Ogden mu_p (Eq. 39)");
  p.addParam<std::vector<Real>>("ogden_alpha", {}, "Ogden alpha_p (Eq. 39)");
  p.addParam<Real>("k1_a", 0.0, "Fiber stress-like parameter, family a0");
  p.addParam<Real>("k2_a", 1.0, "Fiber dimensionless parameter, family a0");
  p.addParam<Real>("k1_b", 0.0, "Fiber stress-like parameter, family b0");
  p.addParam<Real>("k2_b", 1.0, "Fiber dimensionless parameter, family b0");
  p.addParam<RealVectorValue>("fiber_a", RealVectorValue(0, 0, 0), "Reference direction a0 (normalized internally); zero = no fibers)");
  p.addParam<RealVectorValue>("fiber_b", RealVectorValue(0, 0, 0), "Reference direction b0 (normalized internally); zero = no fibers)");
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
    paramError("beta", "beta and tau need 4 entries each (a = 1, 2, 4, 6)");
  if (_ogden && (_mu.empty() || _mu.size() != _alpha.size()))
    paramError("ogden_mu", "Provide equally sized ogden_mu and ogden_alpha");
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

void
HolzapfelGasserViscoelasticStress::initQpStatefulProperties()
{
  ComputeLagrangianStressPK2::initQpStatefulProperties();
  for (unsigned a = 0; a < 4; ++a)
  {
    (*_Q[a])[_qp].zero();
    (*_Siso[a])[_qp].zero();
  }
}

void
HolzapfelGasserViscoelasticStress::evalStress(const RankTwoTensor & C,
                                              const std::array<RankTwoTensor, 4> & H,
                                              const std::array<Real, 4> & x,
                                              RankTwoTensor & S,
                                              std::array<RankTwoTensor, 4> & Siso) const
{
  const RankTwoTensor I = RankTwoTensor::Identity();
  const Real J = std::sqrt(C.det());
  const RankTwoTensor Cinv = C.inverse();
  const Real jm = std::pow(J, -2.0 / 3.0);
  const RankTwoTensor Cb = jm * C; // modified right Cauchy-Green tensor

  std::array<RankTwoTensor, 4> Sh; // fictitious stresses  S^_a = 2 psi_a D_a

  if (_ogden)
  {
    std::vector<Real> ev;
    RankTwoTensor evec;
    Cb.symmetricEigenvaluesEigenvectors(ev, evec); // eigenvectors are columns
    for (unsigned i = 0; i < 3; ++i)
    {
      Real f = 0;
      for (unsigned p = 0; p < _mu.size(); ++p)
        f += _mu[p] * std::pow(ev[i], 0.5 * (_alpha[p] - 2.0)); // mu lam^(alpha-2)
      const RealVectorValue n(evec(0, i), evec(1, i), evec(2, i));
      Sh[0] += f * RankTwoTensor::selfOuterProduct(n);
    }
  }
  else
  {
    const Real I1 = Cb.trace();
    Sh[0] = _c1 * I;
    Sh[1] = _c2 * (I1 * I - Cb);
  }

  // fibers: psi = k1/(2 k2) {exp[k2 (I4-1)^2] - 1}, active only for I4 >= 1
  const Real I4 = _A0.doubleContraction(Cb), I6 = _B0.doubleContraction(Cb);
  if (_A0.L2norm() > 0 && I4 > 1.0)
    Sh[2] = 2.0 * _k1a * (I4 - 1.0) * std::exp(_k2a * (I4 - 1.0) * (I4 - 1.0)) * _A0;
  if (_B0.L2norm() > 0 && I6 > 1.0)
    Sh[3] = 2.0 * _k1b * (I6 - 1.0) * std::exp(_k2b * (I6 - 1.0) * (I6 - 1.0)) * _B0;

  S = _kappa * J * (J - 1.0) * Cinv; // S_vol = J p C^-1
  for (unsigned a = 0; a < 4; ++a)
  {
    Siso[a] = jm * (Sh[a] - (Sh[a].doubleContraction(C) / 3.0) * Cinv); // J^-2/3 P : S^
    S += (1.0 + _beta[a] * x[a]) * Siso[a] + H[a];
  }
}

void
HolzapfelGasserViscoelasticStress::computeQpPK2Stress()
{
  const RankTwoTensor I = RankTwoTensor::Identity();
  const RankTwoTensor C = I + 2.0 * greenLagrangeStrain();

  std::array<RankTwoTensor, 4> H, Siso, Sd;
  std::array<Real, 4> x;
  for (unsigned a = 0; a < 4; ++a)
  {
    x[a] = (_beta[a] > 0 && _tau[a] > 0) ? std::exp(-_dt / (2.0 * _tau[a])) : 0.0;
    if (_beta[a] > 0)
      H[a] = x[a] * (x[a] * (*_Q_old[a])[_qp] - _beta[a] * (*_Siso_old[a])[_qp]);
  }

  RankTwoTensor S;
  evalStress(C, H, x, S, Siso);
  _pk2_stress[_qp] = S;
  for (unsigned a = 0; a < 4; ++a)
  {
    (*_Q[a])[_qp] = H[a] + _beta[a] * x[a] * Siso[a];
    (*_Siso[a])[_qp] = Siso[a];
  }

  // tangent C_ijkl = dS_ij/dE_kl by central differences (history H held fixed)
  const Real eps = 1e-6;
  for (unsigned k = 0; k < 3; ++k)
    for (unsigned l = 0; l <= k; ++l)
    {
      RankTwoTensor dE;
      if (k == l)
        dE(k, k) = eps;
      else
      {
        dE(k, l) = 0.5 * eps;
        dE(l, k) = 0.5 * eps;
      }
      RankTwoTensor Sp, Sm;
      evalStress(C + 2.0 * dE, H, x, Sp, Sd);
      evalStress(C - 2.0 * dE, H, x, Sm, Sd);
      const RankTwoTensor dS = (Sp - Sm) / ((k == l) ? (4.0 * eps) : eps);
      for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j)
        {
          _pk2_jacobian[_qp](i, j, k, l) = dS(i, j);
          _pk2_jacobian[_qp](i, j, l, k) = dS(i, j);
        }
    }
}
