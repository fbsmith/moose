//* Holzapfel & Gasser (2001) CMAME 190:4379-4403 -- finite-strain viscoelastic fiber-reinforced model
#pragma once

#include "ComputeLagrangianStressPK2.h"
#include <array>

/**
 * Second Piola-Kirchhoff stress for a nearly-incompressible matrix reinforced by two fiber
 * families, with one Maxwell-type relaxation process per contribution (Box 1 of the paper).
 *
 * Slots (index -> paper's a):  0 -> a=1 (I1),  1 -> a=2 (I2),  2 -> a=4 (fiber family a0),  3 -> a=6 (fiber family b0)
 * With isotropic_model = OGDEN the whole Ogden isotropic stress sits in slot 0 (relaxes with beta[0], tau[0]).
 * The tangent dS/dE is obtained by central differences of the algorithmic stress (not the closed form of Appendix A).
 */
class HolzapfelGasserViscoelasticStress : public ComputeLagrangianStressPK2
{
public:
  static InputParameters validParams();
  HolzapfelGasserViscoelasticStress(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpPK2Stress() override;

  const RankTwoTensor & greenLagrangeStrain() const { return _E[_qp]; }

  void evalStress(const RankTwoTensor & C,
                  const std::array<RankTwoTensor, 4> & H,
                  const std::array<Real, 4> & x,
                  RankTwoTensor & S,
                  std::array<RankTwoTensor, 4> & Siso) const;

  const bool _ogden;
  const Real _kappa, _c1, _c2, _k1a, _k2a, _k1b, _k2b;
  const std::vector<Real> _mu, _alpha;
  std::array<Real, 4> _beta, _tau;
  RankTwoTensor _A0, _B0;

  std::vector<MaterialProperty<RankTwoTensor> *> _Q, _Siso;
  std::vector<const MaterialProperty<RankTwoTensor> *> _Q_old, _Siso_old;
};
