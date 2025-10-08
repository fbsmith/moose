#pragma once

#include "Material.h"
#include "MooseTypes.h"

class FissionGasMaterial : public Material
{
public:
  static InputParameters validParams();

  FissionGasMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;
}
