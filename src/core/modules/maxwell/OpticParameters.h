/*
 * ExcitonLayer.h
 *
 *  Created on: Sep 26, 2011
 *      Author: paveryan
 */

#ifndef OPTIC_PARAMETERS_H_
#define OPTIC_PARAMETERS_H_

class OpticParameters
{
  public:
    double epsilon;
    double mu;
    double sPML;

    OpticParameters(MaxwellEquations* maxwell, const Elem* elem) {
      ID subdomain = elem->subdomain_id();
      const Material* material = maxwell->get_environment().get_device().get_material(subdomain);

      epsilon = 1.0;
      mu = 1.0;
      sPML = -1.0;

      if (material != NULL) {
        Database dbA = material->get_database();
        dbA.set_section("permittivity");
        epsilon = dbA.get("permittivity", epsilon);
        mu = dbA.get("permeability", mu);

        dbA.set_section("pml");
        sPML = dbA.get("sPML", sPML);
      }
    }
};

#endif /* EXCITONLAYER_H_ */
