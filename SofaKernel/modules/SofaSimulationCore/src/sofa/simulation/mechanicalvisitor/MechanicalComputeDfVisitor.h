/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#pragma once

#include <sofa/simulation/MechanicalVisitor.h>

namespace sofa::simulation::mechanicalvisitor
{

/** Compute the variation of force corresponding to a variation of position.
This action is typically called after a MechanicalPropagateDxVisitor.
*/
class SOFA_SIMULATION_CORE_API MechanicalComputeDfVisitor : public MechanicalVisitor
{
public:
    sofa::core::MultiVecDerivId res;
    bool accumulate; ///< Accumulate everything back to the DOFs through the mappings
    MechanicalComputeDfVisitor(const sofa::core::MechanicalParams* mparams, sofa::core::MultiVecDerivId res)
            : MechanicalVisitor(mparams) , res(res), accumulate(true)
    {
#ifdef SOFA_DUMP_VISITOR_INFO
        setReadWriteVectors();
#endif
    }
    MechanicalComputeDfVisitor(const sofa::core::MechanicalParams* mparams, sofa::core::MultiVecDerivId res, bool accumulate)
            : MechanicalVisitor(mparams) , res(res), accumulate(accumulate)
    {
#ifdef SOFA_DUMP_VISITOR_INFO
        setReadWriteVectors();
#endif
    }
    Result fwdMechanicalState(simulation::Node* /*node*/,sofa::core::behavior::BaseMechanicalState* mm) override;
    Result fwdMappedMechanicalState(simulation::Node* /*node*/,sofa::core::behavior::BaseMechanicalState* mm) override;
    Result fwdForceField(simulation::Node* /*node*/,sofa::core::behavior::BaseForceField* ff) override;
    void bwdMechanicalMapping(simulation::Node* /*node*/, sofa::core::BaseMapping* map) override;
    void bwdMechanicalState(simulation::Node* /*node*/,sofa::core::behavior::BaseMechanicalState* mm) override;

    /// Return a class name for this visitor
    /// Only used for debugging / profiling purposes
    const char* getClassName() const override {return "MechanicalComputeDfVisitor";}
    virtual std::string getInfos() const override
    {
        std::string name="["+res.getName()+"]";
        if (accumulate) name+= " Accumulating";
        else            name+= " Not Accumulating";
        return name;
    }

    /// Specify whether this action can be parallelized.
    bool isThreadSafe() const override
    {
        return true;
    }
#ifdef SOFA_DUMP_VISITOR_INFO
    void setReadWriteVectors() override
    {
        addWriteVector(res);
    }
#endif
};

}