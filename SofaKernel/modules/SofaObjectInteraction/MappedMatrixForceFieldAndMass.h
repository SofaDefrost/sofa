/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 MGH, INRIA, USTL, UJF, CNRS                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_COMPONENT_FORCEFIELD_MappedMatrixForceFieldAndMass_H
#define SOFA_COMPONENT_FORCEFIELD_MappedMatrixForceFieldAndMass_H

#include <sofa/core/behavior/MixedInteractionForceField.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/core/objectmodel/Link.h>
#include <sofa/core/MechanicalParams.h>
#include <SofaBaseLinearSolver/CompressedRowSparseMatrix.h>
#include <SofaBaseLinearSolver/DefaultMultiMatrixAccessor.h>
#include <SofaObjectInteraction/config.h>

#include <sofa/core/topology/BaseMeshTopology.h>


#include <sofa/simulation/MechanicalVisitor.h>
#include <sofa/core/ConstraintParams.h>
#include <sofa/core/MultiVecId.h>
#include <sofa/core/BaseMapping.h>
#include <sofa/defaulttype/BaseMatrix.h>


namespace sofa
{

namespace component
{

namespace interactionforcefield
{


class MechanicalAccumulateJacobian : public simulation::BaseMechanicalVisitor
{
public:
    MechanicalAccumulateJacobian(const core::ConstraintParams* _cparams, core::MultiMatrixDerivId _res)
        : simulation::BaseMechanicalVisitor(_cparams)
        , res(_res)
        , cparams(_cparams)
    {

    }

    virtual void bwdMechanicalMapping(simulation::Node* node, core::BaseMapping* map)
    {
        ctime_t t0 = begin(node, map);
        map->applyJT(cparams, res, res);
        end(node, map, t0);
    }

    /// Return a class name for this visitor
    /// Only used for debugging / profiling purposes
    virtual const char* getClassName() const { return "MechanicalAccumulateJacobian"; }

    virtual bool isThreadSafe() const
    {
        return false;
    }
    // This visitor must go through all mechanical mappings, even if isMechanical flag is disabled
    virtual bool stopAtMechanicalMapping(simulation::Node* /*node*/, core::BaseMapping* /*map*/)
    {
        return false; // !map->isMechanical();
    }

#ifdef SOFA_DUMP_VISITOR_INFO
    void setReadWriteVectors()
    {
    }
#endif

protected:
    core::MultiMatrixDerivId res;
    const sofa::core::ConstraintParams *cparams;
};




using sofa::component::linearsolver::CompressedRowSparseMatrix ;
using sofa::core::behavior::MixedInteractionForceField ;
using sofa::core::behavior::BaseForceField ;
using sofa::core::behavior::BaseMechanicalState ;
using sofa::core::behavior::BaseMass ;
using sofa::core::behavior::MultiMatrixAccessor ;
using sofa::component::linearsolver::DefaultMultiMatrixAccessor ;
using sofa::defaulttype::BaseMatrix ;
using sofa::core::MechanicalParams ;

template<typename TDataTypes1, typename TDataTypes2>
class MappedMatrixForceFieldAndMass : public MixedInteractionForceField<TDataTypes1, TDataTypes2>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(MappedMatrixForceFieldAndMass, TDataTypes1, TDataTypes2),
               SOFA_TEMPLATE2(MixedInteractionForceField, TDataTypes1, TDataTypes2));

    typedef MixedInteractionForceField<TDataTypes1, TDataTypes2> Inherit;
    // Vec3
    typedef TDataTypes1 DataTypes1;
    typedef typename DataTypes1::VecCoord VecCoord1;
    typedef typename DataTypes1::VecDeriv VecDeriv1;
    typedef typename DataTypes1::Coord    Coord1;
    typedef typename DataTypes1::Deriv    Deriv1;
    typedef typename DataTypes1::Real     Real1;
    typedef typename DataTypes1::MatrixDeriv MatrixDeriv1;
    typedef Data<MatrixDeriv1>  DataMatrixDeriv1;
    typedef typename DataTypes1::MatrixDeriv::RowConstIterator MatrixDeriv1RowConstIterator;
    typedef typename DataTypes1::MatrixDeriv::ColConstIterator MatrixDeriv1ColConstIterator;
    static const unsigned int DerivSize1 = Deriv1::total_size;


    // Rigid
    typedef TDataTypes2 DataTypes2;
    typedef typename DataTypes2::VecCoord VecCoord2;
    typedef typename DataTypes2::VecDeriv VecDeriv2;
    typedef typename DataTypes2::Coord    Coord2;
    typedef typename DataTypes2::Deriv    Deriv2;
    typedef typename DataTypes2::Real     Real2;
    typedef typename DataTypes2::MatrixDeriv MatrixDeriv2;
    typedef Data<MatrixDeriv2>  DataMatrixDeriv2;
    typedef typename DataTypes2::MatrixDeriv::RowConstIterator MatrixDeriv2RowConstIterator;
    typedef typename DataTypes2::MatrixDeriv::ColConstIterator MatrixDeriv2ColConstIterator;
    static const unsigned int DerivSize2 = Deriv2::total_size;

    typedef Data<VecCoord1>    DataVecCoord1;
    typedef Data<VecDeriv1>    DataVecDeriv1;
    typedef Data<VecCoord2>    DataVecCoord2;
    typedef Data<VecDeriv2>    DataVecDeriv2;

    typedef sofa::defaulttype::BaseVector::Index  Index;

    typedef typename CompressedRowSparseMatrix<Real1>::Range  Range;


protected:
    SingleLink < MappedMatrixForceFieldAndMass<DataTypes1, DataTypes2>,
                 BaseForceField, BaseLink::FLAG_STOREPATH | BaseLink::FLAG_STRONGLINK > l_mappedForceField;
    SingleLink < MappedMatrixForceFieldAndMass<DataTypes1, DataTypes2>,
                 BaseForceField, BaseLink::FLAG_STOREPATH | BaseLink::FLAG_STRONGLINK > l_mappedForceField2;
    SingleLink < MappedMatrixForceFieldAndMass<DataTypes1, DataTypes2>,
                 BaseMass, BaseLink::FLAG_STOREPATH | BaseLink::FLAG_STRONGLINK > l_mappedMass;

    MappedMatrixForceFieldAndMass() ;

public:

    virtual void init();

    virtual void addForce(const MechanicalParams* mparams,
                          DataVecDeriv1& f1,
                          DataVecDeriv2& f2,
                          const DataVecCoord1& x1,
                          const DataVecCoord2& x2,
                          const DataVecDeriv1& v1,
                          const DataVecDeriv2& v2) ;

    virtual void addDForce(const MechanicalParams* mparams,
                           DataVecDeriv1& df1,
                           DataVecDeriv2& df2,
                           const DataVecDeriv1& dx1,
                           const DataVecDeriv2& dx2) ;

    virtual void addKToMatrix(const MechanicalParams* mparams,
                              const MultiMatrixAccessor* matrix ) ;

    virtual double getPotentialEnergy(const MechanicalParams* mparams,
                                      const DataVecCoord1& x1, const DataVecCoord2& x2) const ;
protected:
    virtual void buildIdentityBlocksInJacobian(core::behavior::BaseMechanicalState* mstate, sofa::core::MatrixDerivId Id);
    virtual void accumulateJacobiansOptimized(const MechanicalParams* mparams);
    virtual void addMassToSystem(const MechanicalParams* mparams, const DefaultMultiMatrixAccessor* KAccessor);
    virtual void addPrecomputedMassToSystem(const MechanicalParams* mparams,const unsigned int mstateSize,const Eigen::SparseMatrix<double> &Jeig, Eigen::SparseMatrix<double>& JtKJeig);
    void accumulateJacobians(const MechanicalParams* mparams);
    virtual void optimizeAndCopyMappingJacobianToEigenFormat1(const typename DataTypes1::MatrixDeriv& J, Eigen::SparseMatrix<double>& Jeig);
    virtual void optimizeAndCopyMappingJacobianToEigenFormat2(const typename DataTypes2::MatrixDeriv& J, Eigen::SparseMatrix<double>& Jeig);

    ////////////////////////// Inherited attributes ////////////////////////////
    /// https://gcc.gnu.org/onlinedocs/gcc/Name-lookup.html
    /// Bring inherited attributes and function in the current lookup context.
    /// otherwise any access to the base::attribute would require
    /// the "this->" approach.
    using MixedInteractionForceField<TDataTypes1, TDataTypes2>::f_printLog ;
    using MixedInteractionForceField<TDataTypes1, TDataTypes2>::mstate1 ;
    using MixedInteractionForceField<TDataTypes1, TDataTypes2>::mstate2 ;
    using MixedInteractionForceField<TDataTypes1, TDataTypes2>::getContext ;
    using MixedInteractionForceField<TDataTypes1, TDataTypes2>::m_componentstate ;
    ////////////////////////////////////////////////////////////////////////////

    sofa::core::behavior::BaseMechanicalState* m_childState ;
};

} // namespace interactionforcefield

} // namespace component

} // namespace sofa


#endif
