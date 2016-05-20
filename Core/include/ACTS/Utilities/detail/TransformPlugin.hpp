// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// GeoPrimitivesHelpers.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef EVENTPRIMITIVES_AMGTRANSFORMPLUGIN_H
#define EVENTPRIMITIVES_AMGTRANSFORMPLUGIN_H


    inline explicit Transform(const Vector3d& rotationMatrixCol0, const Vector3d& rotationMatrixCol1, const Vector3d& rotationMatrixCol2){
        check_template_params();
        m_matrix.block(0,0,3,1) = rotationMatrixCol0;
        m_matrix.block(0,1,3,1) = rotationMatrixCol1;
        m_matrix.block(0,2,3,1) = rotationMatrixCol2;
        if (int(Mode)==Affine)
          makeAffine();
    }

    inline explicit Transform(const Vector3d& translation){
        check_template_params();
        m_matrix.block(0,3,3,1) = translation;
        m_matrix.block(0,0,3,3).setIdentity();
        if (int(Mode)==Affine)
          makeAffine();
    }

    inline explicit Transform(const Matrix<double, 3,3>& rotation, const Vector3d& translation){
        check_template_params();
        m_matrix.block(0,0,3,3) = rotation;
        m_matrix.block(0,3,3,1) = translation;
        if (int(Mode)==Affine)
          makeAffine();
    }
    inline explicit Transform(const Matrix<double, 3,3>& rotation, const TranslationType& translation){
        check_template_params();
        m_matrix.block(0,0,3,3) = rotation;
        m_matrix.block(0,3,3,1) = translation.vector();
        if (int(Mode)==Affine)
          makeAffine();
    }

    inline explicit Transform(const Vector3d& rotationMatrixCol0, const Vector3d& rotationMatrixCol1, const Vector3d& rotationMatrixCol2, const Vector3d& translation){
        check_template_params();
        m_matrix.block(0,0,3,1) = rotationMatrixCol0;
        m_matrix.block(0,1,3,1) = rotationMatrixCol1;
        m_matrix.block(0,2,3,1) = rotationMatrixCol2;
        m_matrix.block(0,3,3,1) = translation;
        if (int(Mode)==Affine)
          makeAffine();
    }

#endif