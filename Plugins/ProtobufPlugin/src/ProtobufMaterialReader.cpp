// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "ACTS/Plugins/ProtobufPlugin/ProtobufMaterialMapReader.hpp"
#include "ACTS/Utilities/GeometryID.hpp"

#include "ACTS/Material/SurfaceMaterial.hpp"
#include "ACTS/Material/BinnedSurfaceMaterial.hpp"
#include "ACTS/Material/HomogeneousSurfaceMaterial.hpp"
#include "ACTS/Material/MaterialProperties.hpp"
#include "ACTS/Utilities/ThrowAssert.hpp"

// this is a private include
#include "gen/MaterialMap.pb.h"

// google/protobuf
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/message_lite.h>

#include <vector>
#include <cmath>


bool
Acts::ProtobufMaterialMapReader::readDelimitedFrom(
    google::protobuf::io::ZeroCopyInputStream* rawInput,
    google::protobuf::MessageLite*             message) const
{
  // We create a new coded stream for each message.  Don't worry, this is fast,
  // and it makes sure the 64MB total size limit is imposed per-message rather
  // than on the whole stream.  (See the CodedInputStream interface for more
  // info on this limit.)
  google::protobuf::io::CodedInputStream input(rawInput);

  // Read the size.
  uint32_t size;
  if (!input.ReadVarint32(&size)) return false;

  // Tell the stream not to read beyond that size.
  google::protobuf::io::CodedInputStream::Limit limit = input.PushLimit(size);

  // Parse the message.
  if (!message->MergeFromCodedStream(&input)) return false;
  if (!input.ConsumedEntireMessage()) return false;

  // Release the limit.
  input.PopLimit(limit);

  return true;
}

bool
Acts::ProtobufMaterialMapReader::checkMaterialProperties(
    const Acts::MaterialProperties& matProp) const 
{
  auto valid = [](const auto& value) {
    return !std::isinf(value) && !std::isnan(value);
  };
  const Material& mat = matProp.material();
  return valid(mat.X0()) && valid(mat.L0()) && valid(mat.A()) && valid(mat.Z())
      && valid(mat.rho()) && valid(matProp.thickness());
}

std::map<Acts::GeometryID, Acts::BinnedSurfaceMaterial*>
Acts::ProtobufMaterialMapReader::read() const 
{

  // return value
  std::map<Acts::GeometryID, Acts::BinnedSurfaceMaterial*> surfaceMaterialMap;

  auto infd = fopen(m_cfg.infile.c_str(), "r");
  if (infd == nullptr) {
    throw std::runtime_error("Material Map Input file does not exist.");
  }
  auto instream
      = std::make_unique<google::protobuf::io::FileInputStream>(fileno(infd));
  
  Acts::protobuf::MaterialMap matMapMsg;

  while (readDelimitedFrom(instream.get(), &matMapMsg)) {
    //std::cout << "read material map" << std::endl;

    //size_t nBins = rows * cols;
    //const Acts::protobuf::MaterialMap::Dimension& l0 = matMapMsg.l0();
    //const Acts::protobuf::MaterialMap::Dimension& l1 = matMapMsg.l1();
    
    size_t rows  = matMapMsg.rows();
    size_t cols  = matMapMsg.cols();
    
    int nBins = rows * cols;

    throw_assert(nBins == matMapMsg.bins_size(), "Mismatch between given dimensions and stored bins");

    //std::cout << "l0: " << rows << ", " << l0.min() << " - " << l0.max() << std::endl;
    //std::cout << "l1: " << cols << ", " << l1.min() << " - " << l1.max() << std::endl;
    //std::cout << "geo_id: " << matMapMsg.geo_id() << std::endl;
    //std::cout << "sen_id: " << matMapMsg.sen_id() << std::endl;

    //std::cout << "nbins: " << matMapMsg.bins_size() << std::endl;

    GeometryID geo_id(matMapMsg.geo_id());

    // binned material
    std::vector<MaterialProperties*> matPropVector(rows, nullptr);
    std::vector<std::vector<MaterialProperties*>> matPropMatrix(cols, matPropVector);
    
    for (int i = 0; i < matMapMsg.bins_size(); i++) {

      int col = i / rows;
      int row = i % rows;

      const Acts::protobuf::MaterialMap::MaterialProperties& matPropMsg = matMapMsg.bins(i);
      if (matPropMsg.valid()) {
        MaterialProperties* matProp = new MaterialProperties(matPropMsg.x0(),
                                                            matPropMsg.l0(),
                                                            matPropMsg.a(),
                                                            matPropMsg.z(),
                                                            matPropMsg.rho(),
                                                            matPropMsg.thickness());
        if (!checkMaterialProperties(*matProp)) {
          throw std::domain_error("Got invalid material properties when reading");
        }
        matPropMatrix[col][row] = matProp;
      }
      else {
        matPropMatrix[col][row] = nullptr;
      }

    }

    // at this stage we have insufficient information to construct the
    // bin utility. construct dummy binutility here.
    // this needs to be reset later!
    Acts::BinUtility binUtility(
        rows, 0, 1, Acts::open, Acts::binX);
    binUtility += Acts::BinUtility(
        cols, 0, 1, Acts::open, Acts::binY);

    surfaceMaterialMap[geo_id] = new BinnedSurfaceMaterial(binUtility, matPropMatrix);


    matMapMsg.Clear();
  }

  return surfaceMaterialMap;
}
