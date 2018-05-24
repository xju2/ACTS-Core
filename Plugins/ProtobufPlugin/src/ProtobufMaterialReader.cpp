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

#include <vector>

bool
Acts::ProtobufMaterialReder::readDelimitedFrom(google::protobuf::io::ZeroCopyInputStream* rawInput,
                  google::protobuf::MessageLite*             message)
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

std::map<GeometryID, SurfaceMaterial*>
read() const 
{
  auto infd = fopen(m_cfg.infile.c_str(), "r");
  auto instream
      = std::make_unique<google::protobuf::io::FileInputStream>(fileno(infd));
  
  Acts::protobuf::MaterialMap matMapMsg;

  while (readDelimitedFrom(instream.get(), &matMapMsg)) {
    std::cout << "read material map" << std::endl;

    size_t rows  = matMapMsg.rows();
    size_t cols  = matMapMsg.cols();
    size_t nBins = rows * cols;

    std::cout << "nrows: " << matMapMsg.rows() << std::endl;
    std::cout << "ncols: " << matMapMsg.cols() << std::endl;
    std::cout << "geo_id: " << matMapMsg.geo_id() << std::endl;
    std::cout << "sen_id: " << matMapMsg.sen_id() << std::endl;

    std::cout << "nbins: " << matMapMsg.bins_size() << std::endl;

    GeometryID geo_id(matMapMsg.geo_id());

    std::vector<std::vector<MaterialProperties*>> matPropMatrix;


    for (int i = 0; i < matMapMsg.bins_size(); i++) {

      int row = i / rows;
      int col = i % rows;

      const Acts::protobuf::MaterialMap::MaterialProperties& bin = matMapMsg.bins(i);
    }

    matMap.Clear();
  }
}
