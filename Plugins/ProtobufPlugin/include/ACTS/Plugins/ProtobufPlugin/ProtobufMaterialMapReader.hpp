// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ACTS/Material/IMaterialMapReader.hpp"

#include <map>
#include <string>

namespace google {
namespace protobuf {
class MessageLite;

namespace io {
class ZeroCopyInputStream;
}
}
}

namespace Acts {

class BinnedSurfaceMaterial;
class GeometryID;
class MaterialProperties;

class ProtobufMaterialMapReader : IMaterialMapReader {

public:
  struct Config 
  {
    std::string infile;
  };

  ProtobufMaterialMapReader(const Config& cfg) 
    : m_cfg(cfg)
  {
  }

  std::map<GeometryID, BinnedSurfaceMaterial*>
  read() const override;

private:
  Config m_cfg;

  // see
  // https://stackoverflow.com/questions/2340730/are-there-c-equivalents-for-the-protocol-buffers-delimited-i-o-functions-in-ja/22927149
  bool
  readDelimitedFrom(google::protobuf::io::ZeroCopyInputStream* rawInput,
                    google::protobuf::MessageLite*             message) const;

  bool
  checkMaterialProperties(const Acts::MaterialProperties& matProp) const;
};

}
