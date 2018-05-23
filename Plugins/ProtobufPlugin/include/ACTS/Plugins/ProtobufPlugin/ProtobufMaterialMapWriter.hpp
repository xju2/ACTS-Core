// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ACTS/Material/IMaterialMapWriter.hpp"

#include <string>
#include <map>

namespace google {
namespace protobuf {
class MessageLite;

namespace io {
class ZeroCopyOutputStream;
}
}
}

namespace Acts {

class SurfaceMaterial;
class GeometryID;

class ProtobufMaterialMapWriter : IMaterialMapWriter {

public:

  struct Config 
  {
    std::string outfile;
  };

  ProtobufMaterialMapWriter(const Config& cfg) 
    : m_cfg(cfg)
  {
  }

  void
  write(const std::map<GeometryID, SurfaceMaterial *> &surfaceMaterialMap) const
      override;

private:
  
  Config m_cfg;

  // see
  // https://stackoverflow.com/questions/2340730/are-there-c-equivalents-for-the-protocol-buffers-delimited-i-o-functions-in-ja/22927149
  bool
  writeDelimitedTo(const google::protobuf::MessageLite &message,
                   google::protobuf::io::ZeroCopyOutputStream *rawOutput) const;
};

}
