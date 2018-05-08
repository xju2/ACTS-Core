#include <iostream>
#include "gen/MaterialMap.pb.h"
#include <google/protobuf/text_format.h>
#include <fstream>
#include <memory>

#include <cstdlib>

//#include <fcntl.h>
#include <stdio.h>

#include <google/protobuf/message_lite.h>
#include <google/protobuf/io/zero_copy_stream.h>
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>

bool writeDelimitedTo(
    const google::protobuf::MessageLite& message,
    google::protobuf::io::ZeroCopyOutputStream* rawOutput) {
  // We create a new coded stream for each message.  Don't worry, this is fast.
  google::protobuf::io::CodedOutputStream output(rawOutput);

  // Write the size.
  const int size = message.ByteSize();
  output.WriteVarint32(size);

  uint8_t* buffer = output.GetDirectBufferForNBytesAndAdvance(size);
  if (buffer != NULL) {
    // Optimization:  The message fits in one buffer, so use the faster
    // direct-to-array serialization path.
    message.SerializeWithCachedSizesToArray(buffer);
  } else {
    // Slightly-slower path when the message is multiple buffers.
    message.SerializeWithCachedSizes(&output);
    if (output.HadError()) return false;
  }

  return true;
}

bool readDelimitedFrom(
    google::protobuf::io::ZeroCopyInputStream* rawInput,
    google::protobuf::MessageLite* message) {
  // We create a new coded stream for each message.  Don't worry, this is fast,
  // and it makes sure the 64MB total size limit is imposed per-message rather
  // than on the whole stream.  (See the CodedInputStream interface for more
  // info on this limit.)
  google::protobuf::io::CodedInputStream input(rawInput);

  // Read the size.
  uint32_t size;
  if (!input.ReadVarint32(&size)) return false;

  // Tell the stream not to read beyond that size.
  google::protobuf::io::CodedInputStream::Limit limit =
      input.PushLimit(size);

  // Parse the message.
  if (!message->MergeFromCodedStream(&input)) return false;
  if (!input.ConsumedEntireMessage()) return false;

  // Release the limit.
  input.PopLimit(limit);

  return true;
}

void write_map(std::string fname) {
  auto outfd = fopen(fname.c_str(), "w+");

  auto outstream = std::make_unique<google::protobuf::io::FileOutputStream>(fileno(outfd));

  size_t rows = 3;
  size_t cols = 3;

  size_t globidx = 0;

  for(size_t j=0;j<50;j++) {

    Acts::protobuf::MaterialMap matMap;

    matMap.set_rows(rows);
    matMap.set_cols(cols);
    matMap.set_geo_id(j);

    matMap.set_sen_id(j);


    for(size_t i=0;i<rows*cols;i++) {
      Acts::protobuf::MaterialMap::MaterialBin* matBin = matMap.add_bins();
      
      double val = globidx;
      globidx++;
      
      //std::cout << "add bin: " << i << " => " << val << " size: " << matMap.bins_size() << std::endl;

      matBin->set_thickness(val);
      matBin->set_x0(0);
      matBin->set_l0(0);
      matBin->set_a(0);
      matBin->set_z(0);
      matBin->set_rho(0);
    }

    writeDelimitedTo(matMap, outstream.get());
    //matMap.Clear();
  }
}

void read_map(std::string fname) {
  auto infd = fopen(fname.c_str(), "r");
  auto instream = std::make_unique<google::protobuf::io::FileInputStream>(fileno(infd));

  Acts::protobuf::MaterialMap matMap;

  while(readDelimitedFrom(instream.get(), &matMap)) {
    std::cout << "read material map" << std::endl;

    size_t rows = matMap.rows();
    size_t cols = matMap.cols();
    size_t nBins = rows * cols;

    std::cout << "nrows: " << matMap.rows() << std::endl;
    std::cout << "ncols: " << matMap.cols() << std::endl;
    std::cout << "geo_id: " << matMap.geo_id() << std::endl;
    std::cout << "sen_id: " << matMap.sen_id() << std::endl;

    std::cout << "nbins: " << matMap.bins_size() << std::endl;

    for(int i=0;i<matMap.bins_size();i++) {

      int row = i / rows;
      int col = i % rows;


      const Acts::protobuf::MaterialMap::MaterialBin &bin = matMap.bins(i);
      std::cout << " - " << row << "x" << col << " Material(X0=" << bin.x0() << " L0=" << bin.l0() << " A=" << bin.a() << " Z=" << bin.z() << " rho=" << bin.rho() << " t=" << bin.thickness() << ")" << std::endl;
    }

    matMap.Clear();
  }

}

int main() {
  std::cout << "Hallo welt" << std::endl;

  GOOGLE_PROTOBUF_VERIFY_VERSION;
  
  write_map("out.pb");


  read_map("out.pb");

}
