// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: MaterialMap.proto

#include "MaterialMap.pb.h"

#include <algorithm>

#include <google/protobuf/stubs/common.h>
#include <google/protobuf/stubs/port.h>
#include <google/protobuf/stubs/once.h>
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/wire_format_lite_inl.h>
#include <google/protobuf/descriptor.h>
#include <google/protobuf/generated_message_reflection.h>
#include <google/protobuf/reflection_ops.h>
#include <google/protobuf/wire_format.h>
// This is a temporary google only hack
#ifdef GOOGLE_PROTOBUF_ENFORCE_UNIQUENESS
#include "third_party/protobuf/version.h"
#endif
// @@protoc_insertion_point(includes)
namespace Acts {
namespace protobuf {
class MaterialMap_MaterialBinDefaultTypeInternal {
 public:
  ::google::protobuf::internal::ExplicitlyConstructed<MaterialMap_MaterialBin>
      _instance;
} _MaterialMap_MaterialBin_default_instance_;
class MaterialMapDefaultTypeInternal {
 public:
  ::google::protobuf::internal::ExplicitlyConstructed<MaterialMap>
      _instance;
} _MaterialMap_default_instance_;
}  // namespace protobuf
}  // namespace Acts
namespace protobuf_MaterialMap_2eproto {
void InitDefaultsMaterialMap_MaterialBinImpl() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

#ifdef GOOGLE_PROTOBUF_ENFORCE_UNIQUENESS
  ::google::protobuf::internal::InitProtobufDefaultsForceUnique();
#else
  ::google::protobuf::internal::InitProtobufDefaults();
#endif  // GOOGLE_PROTOBUF_ENFORCE_UNIQUENESS
  {
    void* ptr = &::Acts::protobuf::_MaterialMap_MaterialBin_default_instance_;
    new (ptr) ::Acts::protobuf::MaterialMap_MaterialBin();
    ::google::protobuf::internal::OnShutdownDestroyMessage(ptr);
  }
  ::Acts::protobuf::MaterialMap_MaterialBin::InitAsDefaultInstance();
}

void InitDefaultsMaterialMap_MaterialBin() {
  static GOOGLE_PROTOBUF_DECLARE_ONCE(once);
  ::google::protobuf::GoogleOnceInit(&once, &InitDefaultsMaterialMap_MaterialBinImpl);
}

void InitDefaultsMaterialMapImpl() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

#ifdef GOOGLE_PROTOBUF_ENFORCE_UNIQUENESS
  ::google::protobuf::internal::InitProtobufDefaultsForceUnique();
#else
  ::google::protobuf::internal::InitProtobufDefaults();
#endif  // GOOGLE_PROTOBUF_ENFORCE_UNIQUENESS
  protobuf_MaterialMap_2eproto::InitDefaultsMaterialMap_MaterialBin();
  {
    void* ptr = &::Acts::protobuf::_MaterialMap_default_instance_;
    new (ptr) ::Acts::protobuf::MaterialMap();
    ::google::protobuf::internal::OnShutdownDestroyMessage(ptr);
  }
  ::Acts::protobuf::MaterialMap::InitAsDefaultInstance();
}

void InitDefaultsMaterialMap() {
  static GOOGLE_PROTOBUF_DECLARE_ONCE(once);
  ::google::protobuf::GoogleOnceInit(&once, &InitDefaultsMaterialMapImpl);
}

::google::protobuf::Metadata file_level_metadata[2];

const ::google::protobuf::uint32 TableStruct::offsets[] GOOGLE_PROTOBUF_ATTRIBUTE_SECTION_VARIABLE(protodesc_cold) = {
  ~0u,  // no _has_bits_
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::Acts::protobuf::MaterialMap_MaterialBin, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::Acts::protobuf::MaterialMap_MaterialBin, thickness_),
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::Acts::protobuf::MaterialMap_MaterialBin, x0_),
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::Acts::protobuf::MaterialMap_MaterialBin, l0_),
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::Acts::protobuf::MaterialMap_MaterialBin, a_),
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::Acts::protobuf::MaterialMap_MaterialBin, z_),
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::Acts::protobuf::MaterialMap_MaterialBin, rho_),
  ~0u,  // no _has_bits_
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::Acts::protobuf::MaterialMap, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::Acts::protobuf::MaterialMap, rows_),
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::Acts::protobuf::MaterialMap, cols_),
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::Acts::protobuf::MaterialMap, geo_id_),
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::Acts::protobuf::MaterialMap, vol_id_),
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::Acts::protobuf::MaterialMap, lay_id_),
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::Acts::protobuf::MaterialMap, app_id_),
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::Acts::protobuf::MaterialMap, sen_id_),
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::Acts::protobuf::MaterialMap, bins_),
};
static const ::google::protobuf::internal::MigrationSchema schemas[] GOOGLE_PROTOBUF_ATTRIBUTE_SECTION_VARIABLE(protodesc_cold) = {
  { 0, -1, sizeof(::Acts::protobuf::MaterialMap_MaterialBin)},
  { 11, -1, sizeof(::Acts::protobuf::MaterialMap)},
};

static ::google::protobuf::Message const * const file_default_instances[] = {
  reinterpret_cast<const ::google::protobuf::Message*>(&::Acts::protobuf::_MaterialMap_MaterialBin_default_instance_),
  reinterpret_cast<const ::google::protobuf::Message*>(&::Acts::protobuf::_MaterialMap_default_instance_),
};

void protobuf_AssignDescriptors() {
  AddDescriptors();
  ::google::protobuf::MessageFactory* factory = NULL;
  AssignDescriptors(
      "MaterialMap.proto", schemas, file_default_instances, TableStruct::offsets, factory,
      file_level_metadata, NULL, NULL);
}

void protobuf_AssignDescriptorsOnce() {
  static GOOGLE_PROTOBUF_DECLARE_ONCE(once);
  ::google::protobuf::GoogleOnceInit(&once, &protobuf_AssignDescriptors);
}

void protobuf_RegisterTypes(const ::std::string&) GOOGLE_PROTOBUF_ATTRIBUTE_COLD;
void protobuf_RegisterTypes(const ::std::string&) {
  protobuf_AssignDescriptorsOnce();
  ::google::protobuf::internal::RegisterAllTypes(file_level_metadata, 2);
}

void AddDescriptorsImpl() {
  InitDefaults();
  static const char descriptor[] GOOGLE_PROTOBUF_ATTRIBUTE_SECTION_VARIABLE(protodesc_cold) = {
      "\n\021MaterialMap.proto\022\rActs.protobuf\"\214\002\n\013M"
      "aterialMap\022\014\n\004rows\030\001 \001(\r\022\014\n\004cols\030\002 \001(\r\022\016"
      "\n\006geo_id\030\003 \001(\004\022\016\n\006vol_id\030\004 \001(\005\022\016\n\006lay_id"
      "\030\005 \001(\005\022\016\n\006app_id\030\006 \001(\005\022\016\n\006sen_id\030\007 \001(\005\0224"
      "\n\004bins\030\010 \003(\0132&.Acts.protobuf.MaterialMap"
      ".MaterialBin\032[\n\013MaterialBin\022\021\n\tthickness"
      "\030\001 \001(\001\022\n\n\002X0\030\002 \001(\001\022\n\n\002L0\030\003 \001(\001\022\t\n\001A\030\004 \001("
      "\001\022\t\n\001Z\030\005 \001(\001\022\013\n\003rho\030\006 \001(\001b\006proto3"
  };
  ::google::protobuf::DescriptorPool::InternalAddGeneratedFile(
      descriptor, 313);
  ::google::protobuf::MessageFactory::InternalRegisterGeneratedFile(
    "MaterialMap.proto", &protobuf_RegisterTypes);
}

void AddDescriptors() {
  static GOOGLE_PROTOBUF_DECLARE_ONCE(once);
  ::google::protobuf::GoogleOnceInit(&once, &AddDescriptorsImpl);
}
// Force AddDescriptors() to be called at dynamic initialization time.
struct StaticDescriptorInitializer {
  StaticDescriptorInitializer() {
    AddDescriptors();
  }
} static_descriptor_initializer;
}  // namespace protobuf_MaterialMap_2eproto
namespace Acts {
namespace protobuf {

// ===================================================================

void MaterialMap_MaterialBin::InitAsDefaultInstance() {
}
#if !defined(_MSC_VER) || _MSC_VER >= 1900
const int MaterialMap_MaterialBin::kThicknessFieldNumber;
const int MaterialMap_MaterialBin::kX0FieldNumber;
const int MaterialMap_MaterialBin::kL0FieldNumber;
const int MaterialMap_MaterialBin::kAFieldNumber;
const int MaterialMap_MaterialBin::kZFieldNumber;
const int MaterialMap_MaterialBin::kRhoFieldNumber;
#endif  // !defined(_MSC_VER) || _MSC_VER >= 1900

MaterialMap_MaterialBin::MaterialMap_MaterialBin()
  : ::google::protobuf::Message(), _internal_metadata_(NULL) {
  if (GOOGLE_PREDICT_TRUE(this != internal_default_instance())) {
    ::protobuf_MaterialMap_2eproto::InitDefaultsMaterialMap_MaterialBin();
  }
  SharedCtor();
  // @@protoc_insertion_point(constructor:Acts.protobuf.MaterialMap.MaterialBin)
}
MaterialMap_MaterialBin::MaterialMap_MaterialBin(const MaterialMap_MaterialBin& from)
  : ::google::protobuf::Message(),
      _internal_metadata_(NULL),
      _cached_size_(0) {
  _internal_metadata_.MergeFrom(from._internal_metadata_);
  ::memcpy(&thickness_, &from.thickness_,
    static_cast<size_t>(reinterpret_cast<char*>(&rho_) -
    reinterpret_cast<char*>(&thickness_)) + sizeof(rho_));
  // @@protoc_insertion_point(copy_constructor:Acts.protobuf.MaterialMap.MaterialBin)
}

void MaterialMap_MaterialBin::SharedCtor() {
  ::memset(&thickness_, 0, static_cast<size_t>(
      reinterpret_cast<char*>(&rho_) -
      reinterpret_cast<char*>(&thickness_)) + sizeof(rho_));
  _cached_size_ = 0;
}

MaterialMap_MaterialBin::~MaterialMap_MaterialBin() {
  // @@protoc_insertion_point(destructor:Acts.protobuf.MaterialMap.MaterialBin)
  SharedDtor();
}

void MaterialMap_MaterialBin::SharedDtor() {
}

void MaterialMap_MaterialBin::SetCachedSize(int size) const {
  GOOGLE_SAFE_CONCURRENT_WRITES_BEGIN();
  _cached_size_ = size;
  GOOGLE_SAFE_CONCURRENT_WRITES_END();
}
const ::google::protobuf::Descriptor* MaterialMap_MaterialBin::descriptor() {
  ::protobuf_MaterialMap_2eproto::protobuf_AssignDescriptorsOnce();
  return ::protobuf_MaterialMap_2eproto::file_level_metadata[kIndexInFileMessages].descriptor;
}

const MaterialMap_MaterialBin& MaterialMap_MaterialBin::default_instance() {
  ::protobuf_MaterialMap_2eproto::InitDefaultsMaterialMap_MaterialBin();
  return *internal_default_instance();
}

MaterialMap_MaterialBin* MaterialMap_MaterialBin::New(::google::protobuf::Arena* arena) const {
  MaterialMap_MaterialBin* n = new MaterialMap_MaterialBin;
  if (arena != NULL) {
    arena->Own(n);
  }
  return n;
}

void MaterialMap_MaterialBin::Clear() {
// @@protoc_insertion_point(message_clear_start:Acts.protobuf.MaterialMap.MaterialBin)
  ::google::protobuf::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  ::memset(&thickness_, 0, static_cast<size_t>(
      reinterpret_cast<char*>(&rho_) -
      reinterpret_cast<char*>(&thickness_)) + sizeof(rho_));
  _internal_metadata_.Clear();
}

bool MaterialMap_MaterialBin::MergePartialFromCodedStream(
    ::google::protobuf::io::CodedInputStream* input) {
#define DO_(EXPRESSION) if (!GOOGLE_PREDICT_TRUE(EXPRESSION)) goto failure
  ::google::protobuf::uint32 tag;
  // @@protoc_insertion_point(parse_start:Acts.protobuf.MaterialMap.MaterialBin)
  for (;;) {
    ::std::pair< ::google::protobuf::uint32, bool> p = input->ReadTagWithCutoffNoLastTag(127u);
    tag = p.first;
    if (!p.second) goto handle_unusual;
    switch (::google::protobuf::internal::WireFormatLite::GetTagFieldNumber(tag)) {
      // double thickness = 1;
      case 1: {
        if (static_cast< ::google::protobuf::uint8>(tag) ==
            static_cast< ::google::protobuf::uint8>(9u /* 9 & 0xFF */)) {

          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   double, ::google::protobuf::internal::WireFormatLite::TYPE_DOUBLE>(
                 input, &thickness_)));
        } else {
          goto handle_unusual;
        }
        break;
      }

      // double X0 = 2;
      case 2: {
        if (static_cast< ::google::protobuf::uint8>(tag) ==
            static_cast< ::google::protobuf::uint8>(17u /* 17 & 0xFF */)) {

          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   double, ::google::protobuf::internal::WireFormatLite::TYPE_DOUBLE>(
                 input, &x0_)));
        } else {
          goto handle_unusual;
        }
        break;
      }

      // double L0 = 3;
      case 3: {
        if (static_cast< ::google::protobuf::uint8>(tag) ==
            static_cast< ::google::protobuf::uint8>(25u /* 25 & 0xFF */)) {

          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   double, ::google::protobuf::internal::WireFormatLite::TYPE_DOUBLE>(
                 input, &l0_)));
        } else {
          goto handle_unusual;
        }
        break;
      }

      // double A = 4;
      case 4: {
        if (static_cast< ::google::protobuf::uint8>(tag) ==
            static_cast< ::google::protobuf::uint8>(33u /* 33 & 0xFF */)) {

          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   double, ::google::protobuf::internal::WireFormatLite::TYPE_DOUBLE>(
                 input, &a_)));
        } else {
          goto handle_unusual;
        }
        break;
      }

      // double Z = 5;
      case 5: {
        if (static_cast< ::google::protobuf::uint8>(tag) ==
            static_cast< ::google::protobuf::uint8>(41u /* 41 & 0xFF */)) {

          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   double, ::google::protobuf::internal::WireFormatLite::TYPE_DOUBLE>(
                 input, &z_)));
        } else {
          goto handle_unusual;
        }
        break;
      }

      // double rho = 6;
      case 6: {
        if (static_cast< ::google::protobuf::uint8>(tag) ==
            static_cast< ::google::protobuf::uint8>(49u /* 49 & 0xFF */)) {

          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   double, ::google::protobuf::internal::WireFormatLite::TYPE_DOUBLE>(
                 input, &rho_)));
        } else {
          goto handle_unusual;
        }
        break;
      }

      default: {
      handle_unusual:
        if (tag == 0) {
          goto success;
        }
        DO_(::google::protobuf::internal::WireFormat::SkipField(
              input, tag, _internal_metadata_.mutable_unknown_fields()));
        break;
      }
    }
  }
success:
  // @@protoc_insertion_point(parse_success:Acts.protobuf.MaterialMap.MaterialBin)
  return true;
failure:
  // @@protoc_insertion_point(parse_failure:Acts.protobuf.MaterialMap.MaterialBin)
  return false;
#undef DO_
}

void MaterialMap_MaterialBin::SerializeWithCachedSizes(
    ::google::protobuf::io::CodedOutputStream* output) const {
  // @@protoc_insertion_point(serialize_start:Acts.protobuf.MaterialMap.MaterialBin)
  ::google::protobuf::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  // double thickness = 1;
  if (this->thickness() != 0) {
    ::google::protobuf::internal::WireFormatLite::WriteDouble(1, this->thickness(), output);
  }

  // double X0 = 2;
  if (this->x0() != 0) {
    ::google::protobuf::internal::WireFormatLite::WriteDouble(2, this->x0(), output);
  }

  // double L0 = 3;
  if (this->l0() != 0) {
    ::google::protobuf::internal::WireFormatLite::WriteDouble(3, this->l0(), output);
  }

  // double A = 4;
  if (this->a() != 0) {
    ::google::protobuf::internal::WireFormatLite::WriteDouble(4, this->a(), output);
  }

  // double Z = 5;
  if (this->z() != 0) {
    ::google::protobuf::internal::WireFormatLite::WriteDouble(5, this->z(), output);
  }

  // double rho = 6;
  if (this->rho() != 0) {
    ::google::protobuf::internal::WireFormatLite::WriteDouble(6, this->rho(), output);
  }

  if ((_internal_metadata_.have_unknown_fields() &&  ::google::protobuf::internal::GetProto3PreserveUnknownsDefault())) {
    ::google::protobuf::internal::WireFormat::SerializeUnknownFields(
        (::google::protobuf::internal::GetProto3PreserveUnknownsDefault()   ? _internal_metadata_.unknown_fields()   : _internal_metadata_.default_instance()), output);
  }
  // @@protoc_insertion_point(serialize_end:Acts.protobuf.MaterialMap.MaterialBin)
}

::google::protobuf::uint8* MaterialMap_MaterialBin::InternalSerializeWithCachedSizesToArray(
    bool deterministic, ::google::protobuf::uint8* target) const {
  (void)deterministic; // Unused
  // @@protoc_insertion_point(serialize_to_array_start:Acts.protobuf.MaterialMap.MaterialBin)
  ::google::protobuf::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  // double thickness = 1;
  if (this->thickness() != 0) {
    target = ::google::protobuf::internal::WireFormatLite::WriteDoubleToArray(1, this->thickness(), target);
  }

  // double X0 = 2;
  if (this->x0() != 0) {
    target = ::google::protobuf::internal::WireFormatLite::WriteDoubleToArray(2, this->x0(), target);
  }

  // double L0 = 3;
  if (this->l0() != 0) {
    target = ::google::protobuf::internal::WireFormatLite::WriteDoubleToArray(3, this->l0(), target);
  }

  // double A = 4;
  if (this->a() != 0) {
    target = ::google::protobuf::internal::WireFormatLite::WriteDoubleToArray(4, this->a(), target);
  }

  // double Z = 5;
  if (this->z() != 0) {
    target = ::google::protobuf::internal::WireFormatLite::WriteDoubleToArray(5, this->z(), target);
  }

  // double rho = 6;
  if (this->rho() != 0) {
    target = ::google::protobuf::internal::WireFormatLite::WriteDoubleToArray(6, this->rho(), target);
  }

  if ((_internal_metadata_.have_unknown_fields() &&  ::google::protobuf::internal::GetProto3PreserveUnknownsDefault())) {
    target = ::google::protobuf::internal::WireFormat::SerializeUnknownFieldsToArray(
        (::google::protobuf::internal::GetProto3PreserveUnknownsDefault()   ? _internal_metadata_.unknown_fields()   : _internal_metadata_.default_instance()), target);
  }
  // @@protoc_insertion_point(serialize_to_array_end:Acts.protobuf.MaterialMap.MaterialBin)
  return target;
}

size_t MaterialMap_MaterialBin::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:Acts.protobuf.MaterialMap.MaterialBin)
  size_t total_size = 0;

  if ((_internal_metadata_.have_unknown_fields() &&  ::google::protobuf::internal::GetProto3PreserveUnknownsDefault())) {
    total_size +=
      ::google::protobuf::internal::WireFormat::ComputeUnknownFieldsSize(
        (::google::protobuf::internal::GetProto3PreserveUnknownsDefault()   ? _internal_metadata_.unknown_fields()   : _internal_metadata_.default_instance()));
  }
  // double thickness = 1;
  if (this->thickness() != 0) {
    total_size += 1 + 8;
  }

  // double X0 = 2;
  if (this->x0() != 0) {
    total_size += 1 + 8;
  }

  // double L0 = 3;
  if (this->l0() != 0) {
    total_size += 1 + 8;
  }

  // double A = 4;
  if (this->a() != 0) {
    total_size += 1 + 8;
  }

  // double Z = 5;
  if (this->z() != 0) {
    total_size += 1 + 8;
  }

  // double rho = 6;
  if (this->rho() != 0) {
    total_size += 1 + 8;
  }

  int cached_size = ::google::protobuf::internal::ToCachedSize(total_size);
  GOOGLE_SAFE_CONCURRENT_WRITES_BEGIN();
  _cached_size_ = cached_size;
  GOOGLE_SAFE_CONCURRENT_WRITES_END();
  return total_size;
}

void MaterialMap_MaterialBin::MergeFrom(const ::google::protobuf::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:Acts.protobuf.MaterialMap.MaterialBin)
  GOOGLE_DCHECK_NE(&from, this);
  const MaterialMap_MaterialBin* source =
      ::google::protobuf::internal::DynamicCastToGenerated<const MaterialMap_MaterialBin>(
          &from);
  if (source == NULL) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:Acts.protobuf.MaterialMap.MaterialBin)
    ::google::protobuf::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:Acts.protobuf.MaterialMap.MaterialBin)
    MergeFrom(*source);
  }
}

void MaterialMap_MaterialBin::MergeFrom(const MaterialMap_MaterialBin& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:Acts.protobuf.MaterialMap.MaterialBin)
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom(from._internal_metadata_);
  ::google::protobuf::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  if (from.thickness() != 0) {
    set_thickness(from.thickness());
  }
  if (from.x0() != 0) {
    set_x0(from.x0());
  }
  if (from.l0() != 0) {
    set_l0(from.l0());
  }
  if (from.a() != 0) {
    set_a(from.a());
  }
  if (from.z() != 0) {
    set_z(from.z());
  }
  if (from.rho() != 0) {
    set_rho(from.rho());
  }
}

void MaterialMap_MaterialBin::CopyFrom(const ::google::protobuf::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:Acts.protobuf.MaterialMap.MaterialBin)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void MaterialMap_MaterialBin::CopyFrom(const MaterialMap_MaterialBin& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:Acts.protobuf.MaterialMap.MaterialBin)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool MaterialMap_MaterialBin::IsInitialized() const {
  return true;
}

void MaterialMap_MaterialBin::Swap(MaterialMap_MaterialBin* other) {
  if (other == this) return;
  InternalSwap(other);
}
void MaterialMap_MaterialBin::InternalSwap(MaterialMap_MaterialBin* other) {
  using std::swap;
  swap(thickness_, other->thickness_);
  swap(x0_, other->x0_);
  swap(l0_, other->l0_);
  swap(a_, other->a_);
  swap(z_, other->z_);
  swap(rho_, other->rho_);
  _internal_metadata_.Swap(&other->_internal_metadata_);
  swap(_cached_size_, other->_cached_size_);
}

::google::protobuf::Metadata MaterialMap_MaterialBin::GetMetadata() const {
  protobuf_MaterialMap_2eproto::protobuf_AssignDescriptorsOnce();
  return ::protobuf_MaterialMap_2eproto::file_level_metadata[kIndexInFileMessages];
}


// ===================================================================

void MaterialMap::InitAsDefaultInstance() {
}
#if !defined(_MSC_VER) || _MSC_VER >= 1900
const int MaterialMap::kRowsFieldNumber;
const int MaterialMap::kColsFieldNumber;
const int MaterialMap::kGeoIdFieldNumber;
const int MaterialMap::kVolIdFieldNumber;
const int MaterialMap::kLayIdFieldNumber;
const int MaterialMap::kAppIdFieldNumber;
const int MaterialMap::kSenIdFieldNumber;
const int MaterialMap::kBinsFieldNumber;
#endif  // !defined(_MSC_VER) || _MSC_VER >= 1900

MaterialMap::MaterialMap()
  : ::google::protobuf::Message(), _internal_metadata_(NULL) {
  if (GOOGLE_PREDICT_TRUE(this != internal_default_instance())) {
    ::protobuf_MaterialMap_2eproto::InitDefaultsMaterialMap();
  }
  SharedCtor();
  // @@protoc_insertion_point(constructor:Acts.protobuf.MaterialMap)
}
MaterialMap::MaterialMap(const MaterialMap& from)
  : ::google::protobuf::Message(),
      _internal_metadata_(NULL),
      bins_(from.bins_),
      _cached_size_(0) {
  _internal_metadata_.MergeFrom(from._internal_metadata_);
  ::memcpy(&rows_, &from.rows_,
    static_cast<size_t>(reinterpret_cast<char*>(&sen_id_) -
    reinterpret_cast<char*>(&rows_)) + sizeof(sen_id_));
  // @@protoc_insertion_point(copy_constructor:Acts.protobuf.MaterialMap)
}

void MaterialMap::SharedCtor() {
  ::memset(&rows_, 0, static_cast<size_t>(
      reinterpret_cast<char*>(&sen_id_) -
      reinterpret_cast<char*>(&rows_)) + sizeof(sen_id_));
  _cached_size_ = 0;
}

MaterialMap::~MaterialMap() {
  // @@protoc_insertion_point(destructor:Acts.protobuf.MaterialMap)
  SharedDtor();
}

void MaterialMap::SharedDtor() {
}

void MaterialMap::SetCachedSize(int size) const {
  GOOGLE_SAFE_CONCURRENT_WRITES_BEGIN();
  _cached_size_ = size;
  GOOGLE_SAFE_CONCURRENT_WRITES_END();
}
const ::google::protobuf::Descriptor* MaterialMap::descriptor() {
  ::protobuf_MaterialMap_2eproto::protobuf_AssignDescriptorsOnce();
  return ::protobuf_MaterialMap_2eproto::file_level_metadata[kIndexInFileMessages].descriptor;
}

const MaterialMap& MaterialMap::default_instance() {
  ::protobuf_MaterialMap_2eproto::InitDefaultsMaterialMap();
  return *internal_default_instance();
}

MaterialMap* MaterialMap::New(::google::protobuf::Arena* arena) const {
  MaterialMap* n = new MaterialMap;
  if (arena != NULL) {
    arena->Own(n);
  }
  return n;
}

void MaterialMap::Clear() {
// @@protoc_insertion_point(message_clear_start:Acts.protobuf.MaterialMap)
  ::google::protobuf::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  bins_.Clear();
  ::memset(&rows_, 0, static_cast<size_t>(
      reinterpret_cast<char*>(&sen_id_) -
      reinterpret_cast<char*>(&rows_)) + sizeof(sen_id_));
  _internal_metadata_.Clear();
}

bool MaterialMap::MergePartialFromCodedStream(
    ::google::protobuf::io::CodedInputStream* input) {
#define DO_(EXPRESSION) if (!GOOGLE_PREDICT_TRUE(EXPRESSION)) goto failure
  ::google::protobuf::uint32 tag;
  // @@protoc_insertion_point(parse_start:Acts.protobuf.MaterialMap)
  for (;;) {
    ::std::pair< ::google::protobuf::uint32, bool> p = input->ReadTagWithCutoffNoLastTag(127u);
    tag = p.first;
    if (!p.second) goto handle_unusual;
    switch (::google::protobuf::internal::WireFormatLite::GetTagFieldNumber(tag)) {
      // uint32 rows = 1;
      case 1: {
        if (static_cast< ::google::protobuf::uint8>(tag) ==
            static_cast< ::google::protobuf::uint8>(8u /* 8 & 0xFF */)) {

          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   ::google::protobuf::uint32, ::google::protobuf::internal::WireFormatLite::TYPE_UINT32>(
                 input, &rows_)));
        } else {
          goto handle_unusual;
        }
        break;
      }

      // uint32 cols = 2;
      case 2: {
        if (static_cast< ::google::protobuf::uint8>(tag) ==
            static_cast< ::google::protobuf::uint8>(16u /* 16 & 0xFF */)) {

          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   ::google::protobuf::uint32, ::google::protobuf::internal::WireFormatLite::TYPE_UINT32>(
                 input, &cols_)));
        } else {
          goto handle_unusual;
        }
        break;
      }

      // uint64 geo_id = 3;
      case 3: {
        if (static_cast< ::google::protobuf::uint8>(tag) ==
            static_cast< ::google::protobuf::uint8>(24u /* 24 & 0xFF */)) {

          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   ::google::protobuf::uint64, ::google::protobuf::internal::WireFormatLite::TYPE_UINT64>(
                 input, &geo_id_)));
        } else {
          goto handle_unusual;
        }
        break;
      }

      // int32 vol_id = 4;
      case 4: {
        if (static_cast< ::google::protobuf::uint8>(tag) ==
            static_cast< ::google::protobuf::uint8>(32u /* 32 & 0xFF */)) {

          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   ::google::protobuf::int32, ::google::protobuf::internal::WireFormatLite::TYPE_INT32>(
                 input, &vol_id_)));
        } else {
          goto handle_unusual;
        }
        break;
      }

      // int32 lay_id = 5;
      case 5: {
        if (static_cast< ::google::protobuf::uint8>(tag) ==
            static_cast< ::google::protobuf::uint8>(40u /* 40 & 0xFF */)) {

          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   ::google::protobuf::int32, ::google::protobuf::internal::WireFormatLite::TYPE_INT32>(
                 input, &lay_id_)));
        } else {
          goto handle_unusual;
        }
        break;
      }

      // int32 app_id = 6;
      case 6: {
        if (static_cast< ::google::protobuf::uint8>(tag) ==
            static_cast< ::google::protobuf::uint8>(48u /* 48 & 0xFF */)) {

          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   ::google::protobuf::int32, ::google::protobuf::internal::WireFormatLite::TYPE_INT32>(
                 input, &app_id_)));
        } else {
          goto handle_unusual;
        }
        break;
      }

      // int32 sen_id = 7;
      case 7: {
        if (static_cast< ::google::protobuf::uint8>(tag) ==
            static_cast< ::google::protobuf::uint8>(56u /* 56 & 0xFF */)) {

          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   ::google::protobuf::int32, ::google::protobuf::internal::WireFormatLite::TYPE_INT32>(
                 input, &sen_id_)));
        } else {
          goto handle_unusual;
        }
        break;
      }

      // repeated .Acts.protobuf.MaterialMap.MaterialBin bins = 8;
      case 8: {
        if (static_cast< ::google::protobuf::uint8>(tag) ==
            static_cast< ::google::protobuf::uint8>(66u /* 66 & 0xFF */)) {
          DO_(::google::protobuf::internal::WireFormatLite::ReadMessage(input, add_bins()));
        } else {
          goto handle_unusual;
        }
        break;
      }

      default: {
      handle_unusual:
        if (tag == 0) {
          goto success;
        }
        DO_(::google::protobuf::internal::WireFormat::SkipField(
              input, tag, _internal_metadata_.mutable_unknown_fields()));
        break;
      }
    }
  }
success:
  // @@protoc_insertion_point(parse_success:Acts.protobuf.MaterialMap)
  return true;
failure:
  // @@protoc_insertion_point(parse_failure:Acts.protobuf.MaterialMap)
  return false;
#undef DO_
}

void MaterialMap::SerializeWithCachedSizes(
    ::google::protobuf::io::CodedOutputStream* output) const {
  // @@protoc_insertion_point(serialize_start:Acts.protobuf.MaterialMap)
  ::google::protobuf::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  // uint32 rows = 1;
  if (this->rows() != 0) {
    ::google::protobuf::internal::WireFormatLite::WriteUInt32(1, this->rows(), output);
  }

  // uint32 cols = 2;
  if (this->cols() != 0) {
    ::google::protobuf::internal::WireFormatLite::WriteUInt32(2, this->cols(), output);
  }

  // uint64 geo_id = 3;
  if (this->geo_id() != 0) {
    ::google::protobuf::internal::WireFormatLite::WriteUInt64(3, this->geo_id(), output);
  }

  // int32 vol_id = 4;
  if (this->vol_id() != 0) {
    ::google::protobuf::internal::WireFormatLite::WriteInt32(4, this->vol_id(), output);
  }

  // int32 lay_id = 5;
  if (this->lay_id() != 0) {
    ::google::protobuf::internal::WireFormatLite::WriteInt32(5, this->lay_id(), output);
  }

  // int32 app_id = 6;
  if (this->app_id() != 0) {
    ::google::protobuf::internal::WireFormatLite::WriteInt32(6, this->app_id(), output);
  }

  // int32 sen_id = 7;
  if (this->sen_id() != 0) {
    ::google::protobuf::internal::WireFormatLite::WriteInt32(7, this->sen_id(), output);
  }

  // repeated .Acts.protobuf.MaterialMap.MaterialBin bins = 8;
  for (unsigned int i = 0,
      n = static_cast<unsigned int>(this->bins_size()); i < n; i++) {
    ::google::protobuf::internal::WireFormatLite::WriteMessageMaybeToArray(
      8, this->bins(static_cast<int>(i)), output);
  }

  if ((_internal_metadata_.have_unknown_fields() &&  ::google::protobuf::internal::GetProto3PreserveUnknownsDefault())) {
    ::google::protobuf::internal::WireFormat::SerializeUnknownFields(
        (::google::protobuf::internal::GetProto3PreserveUnknownsDefault()   ? _internal_metadata_.unknown_fields()   : _internal_metadata_.default_instance()), output);
  }
  // @@protoc_insertion_point(serialize_end:Acts.protobuf.MaterialMap)
}

::google::protobuf::uint8* MaterialMap::InternalSerializeWithCachedSizesToArray(
    bool deterministic, ::google::protobuf::uint8* target) const {
  (void)deterministic; // Unused
  // @@protoc_insertion_point(serialize_to_array_start:Acts.protobuf.MaterialMap)
  ::google::protobuf::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  // uint32 rows = 1;
  if (this->rows() != 0) {
    target = ::google::protobuf::internal::WireFormatLite::WriteUInt32ToArray(1, this->rows(), target);
  }

  // uint32 cols = 2;
  if (this->cols() != 0) {
    target = ::google::protobuf::internal::WireFormatLite::WriteUInt32ToArray(2, this->cols(), target);
  }

  // uint64 geo_id = 3;
  if (this->geo_id() != 0) {
    target = ::google::protobuf::internal::WireFormatLite::WriteUInt64ToArray(3, this->geo_id(), target);
  }

  // int32 vol_id = 4;
  if (this->vol_id() != 0) {
    target = ::google::protobuf::internal::WireFormatLite::WriteInt32ToArray(4, this->vol_id(), target);
  }

  // int32 lay_id = 5;
  if (this->lay_id() != 0) {
    target = ::google::protobuf::internal::WireFormatLite::WriteInt32ToArray(5, this->lay_id(), target);
  }

  // int32 app_id = 6;
  if (this->app_id() != 0) {
    target = ::google::protobuf::internal::WireFormatLite::WriteInt32ToArray(6, this->app_id(), target);
  }

  // int32 sen_id = 7;
  if (this->sen_id() != 0) {
    target = ::google::protobuf::internal::WireFormatLite::WriteInt32ToArray(7, this->sen_id(), target);
  }

  // repeated .Acts.protobuf.MaterialMap.MaterialBin bins = 8;
  for (unsigned int i = 0,
      n = static_cast<unsigned int>(this->bins_size()); i < n; i++) {
    target = ::google::protobuf::internal::WireFormatLite::
      InternalWriteMessageToArray(
        8, this->bins(static_cast<int>(i)), deterministic, target);
  }

  if ((_internal_metadata_.have_unknown_fields() &&  ::google::protobuf::internal::GetProto3PreserveUnknownsDefault())) {
    target = ::google::protobuf::internal::WireFormat::SerializeUnknownFieldsToArray(
        (::google::protobuf::internal::GetProto3PreserveUnknownsDefault()   ? _internal_metadata_.unknown_fields()   : _internal_metadata_.default_instance()), target);
  }
  // @@protoc_insertion_point(serialize_to_array_end:Acts.protobuf.MaterialMap)
  return target;
}

size_t MaterialMap::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:Acts.protobuf.MaterialMap)
  size_t total_size = 0;

  if ((_internal_metadata_.have_unknown_fields() &&  ::google::protobuf::internal::GetProto3PreserveUnknownsDefault())) {
    total_size +=
      ::google::protobuf::internal::WireFormat::ComputeUnknownFieldsSize(
        (::google::protobuf::internal::GetProto3PreserveUnknownsDefault()   ? _internal_metadata_.unknown_fields()   : _internal_metadata_.default_instance()));
  }
  // repeated .Acts.protobuf.MaterialMap.MaterialBin bins = 8;
  {
    unsigned int count = static_cast<unsigned int>(this->bins_size());
    total_size += 1UL * count;
    for (unsigned int i = 0; i < count; i++) {
      total_size +=
        ::google::protobuf::internal::WireFormatLite::MessageSize(
          this->bins(static_cast<int>(i)));
    }
  }

  // uint32 rows = 1;
  if (this->rows() != 0) {
    total_size += 1 +
      ::google::protobuf::internal::WireFormatLite::UInt32Size(
        this->rows());
  }

  // uint32 cols = 2;
  if (this->cols() != 0) {
    total_size += 1 +
      ::google::protobuf::internal::WireFormatLite::UInt32Size(
        this->cols());
  }

  // uint64 geo_id = 3;
  if (this->geo_id() != 0) {
    total_size += 1 +
      ::google::protobuf::internal::WireFormatLite::UInt64Size(
        this->geo_id());
  }

  // int32 vol_id = 4;
  if (this->vol_id() != 0) {
    total_size += 1 +
      ::google::protobuf::internal::WireFormatLite::Int32Size(
        this->vol_id());
  }

  // int32 lay_id = 5;
  if (this->lay_id() != 0) {
    total_size += 1 +
      ::google::protobuf::internal::WireFormatLite::Int32Size(
        this->lay_id());
  }

  // int32 app_id = 6;
  if (this->app_id() != 0) {
    total_size += 1 +
      ::google::protobuf::internal::WireFormatLite::Int32Size(
        this->app_id());
  }

  // int32 sen_id = 7;
  if (this->sen_id() != 0) {
    total_size += 1 +
      ::google::protobuf::internal::WireFormatLite::Int32Size(
        this->sen_id());
  }

  int cached_size = ::google::protobuf::internal::ToCachedSize(total_size);
  GOOGLE_SAFE_CONCURRENT_WRITES_BEGIN();
  _cached_size_ = cached_size;
  GOOGLE_SAFE_CONCURRENT_WRITES_END();
  return total_size;
}

void MaterialMap::MergeFrom(const ::google::protobuf::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:Acts.protobuf.MaterialMap)
  GOOGLE_DCHECK_NE(&from, this);
  const MaterialMap* source =
      ::google::protobuf::internal::DynamicCastToGenerated<const MaterialMap>(
          &from);
  if (source == NULL) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:Acts.protobuf.MaterialMap)
    ::google::protobuf::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:Acts.protobuf.MaterialMap)
    MergeFrom(*source);
  }
}

void MaterialMap::MergeFrom(const MaterialMap& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:Acts.protobuf.MaterialMap)
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom(from._internal_metadata_);
  ::google::protobuf::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  bins_.MergeFrom(from.bins_);
  if (from.rows() != 0) {
    set_rows(from.rows());
  }
  if (from.cols() != 0) {
    set_cols(from.cols());
  }
  if (from.geo_id() != 0) {
    set_geo_id(from.geo_id());
  }
  if (from.vol_id() != 0) {
    set_vol_id(from.vol_id());
  }
  if (from.lay_id() != 0) {
    set_lay_id(from.lay_id());
  }
  if (from.app_id() != 0) {
    set_app_id(from.app_id());
  }
  if (from.sen_id() != 0) {
    set_sen_id(from.sen_id());
  }
}

void MaterialMap::CopyFrom(const ::google::protobuf::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:Acts.protobuf.MaterialMap)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void MaterialMap::CopyFrom(const MaterialMap& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:Acts.protobuf.MaterialMap)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool MaterialMap::IsInitialized() const {
  return true;
}

void MaterialMap::Swap(MaterialMap* other) {
  if (other == this) return;
  InternalSwap(other);
}
void MaterialMap::InternalSwap(MaterialMap* other) {
  using std::swap;
  bins_.InternalSwap(&other->bins_);
  swap(rows_, other->rows_);
  swap(cols_, other->cols_);
  swap(geo_id_, other->geo_id_);
  swap(vol_id_, other->vol_id_);
  swap(lay_id_, other->lay_id_);
  swap(app_id_, other->app_id_);
  swap(sen_id_, other->sen_id_);
  _internal_metadata_.Swap(&other->_internal_metadata_);
  swap(_cached_size_, other->_cached_size_);
}

::google::protobuf::Metadata MaterialMap::GetMetadata() const {
  protobuf_MaterialMap_2eproto::protobuf_AssignDescriptorsOnce();
  return ::protobuf_MaterialMap_2eproto::file_level_metadata[kIndexInFileMessages];
}


// @@protoc_insertion_point(namespace_scope)
}  // namespace protobuf
}  // namespace Acts

// @@protoc_insertion_point(global_scope)
