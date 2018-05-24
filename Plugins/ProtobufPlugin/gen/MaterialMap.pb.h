// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: MaterialMap.proto

#ifndef PROTOBUF_MaterialMap_2eproto__INCLUDED
#define PROTOBUF_MaterialMap_2eproto__INCLUDED

#include <string>

#include <google/protobuf/stubs/common.h>

#if GOOGLE_PROTOBUF_VERSION < 3005000
#error This file was generated by a newer version of protoc which is
#error incompatible with your Protocol Buffer headers.  Please update
#error your headers.
#endif
#if 3005001 < GOOGLE_PROTOBUF_MIN_PROTOC_VERSION
#error This file was generated by an older version of protoc which is
#error incompatible with your Protocol Buffer headers.  Please
#error regenerate this file with a newer version of protoc.
#endif

#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/arena.h>
#include <google/protobuf/arenastring.h>
#include <google/protobuf/generated_message_table_driven.h>
#include <google/protobuf/generated_message_util.h>
#include <google/protobuf/metadata.h>
#include <google/protobuf/message.h>
#include <google/protobuf/repeated_field.h>  // IWYU pragma: export
#include <google/protobuf/extension_set.h>  // IWYU pragma: export
#include <google/protobuf/unknown_field_set.h>
// @@protoc_insertion_point(includes)

namespace protobuf_MaterialMap_2eproto {
// Internal implementation detail -- do not use these members.
struct TableStruct {
  static const ::google::protobuf::internal::ParseTableField entries[];
  static const ::google::protobuf::internal::AuxillaryParseTableField aux[];
  static const ::google::protobuf::internal::ParseTable schema[3];
  static const ::google::protobuf::internal::FieldMetadata field_metadata[];
  static const ::google::protobuf::internal::SerializationTable serialization_table[];
  static const ::google::protobuf::uint32 offsets[];
};
void AddDescriptors();
void InitDefaultsMaterialMap_DimensionImpl();
void InitDefaultsMaterialMap_Dimension();
void InitDefaultsMaterialMap_MaterialPropertiesImpl();
void InitDefaultsMaterialMap_MaterialProperties();
void InitDefaultsMaterialMapImpl();
void InitDefaultsMaterialMap();
inline void InitDefaults() {
  InitDefaultsMaterialMap_Dimension();
  InitDefaultsMaterialMap_MaterialProperties();
  InitDefaultsMaterialMap();
}
}  // namespace protobuf_MaterialMap_2eproto
namespace Acts {
namespace protobuf {
class MaterialMap;
class MaterialMapDefaultTypeInternal;
extern MaterialMapDefaultTypeInternal _MaterialMap_default_instance_;
class MaterialMap_Dimension;
class MaterialMap_DimensionDefaultTypeInternal;
extern MaterialMap_DimensionDefaultTypeInternal _MaterialMap_Dimension_default_instance_;
class MaterialMap_MaterialProperties;
class MaterialMap_MaterialPropertiesDefaultTypeInternal;
extern MaterialMap_MaterialPropertiesDefaultTypeInternal _MaterialMap_MaterialProperties_default_instance_;
}  // namespace protobuf
}  // namespace Acts
namespace Acts {
namespace protobuf {

// ===================================================================

class MaterialMap_Dimension : public ::google::protobuf::Message /* @@protoc_insertion_point(class_definition:Acts.protobuf.MaterialMap.Dimension) */ {
 public:
  MaterialMap_Dimension();
  virtual ~MaterialMap_Dimension();

  MaterialMap_Dimension(const MaterialMap_Dimension& from);

  inline MaterialMap_Dimension& operator=(const MaterialMap_Dimension& from) {
    CopyFrom(from);
    return *this;
  }
  #if LANG_CXX11
  MaterialMap_Dimension(MaterialMap_Dimension&& from) noexcept
    : MaterialMap_Dimension() {
    *this = ::std::move(from);
  }

  inline MaterialMap_Dimension& operator=(MaterialMap_Dimension&& from) noexcept {
    if (GetArenaNoVirtual() == from.GetArenaNoVirtual()) {
      if (this != &from) InternalSwap(&from);
    } else {
      CopyFrom(from);
    }
    return *this;
  }
  #endif
  static const ::google::protobuf::Descriptor* descriptor();
  static const MaterialMap_Dimension& default_instance();

  static void InitAsDefaultInstance();  // FOR INTERNAL USE ONLY
  static inline const MaterialMap_Dimension* internal_default_instance() {
    return reinterpret_cast<const MaterialMap_Dimension*>(
               &_MaterialMap_Dimension_default_instance_);
  }
  static PROTOBUF_CONSTEXPR int const kIndexInFileMessages =
    0;

  void Swap(MaterialMap_Dimension* other);
  friend void swap(MaterialMap_Dimension& a, MaterialMap_Dimension& b) {
    a.Swap(&b);
  }

  // implements Message ----------------------------------------------

  inline MaterialMap_Dimension* New() const PROTOBUF_FINAL { return New(NULL); }

  MaterialMap_Dimension* New(::google::protobuf::Arena* arena) const PROTOBUF_FINAL;
  void CopyFrom(const ::google::protobuf::Message& from) PROTOBUF_FINAL;
  void MergeFrom(const ::google::protobuf::Message& from) PROTOBUF_FINAL;
  void CopyFrom(const MaterialMap_Dimension& from);
  void MergeFrom(const MaterialMap_Dimension& from);
  void Clear() PROTOBUF_FINAL;
  bool IsInitialized() const PROTOBUF_FINAL;

  size_t ByteSizeLong() const PROTOBUF_FINAL;
  bool MergePartialFromCodedStream(
      ::google::protobuf::io::CodedInputStream* input) PROTOBUF_FINAL;
  void SerializeWithCachedSizes(
      ::google::protobuf::io::CodedOutputStream* output) const PROTOBUF_FINAL;
  ::google::protobuf::uint8* InternalSerializeWithCachedSizesToArray(
      bool deterministic, ::google::protobuf::uint8* target) const PROTOBUF_FINAL;
  int GetCachedSize() const PROTOBUF_FINAL { return _cached_size_; }
  private:
  void SharedCtor();
  void SharedDtor();
  void SetCachedSize(int size) const PROTOBUF_FINAL;
  void InternalSwap(MaterialMap_Dimension* other);
  private:
  inline ::google::protobuf::Arena* GetArenaNoVirtual() const {
    return NULL;
  }
  inline void* MaybeArenaPtr() const {
    return NULL;
  }
  public:

  ::google::protobuf::Metadata GetMetadata() const PROTOBUF_FINAL;

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  // double min = 2;
  void clear_min();
  static const int kMinFieldNumber = 2;
  double min() const;
  void set_min(double value);

  // double max = 3;
  void clear_max();
  static const int kMaxFieldNumber = 3;
  double max() const;
  void set_max(double value);

  // uint32 nBins = 1;
  void clear_nbins();
  static const int kNBinsFieldNumber = 1;
  ::google::protobuf::uint32 nbins() const;
  void set_nbins(::google::protobuf::uint32 value);

  // @@protoc_insertion_point(class_scope:Acts.protobuf.MaterialMap.Dimension)
 private:

  ::google::protobuf::internal::InternalMetadataWithArena _internal_metadata_;
  double min_;
  double max_;
  ::google::protobuf::uint32 nbins_;
  mutable int _cached_size_;
  friend struct ::protobuf_MaterialMap_2eproto::TableStruct;
  friend void ::protobuf_MaterialMap_2eproto::InitDefaultsMaterialMap_DimensionImpl();
};
// -------------------------------------------------------------------

class MaterialMap_MaterialProperties : public ::google::protobuf::Message /* @@protoc_insertion_point(class_definition:Acts.protobuf.MaterialMap.MaterialProperties) */ {
 public:
  MaterialMap_MaterialProperties();
  virtual ~MaterialMap_MaterialProperties();

  MaterialMap_MaterialProperties(const MaterialMap_MaterialProperties& from);

  inline MaterialMap_MaterialProperties& operator=(const MaterialMap_MaterialProperties& from) {
    CopyFrom(from);
    return *this;
  }
  #if LANG_CXX11
  MaterialMap_MaterialProperties(MaterialMap_MaterialProperties&& from) noexcept
    : MaterialMap_MaterialProperties() {
    *this = ::std::move(from);
  }

  inline MaterialMap_MaterialProperties& operator=(MaterialMap_MaterialProperties&& from) noexcept {
    if (GetArenaNoVirtual() == from.GetArenaNoVirtual()) {
      if (this != &from) InternalSwap(&from);
    } else {
      CopyFrom(from);
    }
    return *this;
  }
  #endif
  static const ::google::protobuf::Descriptor* descriptor();
  static const MaterialMap_MaterialProperties& default_instance();

  static void InitAsDefaultInstance();  // FOR INTERNAL USE ONLY
  static inline const MaterialMap_MaterialProperties* internal_default_instance() {
    return reinterpret_cast<const MaterialMap_MaterialProperties*>(
               &_MaterialMap_MaterialProperties_default_instance_);
  }
  static PROTOBUF_CONSTEXPR int const kIndexInFileMessages =
    1;

  void Swap(MaterialMap_MaterialProperties* other);
  friend void swap(MaterialMap_MaterialProperties& a, MaterialMap_MaterialProperties& b) {
    a.Swap(&b);
  }

  // implements Message ----------------------------------------------

  inline MaterialMap_MaterialProperties* New() const PROTOBUF_FINAL { return New(NULL); }

  MaterialMap_MaterialProperties* New(::google::protobuf::Arena* arena) const PROTOBUF_FINAL;
  void CopyFrom(const ::google::protobuf::Message& from) PROTOBUF_FINAL;
  void MergeFrom(const ::google::protobuf::Message& from) PROTOBUF_FINAL;
  void CopyFrom(const MaterialMap_MaterialProperties& from);
  void MergeFrom(const MaterialMap_MaterialProperties& from);
  void Clear() PROTOBUF_FINAL;
  bool IsInitialized() const PROTOBUF_FINAL;

  size_t ByteSizeLong() const PROTOBUF_FINAL;
  bool MergePartialFromCodedStream(
      ::google::protobuf::io::CodedInputStream* input) PROTOBUF_FINAL;
  void SerializeWithCachedSizes(
      ::google::protobuf::io::CodedOutputStream* output) const PROTOBUF_FINAL;
  ::google::protobuf::uint8* InternalSerializeWithCachedSizesToArray(
      bool deterministic, ::google::protobuf::uint8* target) const PROTOBUF_FINAL;
  int GetCachedSize() const PROTOBUF_FINAL { return _cached_size_; }
  private:
  void SharedCtor();
  void SharedDtor();
  void SetCachedSize(int size) const PROTOBUF_FINAL;
  void InternalSwap(MaterialMap_MaterialProperties* other);
  private:
  inline ::google::protobuf::Arena* GetArenaNoVirtual() const {
    return NULL;
  }
  inline void* MaybeArenaPtr() const {
    return NULL;
  }
  public:

  ::google::protobuf::Metadata GetMetadata() const PROTOBUF_FINAL;

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  // double thickness = 1;
  void clear_thickness();
  static const int kThicknessFieldNumber = 1;
  double thickness() const;
  void set_thickness(double value);

  // double X0 = 2;
  void clear_x0();
  static const int kX0FieldNumber = 2;
  double x0() const;
  void set_x0(double value);

  // double L0 = 3;
  void clear_l0();
  static const int kL0FieldNumber = 3;
  double l0() const;
  void set_l0(double value);

  // double A = 4;
  void clear_a();
  static const int kAFieldNumber = 4;
  double a() const;
  void set_a(double value);

  // double Z = 5;
  void clear_z();
  static const int kZFieldNumber = 5;
  double z() const;
  void set_z(double value);

  // double rho = 6;
  void clear_rho();
  static const int kRhoFieldNumber = 6;
  double rho() const;
  void set_rho(double value);

  // @@protoc_insertion_point(class_scope:Acts.protobuf.MaterialMap.MaterialProperties)
 private:

  ::google::protobuf::internal::InternalMetadataWithArena _internal_metadata_;
  double thickness_;
  double x0_;
  double l0_;
  double a_;
  double z_;
  double rho_;
  mutable int _cached_size_;
  friend struct ::protobuf_MaterialMap_2eproto::TableStruct;
  friend void ::protobuf_MaterialMap_2eproto::InitDefaultsMaterialMap_MaterialPropertiesImpl();
};
// -------------------------------------------------------------------

class MaterialMap : public ::google::protobuf::Message /* @@protoc_insertion_point(class_definition:Acts.protobuf.MaterialMap) */ {
 public:
  MaterialMap();
  virtual ~MaterialMap();

  MaterialMap(const MaterialMap& from);

  inline MaterialMap& operator=(const MaterialMap& from) {
    CopyFrom(from);
    return *this;
  }
  #if LANG_CXX11
  MaterialMap(MaterialMap&& from) noexcept
    : MaterialMap() {
    *this = ::std::move(from);
  }

  inline MaterialMap& operator=(MaterialMap&& from) noexcept {
    if (GetArenaNoVirtual() == from.GetArenaNoVirtual()) {
      if (this != &from) InternalSwap(&from);
    } else {
      CopyFrom(from);
    }
    return *this;
  }
  #endif
  static const ::google::protobuf::Descriptor* descriptor();
  static const MaterialMap& default_instance();

  static void InitAsDefaultInstance();  // FOR INTERNAL USE ONLY
  static inline const MaterialMap* internal_default_instance() {
    return reinterpret_cast<const MaterialMap*>(
               &_MaterialMap_default_instance_);
  }
  static PROTOBUF_CONSTEXPR int const kIndexInFileMessages =
    2;

  void Swap(MaterialMap* other);
  friend void swap(MaterialMap& a, MaterialMap& b) {
    a.Swap(&b);
  }

  // implements Message ----------------------------------------------

  inline MaterialMap* New() const PROTOBUF_FINAL { return New(NULL); }

  MaterialMap* New(::google::protobuf::Arena* arena) const PROTOBUF_FINAL;
  void CopyFrom(const ::google::protobuf::Message& from) PROTOBUF_FINAL;
  void MergeFrom(const ::google::protobuf::Message& from) PROTOBUF_FINAL;
  void CopyFrom(const MaterialMap& from);
  void MergeFrom(const MaterialMap& from);
  void Clear() PROTOBUF_FINAL;
  bool IsInitialized() const PROTOBUF_FINAL;

  size_t ByteSizeLong() const PROTOBUF_FINAL;
  bool MergePartialFromCodedStream(
      ::google::protobuf::io::CodedInputStream* input) PROTOBUF_FINAL;
  void SerializeWithCachedSizes(
      ::google::protobuf::io::CodedOutputStream* output) const PROTOBUF_FINAL;
  ::google::protobuf::uint8* InternalSerializeWithCachedSizesToArray(
      bool deterministic, ::google::protobuf::uint8* target) const PROTOBUF_FINAL;
  int GetCachedSize() const PROTOBUF_FINAL { return _cached_size_; }
  private:
  void SharedCtor();
  void SharedDtor();
  void SetCachedSize(int size) const PROTOBUF_FINAL;
  void InternalSwap(MaterialMap* other);
  private:
  inline ::google::protobuf::Arena* GetArenaNoVirtual() const {
    return NULL;
  }
  inline void* MaybeArenaPtr() const {
    return NULL;
  }
  public:

  ::google::protobuf::Metadata GetMetadata() const PROTOBUF_FINAL;

  // nested types ----------------------------------------------------

  typedef MaterialMap_Dimension Dimension;
  typedef MaterialMap_MaterialProperties MaterialProperties;

  // accessors -------------------------------------------------------

  // repeated .Acts.protobuf.MaterialMap.MaterialProperties bins = 8;
  int bins_size() const;
  void clear_bins();
  static const int kBinsFieldNumber = 8;
  const ::Acts::protobuf::MaterialMap_MaterialProperties& bins(int index) const;
  ::Acts::protobuf::MaterialMap_MaterialProperties* mutable_bins(int index);
  ::Acts::protobuf::MaterialMap_MaterialProperties* add_bins();
  ::google::protobuf::RepeatedPtrField< ::Acts::protobuf::MaterialMap_MaterialProperties >*
      mutable_bins();
  const ::google::protobuf::RepeatedPtrField< ::Acts::protobuf::MaterialMap_MaterialProperties >&
      bins() const;

  // .Acts.protobuf.MaterialMap.Dimension l0 = 1;
  bool has_l0() const;
  void clear_l0();
  static const int kL0FieldNumber = 1;
  const ::Acts::protobuf::MaterialMap_Dimension& l0() const;
  ::Acts::protobuf::MaterialMap_Dimension* release_l0();
  ::Acts::protobuf::MaterialMap_Dimension* mutable_l0();
  void set_allocated_l0(::Acts::protobuf::MaterialMap_Dimension* l0);

  // .Acts.protobuf.MaterialMap.Dimension l1 = 2;
  bool has_l1() const;
  void clear_l1();
  static const int kL1FieldNumber = 2;
  const ::Acts::protobuf::MaterialMap_Dimension& l1() const;
  ::Acts::protobuf::MaterialMap_Dimension* release_l1();
  ::Acts::protobuf::MaterialMap_Dimension* mutable_l1();
  void set_allocated_l1(::Acts::protobuf::MaterialMap_Dimension* l1);

  // uint64 geo_id = 3;
  void clear_geo_id();
  static const int kGeoIdFieldNumber = 3;
  ::google::protobuf::uint64 geo_id() const;
  void set_geo_id(::google::protobuf::uint64 value);

  // int32 vol_id = 4;
  void clear_vol_id();
  static const int kVolIdFieldNumber = 4;
  ::google::protobuf::int32 vol_id() const;
  void set_vol_id(::google::protobuf::int32 value);

  // int32 lay_id = 5;
  void clear_lay_id();
  static const int kLayIdFieldNumber = 5;
  ::google::protobuf::int32 lay_id() const;
  void set_lay_id(::google::protobuf::int32 value);

  // int32 app_id = 6;
  void clear_app_id();
  static const int kAppIdFieldNumber = 6;
  ::google::protobuf::int32 app_id() const;
  void set_app_id(::google::protobuf::int32 value);

  // int32 sen_id = 7;
  void clear_sen_id();
  static const int kSenIdFieldNumber = 7;
  ::google::protobuf::int32 sen_id() const;
  void set_sen_id(::google::protobuf::int32 value);

  // @@protoc_insertion_point(class_scope:Acts.protobuf.MaterialMap)
 private:

  ::google::protobuf::internal::InternalMetadataWithArena _internal_metadata_;
  ::google::protobuf::RepeatedPtrField< ::Acts::protobuf::MaterialMap_MaterialProperties > bins_;
  ::Acts::protobuf::MaterialMap_Dimension* l0_;
  ::Acts::protobuf::MaterialMap_Dimension* l1_;
  ::google::protobuf::uint64 geo_id_;
  ::google::protobuf::int32 vol_id_;
  ::google::protobuf::int32 lay_id_;
  ::google::protobuf::int32 app_id_;
  ::google::protobuf::int32 sen_id_;
  mutable int _cached_size_;
  friend struct ::protobuf_MaterialMap_2eproto::TableStruct;
  friend void ::protobuf_MaterialMap_2eproto::InitDefaultsMaterialMapImpl();
};
// ===================================================================


// ===================================================================

#ifdef __GNUC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif  // __GNUC__
// MaterialMap_Dimension

// uint32 nBins = 1;
inline void MaterialMap_Dimension::clear_nbins() {
  nbins_ = 0u;
}
inline ::google::protobuf::uint32 MaterialMap_Dimension::nbins() const {
  // @@protoc_insertion_point(field_get:Acts.protobuf.MaterialMap.Dimension.nBins)
  return nbins_;
}
inline void MaterialMap_Dimension::set_nbins(::google::protobuf::uint32 value) {
  
  nbins_ = value;
  // @@protoc_insertion_point(field_set:Acts.protobuf.MaterialMap.Dimension.nBins)
}

// double min = 2;
inline void MaterialMap_Dimension::clear_min() {
  min_ = 0;
}
inline double MaterialMap_Dimension::min() const {
  // @@protoc_insertion_point(field_get:Acts.protobuf.MaterialMap.Dimension.min)
  return min_;
}
inline void MaterialMap_Dimension::set_min(double value) {
  
  min_ = value;
  // @@protoc_insertion_point(field_set:Acts.protobuf.MaterialMap.Dimension.min)
}

// double max = 3;
inline void MaterialMap_Dimension::clear_max() {
  max_ = 0;
}
inline double MaterialMap_Dimension::max() const {
  // @@protoc_insertion_point(field_get:Acts.protobuf.MaterialMap.Dimension.max)
  return max_;
}
inline void MaterialMap_Dimension::set_max(double value) {
  
  max_ = value;
  // @@protoc_insertion_point(field_set:Acts.protobuf.MaterialMap.Dimension.max)
}

// -------------------------------------------------------------------

// MaterialMap_MaterialProperties

// double thickness = 1;
inline void MaterialMap_MaterialProperties::clear_thickness() {
  thickness_ = 0;
}
inline double MaterialMap_MaterialProperties::thickness() const {
  // @@protoc_insertion_point(field_get:Acts.protobuf.MaterialMap.MaterialProperties.thickness)
  return thickness_;
}
inline void MaterialMap_MaterialProperties::set_thickness(double value) {
  
  thickness_ = value;
  // @@protoc_insertion_point(field_set:Acts.protobuf.MaterialMap.MaterialProperties.thickness)
}

// double X0 = 2;
inline void MaterialMap_MaterialProperties::clear_x0() {
  x0_ = 0;
}
inline double MaterialMap_MaterialProperties::x0() const {
  // @@protoc_insertion_point(field_get:Acts.protobuf.MaterialMap.MaterialProperties.X0)
  return x0_;
}
inline void MaterialMap_MaterialProperties::set_x0(double value) {
  
  x0_ = value;
  // @@protoc_insertion_point(field_set:Acts.protobuf.MaterialMap.MaterialProperties.X0)
}

// double L0 = 3;
inline void MaterialMap_MaterialProperties::clear_l0() {
  l0_ = 0;
}
inline double MaterialMap_MaterialProperties::l0() const {
  // @@protoc_insertion_point(field_get:Acts.protobuf.MaterialMap.MaterialProperties.L0)
  return l0_;
}
inline void MaterialMap_MaterialProperties::set_l0(double value) {
  
  l0_ = value;
  // @@protoc_insertion_point(field_set:Acts.protobuf.MaterialMap.MaterialProperties.L0)
}

// double A = 4;
inline void MaterialMap_MaterialProperties::clear_a() {
  a_ = 0;
}
inline double MaterialMap_MaterialProperties::a() const {
  // @@protoc_insertion_point(field_get:Acts.protobuf.MaterialMap.MaterialProperties.A)
  return a_;
}
inline void MaterialMap_MaterialProperties::set_a(double value) {
  
  a_ = value;
  // @@protoc_insertion_point(field_set:Acts.protobuf.MaterialMap.MaterialProperties.A)
}

// double Z = 5;
inline void MaterialMap_MaterialProperties::clear_z() {
  z_ = 0;
}
inline double MaterialMap_MaterialProperties::z() const {
  // @@protoc_insertion_point(field_get:Acts.protobuf.MaterialMap.MaterialProperties.Z)
  return z_;
}
inline void MaterialMap_MaterialProperties::set_z(double value) {
  
  z_ = value;
  // @@protoc_insertion_point(field_set:Acts.protobuf.MaterialMap.MaterialProperties.Z)
}

// double rho = 6;
inline void MaterialMap_MaterialProperties::clear_rho() {
  rho_ = 0;
}
inline double MaterialMap_MaterialProperties::rho() const {
  // @@protoc_insertion_point(field_get:Acts.protobuf.MaterialMap.MaterialProperties.rho)
  return rho_;
}
inline void MaterialMap_MaterialProperties::set_rho(double value) {
  
  rho_ = value;
  // @@protoc_insertion_point(field_set:Acts.protobuf.MaterialMap.MaterialProperties.rho)
}

// -------------------------------------------------------------------

// MaterialMap

// .Acts.protobuf.MaterialMap.Dimension l0 = 1;
inline bool MaterialMap::has_l0() const {
  return this != internal_default_instance() && l0_ != NULL;
}
inline void MaterialMap::clear_l0() {
  if (GetArenaNoVirtual() == NULL && l0_ != NULL) {
    delete l0_;
  }
  l0_ = NULL;
}
inline const ::Acts::protobuf::MaterialMap_Dimension& MaterialMap::l0() const {
  const ::Acts::protobuf::MaterialMap_Dimension* p = l0_;
  // @@protoc_insertion_point(field_get:Acts.protobuf.MaterialMap.l0)
  return p != NULL ? *p : *reinterpret_cast<const ::Acts::protobuf::MaterialMap_Dimension*>(
      &::Acts::protobuf::_MaterialMap_Dimension_default_instance_);
}
inline ::Acts::protobuf::MaterialMap_Dimension* MaterialMap::release_l0() {
  // @@protoc_insertion_point(field_release:Acts.protobuf.MaterialMap.l0)
  
  ::Acts::protobuf::MaterialMap_Dimension* temp = l0_;
  l0_ = NULL;
  return temp;
}
inline ::Acts::protobuf::MaterialMap_Dimension* MaterialMap::mutable_l0() {
  
  if (l0_ == NULL) {
    l0_ = new ::Acts::protobuf::MaterialMap_Dimension;
  }
  // @@protoc_insertion_point(field_mutable:Acts.protobuf.MaterialMap.l0)
  return l0_;
}
inline void MaterialMap::set_allocated_l0(::Acts::protobuf::MaterialMap_Dimension* l0) {
  ::google::protobuf::Arena* message_arena = GetArenaNoVirtual();
  if (message_arena == NULL) {
    delete l0_;
  }
  if (l0) {
    ::google::protobuf::Arena* submessage_arena = NULL;
    if (message_arena != submessage_arena) {
      l0 = ::google::protobuf::internal::GetOwnedMessage(
          message_arena, l0, submessage_arena);
    }
    
  } else {
    
  }
  l0_ = l0;
  // @@protoc_insertion_point(field_set_allocated:Acts.protobuf.MaterialMap.l0)
}

// .Acts.protobuf.MaterialMap.Dimension l1 = 2;
inline bool MaterialMap::has_l1() const {
  return this != internal_default_instance() && l1_ != NULL;
}
inline void MaterialMap::clear_l1() {
  if (GetArenaNoVirtual() == NULL && l1_ != NULL) {
    delete l1_;
  }
  l1_ = NULL;
}
inline const ::Acts::protobuf::MaterialMap_Dimension& MaterialMap::l1() const {
  const ::Acts::protobuf::MaterialMap_Dimension* p = l1_;
  // @@protoc_insertion_point(field_get:Acts.protobuf.MaterialMap.l1)
  return p != NULL ? *p : *reinterpret_cast<const ::Acts::protobuf::MaterialMap_Dimension*>(
      &::Acts::protobuf::_MaterialMap_Dimension_default_instance_);
}
inline ::Acts::protobuf::MaterialMap_Dimension* MaterialMap::release_l1() {
  // @@protoc_insertion_point(field_release:Acts.protobuf.MaterialMap.l1)
  
  ::Acts::protobuf::MaterialMap_Dimension* temp = l1_;
  l1_ = NULL;
  return temp;
}
inline ::Acts::protobuf::MaterialMap_Dimension* MaterialMap::mutable_l1() {
  
  if (l1_ == NULL) {
    l1_ = new ::Acts::protobuf::MaterialMap_Dimension;
  }
  // @@protoc_insertion_point(field_mutable:Acts.protobuf.MaterialMap.l1)
  return l1_;
}
inline void MaterialMap::set_allocated_l1(::Acts::protobuf::MaterialMap_Dimension* l1) {
  ::google::protobuf::Arena* message_arena = GetArenaNoVirtual();
  if (message_arena == NULL) {
    delete l1_;
  }
  if (l1) {
    ::google::protobuf::Arena* submessage_arena = NULL;
    if (message_arena != submessage_arena) {
      l1 = ::google::protobuf::internal::GetOwnedMessage(
          message_arena, l1, submessage_arena);
    }
    
  } else {
    
  }
  l1_ = l1;
  // @@protoc_insertion_point(field_set_allocated:Acts.protobuf.MaterialMap.l1)
}

// uint64 geo_id = 3;
inline void MaterialMap::clear_geo_id() {
  geo_id_ = GOOGLE_ULONGLONG(0);
}
inline ::google::protobuf::uint64 MaterialMap::geo_id() const {
  // @@protoc_insertion_point(field_get:Acts.protobuf.MaterialMap.geo_id)
  return geo_id_;
}
inline void MaterialMap::set_geo_id(::google::protobuf::uint64 value) {
  
  geo_id_ = value;
  // @@protoc_insertion_point(field_set:Acts.protobuf.MaterialMap.geo_id)
}

// int32 vol_id = 4;
inline void MaterialMap::clear_vol_id() {
  vol_id_ = 0;
}
inline ::google::protobuf::int32 MaterialMap::vol_id() const {
  // @@protoc_insertion_point(field_get:Acts.protobuf.MaterialMap.vol_id)
  return vol_id_;
}
inline void MaterialMap::set_vol_id(::google::protobuf::int32 value) {
  
  vol_id_ = value;
  // @@protoc_insertion_point(field_set:Acts.protobuf.MaterialMap.vol_id)
}

// int32 lay_id = 5;
inline void MaterialMap::clear_lay_id() {
  lay_id_ = 0;
}
inline ::google::protobuf::int32 MaterialMap::lay_id() const {
  // @@protoc_insertion_point(field_get:Acts.protobuf.MaterialMap.lay_id)
  return lay_id_;
}
inline void MaterialMap::set_lay_id(::google::protobuf::int32 value) {
  
  lay_id_ = value;
  // @@protoc_insertion_point(field_set:Acts.protobuf.MaterialMap.lay_id)
}

// int32 app_id = 6;
inline void MaterialMap::clear_app_id() {
  app_id_ = 0;
}
inline ::google::protobuf::int32 MaterialMap::app_id() const {
  // @@protoc_insertion_point(field_get:Acts.protobuf.MaterialMap.app_id)
  return app_id_;
}
inline void MaterialMap::set_app_id(::google::protobuf::int32 value) {
  
  app_id_ = value;
  // @@protoc_insertion_point(field_set:Acts.protobuf.MaterialMap.app_id)
}

// int32 sen_id = 7;
inline void MaterialMap::clear_sen_id() {
  sen_id_ = 0;
}
inline ::google::protobuf::int32 MaterialMap::sen_id() const {
  // @@protoc_insertion_point(field_get:Acts.protobuf.MaterialMap.sen_id)
  return sen_id_;
}
inline void MaterialMap::set_sen_id(::google::protobuf::int32 value) {
  
  sen_id_ = value;
  // @@protoc_insertion_point(field_set:Acts.protobuf.MaterialMap.sen_id)
}

// repeated .Acts.protobuf.MaterialMap.MaterialProperties bins = 8;
inline int MaterialMap::bins_size() const {
  return bins_.size();
}
inline void MaterialMap::clear_bins() {
  bins_.Clear();
}
inline const ::Acts::protobuf::MaterialMap_MaterialProperties& MaterialMap::bins(int index) const {
  // @@protoc_insertion_point(field_get:Acts.protobuf.MaterialMap.bins)
  return bins_.Get(index);
}
inline ::Acts::protobuf::MaterialMap_MaterialProperties* MaterialMap::mutable_bins(int index) {
  // @@protoc_insertion_point(field_mutable:Acts.protobuf.MaterialMap.bins)
  return bins_.Mutable(index);
}
inline ::Acts::protobuf::MaterialMap_MaterialProperties* MaterialMap::add_bins() {
  // @@protoc_insertion_point(field_add:Acts.protobuf.MaterialMap.bins)
  return bins_.Add();
}
inline ::google::protobuf::RepeatedPtrField< ::Acts::protobuf::MaterialMap_MaterialProperties >*
MaterialMap::mutable_bins() {
  // @@protoc_insertion_point(field_mutable_list:Acts.protobuf.MaterialMap.bins)
  return &bins_;
}
inline const ::google::protobuf::RepeatedPtrField< ::Acts::protobuf::MaterialMap_MaterialProperties >&
MaterialMap::bins() const {
  // @@protoc_insertion_point(field_list:Acts.protobuf.MaterialMap.bins)
  return bins_;
}

#ifdef __GNUC__
  #pragma GCC diagnostic pop
#endif  // __GNUC__
// -------------------------------------------------------------------

// -------------------------------------------------------------------


// @@protoc_insertion_point(namespace_scope)

}  // namespace protobuf
}  // namespace Acts

// @@protoc_insertion_point(global_scope)

#endif  // PROTOBUF_MaterialMap_2eproto__INCLUDED
