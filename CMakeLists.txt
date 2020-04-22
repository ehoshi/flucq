cmake_minimum_required(VERSION 2.8.12)

project(flucq Fortran)
add_subdirectory(support EXCLUDE_FROM_ALL)

add_executable(flucq
  qdyn.f
  Jstuff.f
  dgesv.f
  ewaldstuff.f
  geniostuff.f
  iostuff.f
  mathstuff.f
  qdynmd.f
  qdynstuff.f
  qpotstuff.f
  ranstuff.f)
set_property(TARGET flucq PROPERTY OUTPUT_NAME_DEBUG dflucq)

add_executable(props
  propsm.f
  Jstuff.f
  dgesv.f
  ewaldstuff.f
  geniostuff.f
  iostuff.f
  mathstuff.f
  qdynmd.f
  qdynstuff.f
  qpotstuff.f
  ranstuff.f)

set_target_properties(flucq props PROPERTIES Fortran_FORMAT FIXED)
