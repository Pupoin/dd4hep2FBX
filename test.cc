/**************************************************************************
 * BASF2 (Belle Analysis Framework 2)                                     *
 * Copyright(C) 2016 - Belle II Collaboration                             *
 *                                                                        *
 * Author: ZhaoYang Yuan                                                   *
 * 2024/02/20                                                                        *
 * This software is provided "as is" without any warranty.                *
 **************************************************************************/

#include "dd4hep2FBXWriter.h"
#include <CLHEP/Geometry/Transform3D.h>

#include <iomanip>


int main()
{
  // Detector &lcdd = Detector::getInstance();
  string filePath="/home/wln/DD4hep_source/DDDetectors/compact/SiD.xml";
  // string filePath="/home/wln/DD4hep/DDJunoDetectors/compact/Juno.xml";
  // Tube waterPool(0, 10 * mm, 10 * mm, 20 * mm, 30 * mm);
  // lcdd.fromCompact("/home/wln/DD4hep_source/DDDetectors/compact/SiD.xml");
  // const HandleMap det_map = lcdd.detectors();

  // dd4hep2FBXWriter test(det_map, true);
  dd4hep2FBXWriter test(filePath, true);
  test.doit("test.fbx");


  std::cout << "----------------- test ----------------------" << std::endl;
  ////** do no comment me!!!!!!!!!!!!!!!!!!!!!!!!!1
  HepGeom::Point3D<double> base(1, 2, 3);
  HepGeom::Vector3D<double> base2(1, 2, 3);
  HepGeom::Transform3D m;
  HepGeom::RotateX3D rz = HepGeom::RotateX3D(1.2);
  HepGeom::RotateX3D rx1 = HepGeom::RotateX3D(1.3);
  HepGeom::RotateX3D rz2 = HepGeom::RotateX3D(1.5);
  //  HepGeom::Vector3D<double>  out=rz2*(rx1*(rz*base));
  std::cout << "||||||" << std::endl;

  HepGeom::Point3D<double> out = rx1 * base;
  HepGeom::Vector3D<double> out2 = rx1 * base2;

  std::cout << out.x() << " " << out.y() << " " << out.z() << std::endl;


  return 0 ;
}