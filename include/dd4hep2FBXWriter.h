/**************************************************************************
 * BASF2 (Belle Analysis Framework 2)                                     *
 * Copyright(C) 2016 - Belle II Collaboration                             *
 *                                                                        *
 * Author: Leo Piilonen                                                   *
 *                                                                        *
 * This software is provided "as is" without any warranty.                *
 **************************************************************************/

#ifndef DD4HEP2FBXWRITER_H
#define DD4HEP2FBXWRITER_H
#include "HepPolyhedron.h"
#include "GPolyhedron.hh"
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Detector.h"
#include <DD4hep/Handle.h>
#include "XML/Layering.h"
#include "XML/XML.h"
#include "TColor.h"

// #include "Math/Vector3D.h"
#include "Math/Transform3D.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;
typedef std::map<std::string, Handle<NamedObject>> HandleMap;
typedef std::map<std::string, DetElement> Children;


/** The FBX-writer module.
 *
 * This module goes through all volumes of the constructed GEANT4
 * geometry and writes an Autodesk FBX (text) file.
 *
 * Prerequisite: This module requires a valid GEANT4 geometry.
 *
 */
class dd4hep2FBXWriter {

public:

  //! Constructor [empty]
  dd4hep2FBXWriter() {}
  dd4hep2FBXWriter(string filePath, bool usePrototypes);

  //! Destructor [empty]
  ~dd4hep2FBXWriter() {}

  //! Write the geometry in the FBX (text) format
  //! @param outputFilename: user-specified output filename [empty string -> "geometry.fbx"]
  //! @param usePrototypes: true => write LogVol and PhysVol prototypes; false => write each instance of LogVol and PhysVol
  //! @return true if the FBX file was written; false otherwise
  bool doit(std::string outputFilename);

private:

  //! Create unique and legal name for each solid, logical volume, physical volume
  // void assignName(std::vector<std::string>*, unsigned int, const G4String&, int);
  std::vector<std::string> assignName(std::vector<std::string> names, string originalName, unsigned int mindex);
  // get all children detectors of det, include children's children
  void getDets(DetElement det);

  // //! Write FBX definition for each solid's polyhedron
  // void writeGeometryNode(G4VSolid*, const std::string&, unsigned long long);
  void writeGeometryNode(Solid solid, const std::string solidName, unsigned long long solidID);

  // //! Write FBX definition for each logical volume's color information
  // void writeMaterialNode(int, const std::string&);
  void writeMaterialNode(int lvIndex, const std::string matName);

  // //! Write FBX definition for each logical volume
  // void writeLVModelNode(G4LogicalVolume*, const std::string&, unsigned long long);
  void writeLVModelNode(TGeoVolume *Vol, const std::string lvName, unsigned long long lvID);


  // //! Write FBX definition for each physical volume
  // void writePVModelNode(G4VPhysicalVolume*, const std::string&, unsigned long long);
  void writePVModelNode(TGeoNode *node, const std::string pvName, unsigned long long pvID);

  // //! Count the physical volumes, logical volumes, materials and solids (recursive)
  // void countEntities(G4VPhysicalVolume*);
  // void countEntities(DetElement world)

  // //! Process one physical volume for FBX-node writing (recursive)
  // void addModels(G4VPhysicalVolume*, int);
  void addModels(TGeoNode *node, int replica, unsigned long long pvIndex);

  // //! Write FBX connections among all of the nodes in the tree (recursive)
  void addConnections(DetElement, int);

  //! Write FBX at the start of the file
  void writePreamble(int, int, int);

  // //! Write FBX definition for the solid's polyhedron
  // void writePolyhedron(Solid, GPolyhedron*, const std::string&, unsigned long long);
  void writePolyhedron(Solid solid, GPolyhedron* polyhedron, const std::string name,
                                      unsigned long long solidID);

  // //! Write FBX connection for each logical volume's solid and color info
  void writeSolidToLV(const std::string, const std::string, bool, unsigned long long, unsigned long long, unsigned long long);

  // //! Write FBX connection for each physical volume's solid and color info (bypass singleton logical volume)
  void writeSolidToPV(const std::string, const std::string, bool, unsigned long long, unsigned long long, unsigned long long);

  // //! Write FBX connection for the (unique) logical volume of a physical volume
  void writeLVToPV(const std::string, const std::string, unsigned long long, unsigned long long);

  // no //! Write FBX connection for each physical-volume daughter of a parent logical volume
  void writePVToParentLV(const std::string, const std::string, unsigned long long, unsigned long long);

  // //! Write FBX connection for each physical-volume daughter of a parent physical volume (bypass singleton logical volume)
  void writePVToParentPV(const std::string, const std::string, unsigned long long, unsigned long long);

  // //! Create polyhedron for a boolean solid (recursive)
  // HepPolyhedron* getBooleanSolidPolyhedron(G4VSolid*);
  HepPolyhedron* getBooleanSolidPolyhedron(Solid solid);



void getVolSolid(TGeoNode *node);


/////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Detector mm_lcdd;
  string  m_filePath;
  // std::vector<DetElement> m_det;
  DetElement m_world;
/////////////////////////////////////////////////////////////////////////////////////////////////////////
  //! User-specified flag to select whether to write and re-use logical- and physical-volume
  //! prototypes once (true) or to write duplicate instances of each such volume (false).
  bool m_UsePrototypes;

  //! Output file
  std::ofstream m_File;

  //! Modified (legal-character and unique) physical-volume name
  std::vector<std::string> m_detName;

  //! Modified (legal-character and unique) logical-volume name
  std::vector<std::string> m_volName;

  //! Modified (legal-character and unique) solid name
  std::vector<std::string> m_solidName;

  std::vector<std::string> m_placedVolName;
  //! Modified (legal-character and unique) physical-volume 
  std::vector<DetElement> m_det;

  //! Modified (legal-character and unique) logical-volume 
  // std::vector<Volume> m_childrenVol;
  std::vector<TGeoVolume*> m_vol;
  std::vector<TGeoShape*> m_solid;
  std::vector<TGeoNode*> m_placedVol;
  std::vector<double> m_volChildrenNum;
  // std::vector<TGeoVolume*> m_Vol;
  //! Flag to indicate that the logical volume is visible
  std::vector<bool> m_visible;
  std::vector<bool> m_assembly;
  std::vector<TColor*> m_color;


  //! Modified (legal-character and unique) solid 
  // std::vector<Solid> m_childrenSolid;


  //! Unique identifiers for physical volumes (Model nodes with transformation information)
  std::vector<unsigned long long>* m_PVID;

  //! Unique identifiers for logical volumes (Model nodes with links to Geometry and Material)
  std::vector<unsigned long long>* m_LVID;

  //! Unique identifiers for logical volumes' color information (Material nodes)
  std::vector<unsigned long long>* m_MatID;

  //! Unique identifiers for solids (Geometry nodes)
  std::vector<unsigned long long>* m_SolidID;


  //! Count of number of instances of each physical volume
  std::vector<unsigned int>* m_PVCount;

  //! Count of number of instances of each logical volume
  std::vector<unsigned int>* m_LVCount;

  //! Count of number of instances of each solid (typically 1)
  std::vector<unsigned int>* m_SolidCount;

  //! Count of number of replicas of each replicated physical volume
  std::vector<unsigned int>* m_PVReplicas;

  //! Count of number of replicas of each logical volume associated with a replicated physical volume
  std::vector<unsigned int>* m_LVReplicas;

  //! Count of number of replicas of each solid (extras for replicas with modified solids)
  std::vector<unsigned int>* m_SolidReplicas;

  //! Flag to indicate that a logical volume is referenced at most once (eligible for bypass)
    std::vector<bool>* m_LVUnique;
  
  //! logical Vol = phyVol + phyVol + phyVol + .... , get the num of logVol children
  std::vector<double> m_VolChildrenNum;


};

#endif
