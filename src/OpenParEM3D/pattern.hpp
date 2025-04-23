////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//    OpenParEM3D - A fullwave 3D electromagnetic simulator.                  //
//    Copyright (C) 2025 Brian Young                                          //
//                                                                            //
//    This program is free software: you can redistribute it and/or modify    //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation, either version 3 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    This program is distributed in the hope that it will be useful,         //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef PATTERN_H
#define PATTERN_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <limits>
#include <filesystem>
#include <unistd.h>
#include <vector>
#include <cmath>
#include <climits>
#include <cfloat>
#include <complex>
#include "mpi.h"
#include "mfem.hpp"
#include "petscsys.h"
#include "project.h"
#include "misc.hpp"
#include "keywordPair.hpp"

using namespace std;
using namespace mfem;

double signed_angle_between_two_lines (struct point, struct point, struct point, struct point);

bool double_compare (double, double, double);

class OPEMpoint
{
   private:
      double x,y,z;
      double theta;                            // spherical coordinates: angle from the z-axis to this point
      double phi;                              // spherical coordinates: angle from the x-axis in the x-y plane to this point
      complex<double> Etheta,Ephi,Htheta,Hphi; // computed results 
      double gain,directivity;                 // post-processed results, convenience for writing a raw file
   public:
      OPEMpoint ();
      OPEMpoint (double, double, double);
      OPEMpoint (OPEMpoint *);
      void copy (OPEMpoint *);
      void set_xyz (double x_, double y_, double z_) {x=x_; y=y_; z=z_;}
      double get_x () {return x;}
      double get_y () {return y;}
      double get_z () {return z;}
      void add (OPEMpoint *);
      void subtract (OPEMpoint *);
      void scale (double);
      double dot (OPEMpoint *);
      OPEMpoint* cross (OPEMpoint *);
      double get_angle (OPEMpoint *);
      double get_angle (OPEMpoint *, OPEMpoint *);
      double get_angle (OPEMpoint, OPEMpoint);
      double get_signed_angle (OPEMpoint, OPEMpoint, OPEMpoint *);
      OPEMpoint* midpoint (OPEMpoint *);
      void midpoint (OPEMpoint *, OPEMpoint *);
      bool is_match (OPEMpoint *);
      double length ();
      void normalize (double);
      double distance (OPEMpoint *);
      void set_theta (OPEMpoint *);
      void set_phi (OPEMpoint *);
      double get_theta () {return theta;}
      double get_phi () {return phi;}
      void set_gain (double gain_) {gain=gain_;}
      void set_directivity (double directivity_) {directivity=directivity_;}
      void addMeshVertex (Mesh *);
      void add_Etheta (complex<double> Etheta_) {Etheta+=Etheta_;}
      void add_Ephi (complex<double> Ephi_) {Ephi+=Ephi_;}
      void add_Htheta (complex<double> Htheta_) {Htheta+=Htheta_;}
      void add_Hphi (complex<double> Hphi_) {Hphi+=Hphi_;}
      void set_fields (complex<double> Ephi_, complex<double> Etheta_, complex<double> Hphi_, complex<double> Htheta_)
         {Etheta=Etheta_; Ephi=Ephi_; Htheta=Htheta_; Hphi=Hphi_;}
      complex<double> get_fieldValue (string);
      void sendFieldsTo (int);
      void recvFieldsFrom (int);
      complex<double> get_power() {return 0.5*(Etheta*conj(Hphi)-Ephi*conj(Htheta));}
      void print ();
      void saveHeader (ostream *, double);
      void save (ostream *, double);
};

class OPEMtriangle
{
   private:
      long unsigned int a,b,c;     // indices to pointList that form the triangle
      long unsigned int ab,bc,ca;  // indices to triangleList indicating neighboring triangles 
      long unsigned int center;    // index to pointList
      double area;                 // the triangle's area
   public:
      OPEMtriangle () {};
      OPEMtriangle (long unsigned int, long unsigned int, long unsigned int);
      void copy (OPEMtriangle *);
      bool has_point_index (long unsigned int);
      OPEMpoint* getCenterOPEMpoint (vector<OPEMpoint *> *, OPEMpoint *, double);
      void set_a (long unsigned int a_) {a=a_;}
      void set_b (long unsigned int b_) {b=b_;}
      void set_c (long unsigned int c_) {c=c_;}
      long unsigned int get_a () {return a;}
      long unsigned int get_b () {return b;}
      long unsigned int get_c () {return c;}
      long unsigned int get_ab () {return ab;}
      long unsigned int get_bc () {return bc;}
      long unsigned int get_ca () {return ca;}
      long unsigned int get_center () {return center;}
      double get_area () {return area;}
      //double get_areaPlanar (vector<OPEMpoint *> *);
      void calculateArea (vector<OPEMpoint *> *, OPEMpoint*, double);
      void set_metrics (OPEMpoint *, double, vector<OPEMpoint *> *);
      bool is_matched ();
      bool has_shared_edge (OPEMtriangle *, long unsigned int, long unsigned int);
      void print ();
      double get_theta (vector<OPEMpoint *> *);
      double get_phi (vector<OPEMpoint *> *);
      void addMeshTriangle (Mesh *);
};

class Sphere
{
   private:
      double radius;
      OPEMpoint center;
      vector<OPEMpoint *> pointList;             // points making up the sphere
      vector<OPEMtriangle *> triangleList;       // triangles formed by the points
      vector<double> areaList;                   // areas of the sphere allocated to each point
      double angularResolution;
      bool hasSaved;                             // convenience data for saving
   public:
      Sphere ();
      Sphere* clone ();
      bool is_match (Sphere *);
      void set_radius (double radius_) {radius=radius_;}
      void set_center (double x, double y, double z) {center.set_xyz(x,y,z);}
      void pushOPEMpoint (double, double, double);
      void pushOPEMtriangle (long unsigned int, long unsigned int, long unsigned int);
      long unsigned int getTriangleCount () {return triangleList.size();}
      void create ();
      long unsigned int get_index (OPEMpoint *);
      void refine (long unsigned int);
      void findNeighbors ();
      void checkNeighbors ();
      void allocateAreasToPoints ();
      void setMetrics ();
      double getTotalArea ();
      double getMaxAbsValue (string);
      complex<double> getMaxValue (string);
      double get_angularResolution () {return angularResolution;}
      complex<double> getRadiatedPower ();
      void createPatternMesh (struct projectData *, double, complex<double>, complex<double>, string, double, int);
      vector<OPEMpoint *>* get_pointList() {return &pointList;}
      vector<OPEMtriangle *>* get_triangleList () {return &triangleList;}
      void calculateAngularResolution ();
      double calculateIsotropicGain (complex<double>, double);
      double calculateDirectivity (complex<double>, double);
      void print ();
      bool get_hasSaved () {return hasSaved;}
      void set_hasSaved (bool hasSaved_) {hasSaved=hasSaved_;}
      void save (ostream *);
      ~Sphere ();
};

class Circle
{
   private:
      double radius;
      OPEMpoint center;
      string plane;
      double theta,phi,latitude,rotation;
      int nPoints;                         // number of points on the circle
      vector<struct point> xyzList;        // unrotated xyz values of the circle (convenience data)
      vector<OPEMpoint *> pointList;
      double angularResolution;
      bool hasSaved;                       // convenience data for saving
   public:
      Circle ();
      bool has_named_plane ();
      Circle* clone ();
      void set_radius (double radius_) {radius=radius_;}
      void set_center (double x, double y, double z) {center.set_xyz(x,y,z);}
      void set_plane (string plane_) {plane=plane_;}
      void set_angles (double theta_, double phi_, double latitude_, double rotation_) {
         theta=theta_*M_PI/180; phi=phi_*M_PI/180; latitude=latitude_*M_PI/180; rotation=rotation_*M_PI/180;
      }
      void set_nPoints (int nPoints_) {nPoints=nPoints_;}
      vector<OPEMpoint *>* get_pointList() {return &pointList;}
      void create ();
      bool is_match (Circle *);
      void createPatternMesh (struct projectData *, double, complex<double>, complex<double>, Sphere *, string, double, int);
      bool createReport (struct projectData *, string, string, complex<double>, complex<double>,
                         Sphere *, double, int, double, double, double, double);
      void calculateAngularResolution ();
      double calculateIsotropicGain (complex<double>, double);
      double calculateDirectivity (complex<double>, double);
      void print ();
      void print (PetscMPIInt);
      bool get_hasSaved () {return hasSaved;}
      void set_hasSaved (bool hasSaved_) {hasSaved=hasSaved_;}
      void save (ostream *);
      ~Circle ();
};

class Current
{
   private:
      double x,y,z;
      double area;
      complex<double> Jx,Jy,Jz,Mx,My,Mz;
   public:
      Current () {};
      Current (double, double, double, complex<double>, complex<double>, complex<double>,
               complex<double>, complex<double>, complex<double>, double);
      double length ();
      double dotproduct (double, double, double);
      double get_area () {return area;}
      double get_x () {return x;}
      double get_y () {return y;}
      double get_z () {return z;}
      complex<double> get_Jx () {return Jx;}
      complex<double> get_Jy () {return Jy;}
      complex<double> get_Jz () {return Jz;}
      complex<double> get_Mx () {return Mx;}
      complex<double> get_My () {return My;}
      complex<double> get_Mz () {return Mz;}
      Current* clone ();
      void sendTo (int);
      void recvFrom (int);
};

class Pattern
{
   private:
      bool active;                  // indicates the latest result from iterative refinement
      int iteration;                // iteration number in adaptive refinement
      double frequency;
      int Sport;                    // driving Sport

      int dim;                      // dim=2 for 2D slice plots; dim=3 for 3D volume plots
      string quantity1;             // required quantity to plot

      // additional data for 3D plots
      double totalArea;
      complex<double> radiatedPower;
      complex<double> acceptedPower;
      double gain;                  // peak gain on the circle for dim=2 or the sphere for dim=3
      double directivity;           // peak directivity on the circle for dim=2 or the sphere for dim=3
      double radiationEfficiency;   // peak radiation efficiency on the circle for dim=2 or the sphere for dim=3

      // additional variables for 2D slice plots
      string quantity2;             // optional 2nd quantity to plot for 2D plots
      string plane;                 // keyword to set phi, theta, and spin_variable for common planes
      double theta;                 // rotation of the x-y plane from the z-axis (theta from spherical coordinates)
      double phi;                   // rotation of the x-y plane from the x-axis (phi from spherical coordinates)
      double rotation;              // rotation of the printable report plot to align axes to user preference if the default alignments are not preferred
      double latitude;              // offset in the direction normal to the plane given as latitude on a sphere

      Circle *circle;               // 2D pattern data for dim=2; nullptr for dim=3
      Sphere *sphere;               // 3D pattern data for dim=3; convenience pointer for dim=2

      double equalMagLimit=1e-12;
   public:
      Pattern (double, int, int, complex<double>, int, string, string, string, double, double, double, double, Circle *, Sphere *);
      void copy_additional_data (Pattern *);
      Pattern* clone ();
      bool is_2D () {if (dim == 2) return true; return false;}
      bool is_3D () {if (dim == 3) return true; return false;}
      bool is_match (Circle *);
      bool is_match (Sphere *);
      bool is_match (double, string);
      bool is_match (double, int);
      bool is_match (double, int , int);
      bool is_match (double, int , int, string);
      bool is_match (double, int , string);
      bool is_prior_version (Pattern *);
      bool is_active () {return active;}
      void set_inactive () {active=false;}
      int get_iteration () {return iteration;}
      int get_Sport () {return Sport;}
      int get_dim () {return dim;}
      string get_quantity1 () {return quantity1;}
      string get_quantity2 () {return quantity2;}
      double get_frequency () {return frequency;}
      double get_totalArea () {return totalArea;}
      double get_gain () {return gain;}
      double get_directivity () {return directivity;}
      double get_radiationEfficiency () {return radiationEfficiency;}
      void set_gain (double gain_) {gain=gain_;}
      void set_directivity (double directivity_) {directivity=directivity_;}
      void set_radiationEfficiency (double radiationEfficiency_) {radiationEfficiency=radiationEfficiency_;}
      void calculateTotalArea ();
      void calculateRadiatedPower ();
      void calculateIsotropicGain ();
      void calculateDirectivity ();
      void calculateRadiationEfficiency ();
      complex<double> get_radiatedPower () {return radiatedPower;}
      void set_acceptedPower (complex<double> acceptedPower_) {acceptedPower=acceptedPower_;}
      Circle* get_circle () {return circle;}
      Sphere * get_sphere () {return sphere;}
      void set_circle (Circle *circle_) {circle=circle_;}
      void set_sphere (Sphere *sphere_) {sphere=sphere_;}
      void save3DParaView (struct projectData *, double, int);
      void save2DParaView (struct projectData *, double, int);
      void copyStats (Pattern *);
      void print (PetscMPIInt);
      void set_hasSavedSphere (bool a) {if (sphere) sphere->set_hasSaved(a);}
      void set_hasSavedCircle (bool a) {if (circle) circle->set_hasSaved(a);}
      bool get_hasSavedSphere () {if (sphere) return sphere->get_hasSaved(); else return true;}
      void save (ostream *);
      void save_as_test (ofstream *, const char *, int, int *);
      ~Pattern ();
};

class PatternDatabase
{
   private:
      vector<Pattern *> patternList;
      vector<double> unique_frequencies;  // sorted in increasing order
      double tol=1e-12;
   public:
      void addPattern (double, int, int, complex<double>, int, string, string, string, double, double, double, double, Circle *, Sphere *);
      Pattern* get_Pattern (double, int);
      Circle* get_circle (double, int, int, Circle *);
      Sphere* get_sphere (double, int, int, Sphere *);
      void saveParaView (struct projectData *, double, int);
      void calculateTotalArea (double, int);
      void calculateRadiatedPower (double, int);
      void calculateIsotropicGain (double, int);
      void calculateDirectivity (double, int);
      void calculateRadiationEfficiency (double, int);
      void populate_sphere (double, int, int);
      double get_gain (double, int, int);
      double get_directivity (double, int, int);
      double get_radiationEfficiency (double, int, int);
      bool saveCSV (struct projectData *, vector<double> *, int, bool);
      bool loadCSV (string, struct projectData *);
      void reset ();
      void print (PetscMPIInt);
      void save (struct projectData *);
      void save_as_test (string, struct projectData *, int, int *);
      ~PatternDatabase ();
};

#endif
