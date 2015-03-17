#include "clean.h"
#include <igl/barycenter.h>
#include <igl/boundary_facets.h>
#include <igl/tetgen/tetrahedralize.h>
#include <igl/tetgen/cdt.h>
#include <igl/winding_number.h>
#include <igl/unique_simplices.h>
#include <igl/remove_unreferenced.h>
#include <igl/writeOBJ.h>
#include <igl/writeMESH.h>
#include <igl/cgal/remesh_self_intersections.h>

bool clean(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & CV,
  Eigen::MatrixXi & CF,
  Eigen::VectorXi & IM)
{
  using namespace igl;
  using namespace Eigen;
  using namespace std;
  //writeOBJ("VF.obj",V,F);
  const auto & validate_IM = [](
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXd & CV,
    const Eigen::VectorXi & IM)
  {
    for(int i = 0;i<IM.rows();i++)
    {
      if(IM(i)<V.rows())
      {
        double diff = (V.row(IM(i))-CV.row(i)).norm();
        if(diff>1e-6)
        {
          cout<<i<<": "<<IM(i)<<" "<<diff<<endl;
        }
      }
    }
  };
  {
    MatrixXi _1;
    VectorXi _2;
#ifdef VERBOSE
    cerr<<"remesh_self_intersections"<<endl;
#endif
    remesh_self_intersections(V,F,{},CV,CF,_1,_2,IM);
    for_each(CF.data(),CF.data()+CF.size(),[&IM](int & a){a=IM(a);});
    validate_IM(V,CV,IM);
    cout<<"remove_unreferenced"<<endl;
    {
      MatrixXi oldCF = CF;
      unique_simplices(oldCF,CF);
    }
    MatrixXd oldCV = CV;
    MatrixXi oldCF = CF;
    VectorXi newIM;
    remove_unreferenced(oldCV,oldCF,CV,CF,newIM);
    for_each(IM.data(),IM.data()+IM.size(),[&newIM](int & a){a=newIM(a);});
    validate_IM(V,CV,IM);
  }
  MatrixXd TV;
  MatrixXi TT;
  {
    MatrixXi _1;
    // c  convex hull
    // Y  no boundary steiners
    // p  polygon input
#ifdef VERBOSE
    cerr<<"tetrahedralize"<<endl;
    //writeOBJ("CVCF.obj",CV,CF);
#endif
    CDTParam params;
    params.flags = "CY";
    params.use_bounding_box = true;
    if(cdt(CV,CF,params,TV,TT,_1) != 0)
    {
      cout<<REDRUM("CDT failed.")<<endl;
      return false;
    }
    //writeMESH("TVTT.mesh",TV,TT,MatrixXi());
  }
  {
    MatrixXd BC;
    barycenter(TV,TT,BC);
    VectorXd W;
#ifdef VERBOSE
    cerr<<"winding_number"<<endl;
#endif
    winding_number(V,F,BC,W);
    W = W.array().abs();
    const double thresh = 0.5;
    const int count = (W.array()>thresh).cast<int>().sum();
    MatrixXi CT(count,TT.cols());
    int c = 0;
    for(int t = 0;t<TT.rows();t++)
    {
      if(W(t)>thresh)
      {
        CT.row(c++) = TT.row(t);
      }
    }
    assert(c==count);
    boundary_facets(CT,CF);
    //writeMESH("CVCTCF.mesh",TV,CT,CF);
    cout<<"remove_unreferenced"<<endl;
    // Force all original vertices to be referenced
    MatrixXi FF = F;
    for_each(FF.data(),FF.data()+FF.size(),[&IM](int & a){a=IM(a);});
    int ncf = CF.rows();
    MatrixXi ref(ncf+FF.rows(),3);
    ref<<CF,FF;
    VectorXi newIM;
    remove_unreferenced(TV,ref,CV,CF,newIM);
    // Only keep boundary faces
    CF.conservativeResize(ncf,3);
    cout<<"IM.minCoeff(): "<<IM.minCoeff()<<endl;
    for_each(IM.data(),IM.data()+IM.size(),[&newIM](int & a){a=newIM(a);});
    cout<<"IM.minCoeff(): "<<IM.minCoeff()<<endl;
    validate_IM(V,CV,IM);
  }
  return true;
}

