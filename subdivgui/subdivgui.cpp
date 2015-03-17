#include "clean.h"
#include "subdiv_weights.h"
#include "Mesh.h"
extern "C"
{
#include "subdivision_skinning_wrapper.h"
}
#include <igl/Camera.h>
#include <igl/MouseController.h>
#include <igl/REDRUM.h>
#include <igl/ReAntTweakBar.h>
#include <igl/STR.h>
#include <igl/bone_parents.h>
#include <igl/centroid.h>
#include <igl/colon.h>
#include <igl/colon.h>
#include <igl/draw_beach_ball.h>
#include <igl/draw_floor.h>
#include <igl/draw_skeleton_3d.h>
#include <igl/draw_skeleton_vector_graphics.h>
#include <igl/forward_kinematics.h>
#include <igl/get_seconds.h>
#include <igl/lbs_matrix.h>
#include <igl/material_colors.h>
#include <igl/next_filename.h>
#include <igl/normalize_row_sums.h>
#include <igl/pathinfo.h>
#include <igl/polygon_mesh_to_triangle_mesh.h>
#include <igl/quat_to_mat.h>
#include <igl/readDMAT.h>
#include <igl/readTGF.h>
#include <igl/readOBJ.h>
#include <igl/list_to_matrix.h>
#include <igl/remove_unreferenced.h>
#include <igl/report_gl_error.h>
#include <igl/snap_to_canonical_view_quat.h>
#include <igl/snap_to_fixed_up.h>
#include <igl/trackball.h>
#include <igl/two_axis_valuator_fixed_up.h>
#include <igl/writeDMAT.h>
#include <igl/writeOBJ.h>
#include <igl/png/render_to_png_async.h>
#include <igl/writeMESH.h>
#include <igl/writeOFF.h>
#include <igl/writeTGF.h>

#include <Eigen/Core>
#include <Eigen/Geometry>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <iostream>
#include <fstream>

#include <string>
#include <vector>
#include <queue>
#include <stack>
#include <iostream>
#include <iomanip>
#include <unistd.h>

#define VERBOSE
#define DEBUG_WEIGHTS

enum SkelStyleType
{
  SKEL_STYLE_TYPE_3D = 0,
  SKEL_STYLE_TYPE_VECTOR_GRAPHICS = 1,
  NUM_SKEL_STYLE_TYPE = 2
}skel_style;
bool show_skeleton = true;

double fps = 0.;
bool force_anim = false;
bool floor_visible = true;
bool render_to_png_on_next = false;
bool render_to_png_on_anim = false;

Eigen::MatrixXd M;
Mesh c(false),r(true),l(false); // coarse, refined, lbs

Eigen::Vector3d Vmid;
double bbd = 1.0;
igl::Camera camera;

Eigen::MatrixXd C;
Eigen::MatrixXi BE;
Eigen::VectorXi P,RP;

struct State
{
  igl::MouseController mouse;
  Eigen::MatrixXf colors;
} s;

bool wireframe = false;

// See README for descriptions
enum RotationType
{
  ROTATION_TYPE_IGL_TRACKBALL = 0,
  ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP = 1,
  NUM_ROTATION_TYPES = 2,
} rotation_type = ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP;

std::stack<State> undo_stack;
std::stack<State> redo_stack;

bool is_rotating = false;
bool centroid_is_visible = false;
int down_x,down_y;
igl::Camera down_camera;
std::string output_prefix;

struct CameraAnimation
{
  bool is_animating = false;
  double DURATION = 0.5;
  double start_time = 0;
  Eigen::Quaterniond from_quat,to_quat;
} canim;

typedef std::vector<
Eigen::Quaterniond,
  Eigen::aligned_allocator<Eigen::Quaterniond> > RotationList;

struct PoseAnimation
{
  bool is_animating = false;
  double DURATION = 2;
  double start_time = 0;
  RotationList pose;
  int count = 0;
} panim;

int width,height;
Eigen::Vector4f light_pos(-0.1,-0.1,0.9,0);

#define REBAR_NAME "temp.rbr"
igl::ReTwBar rebar;

// Subdiv library pointers
void* engine;
void* eval;

bool show_cage_on_drag = true;

void push_undo()
{
  undo_stack.push(s);
  // Clear
  redo_stack = std::stack<State>();
}

// No-op setter, does nothing
void TW_CALL no_op(const void * /*value*/, void * /*clientData*/)
{
}

void TW_CALL set_rotation_type(const void * value, void * clientData)
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;
  const RotationType old_rotation_type = rotation_type;
  rotation_type = *(const RotationType *)(value);
  if(rotation_type == ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP &&
    old_rotation_type != ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP)
  {
    push_undo();
    canim.from_quat = camera.m_rotation_conj;
    snap_to_fixed_up(canim.from_quat,canim.to_quat);
    // start animation
    canim.start_time = get_seconds();
    canim.is_animating = true;
  }
}
void TW_CALL get_rotation_type(void * value, void *clientData)
{
  RotationType * rt = (RotationType *)(value);
  *rt = rotation_type;
}

void reshape(int width, int height)
{
  ::width = width;
  ::height = height;
  glViewport(0,0,width,height);
  // Send the new window size to AntTweakBar
  TwWindowSize(width, height);
  camera.m_aspect = (double)width/(double)height;
  s.mouse.reshape(width,height);
}
    

void push_scene()
{
  using namespace igl;
  using namespace std;
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  gluPerspective(camera.m_angle,camera.m_aspect,camera.m_near,camera.m_far);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  gluLookAt(
    camera.eye()(0), camera.eye()(1), camera.eye()(2),
    camera.at()(0), camera.at()(1), camera.at()(2),
    camera.up()(0), camera.up()(1), camera.up()(2));
}

void push_object()
{
  using namespace igl;
  glPushMatrix();
  glScaled(2./bbd,2./bbd,2./bbd);
  glTranslated(-Vmid(0),-Vmid(1),-Vmid(2));
}

void pop_object()
{
  glPopMatrix();
}

void pop_scene()
{
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}

// Set up double-sided lights
void lights()
{
  using namespace std;
  using namespace Eigen;
  glEnable(GL_LIGHTING);
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);
  float WHITE[4] =  {0.8,0.8,0.8,1.};
  float GREY[4] =  {0.4,0.4,0.4,1.};
  float BLACK[4] =  {0.,0.,0.,1.};
  Vector4f pos = light_pos;
  glLightfv(GL_LIGHT0,GL_AMBIENT,GREY);
  glLightfv(GL_LIGHT0,GL_DIFFUSE,WHITE);
  glLightfv(GL_LIGHT0,GL_SPECULAR,BLACK);
  glLightfv(GL_LIGHT0,GL_POSITION,pos.data());
  pos(0) *= -1;
  pos(1) *= -1;
  pos(2) *= -1;
  glLightfv(GL_LIGHT1,GL_AMBIENT,GREY);
  glLightfv(GL_LIGHT1,GL_DIFFUSE,WHITE);
  glLightfv(GL_LIGHT1,GL_SPECULAR,BLACK);
  glLightfv(GL_LIGHT1,GL_POSITION,pos.data());
}

template <typename DerivedQ, typename DerivedU, typename DerivedUN>
void update_refined(
  const Eigen::Matrix<subdivision_evaluator_real_t,Eigen::Dynamic,Eigen::Dynamic> & TT,
  const Eigen::PlainObjectBase<DerivedQ> & Q,
  subdivision_evaluator_real_t * vs_modified,
  Eigen::PlainObjectBase<DerivedU> & U,
  Eigen::PlainObjectBase<DerivedUN> & UN)
{
  using namespace Eigen;
  using namespace igl;

  // Retrieve refined mesh
  const int num_refined_vs = num_refined_vertices_of_subdivision_evaluator( eval );
  subdivision_evaluator_real_t* refined_vs = new subdivision_evaluator_real_t[ num_refined_vs*3 ];
  get_refined_vertices_of_subdivision_skinning_engine_with_control_vertices( engine, 3, vs_modified, refined_vs );
  const auto & UT = Map<Matrix<subdivision_evaluator_real_t,Dynamic,Dynamic> >
    (refined_vs,3,num_refined_vs);
  U = UT.cast<double>().transpose();
  // Clean up
  delete[] refined_vs;
}

void display()
{
  using namespace igl;
  using namespace std;
  using namespace Eigen;
  //const float back[4] = {0.75, 0.75, 0.75,0};
  const float back[4] = {1.,1.,1.,0};
  glClearColor(back[0],back[1],back[2],0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  if(canim.is_animating)
  {
    double t = (get_seconds() - canim.start_time)/canim.DURATION;
    if(t > 1)
    {
      t = 1;
      canim.is_animating = false;
    }
    Quaterniond q = canim.from_quat.slerp(t,canim.to_quat).normalized();
    camera.orbit(q.conjugate());
  }

  RotationList dQ;
  if(panim.is_animating)
  {
    double t = 
      render_to_png_on_anim ?
      double(panim.count)/(30.*panim.DURATION) :
      (get_seconds() - panim.start_time)/panim.DURATION;
    if(t > 1)
    {
      t = 1;
      panim.is_animating = false;
    }
    const auto & ease = [](const double t)
    {
      return 3.*t*t-2.*t*t*t;
    };
    double f = (t<0.5?ease(2.*t):ease(2.-2.*t));
    dQ.resize(panim.pose.size());
    for(int e = 0;e<(int)panim.pose.size();e++)
    {
      dQ[e] = panim.pose[e].slerp(f,Quaterniond::Identity()).normalized();
    }
    panim.count++;
  }else
  {
    dQ = s.mouse.rotations();
  }
  RotationList vQ;
  vector<Vector3d> vT;
  forward_kinematics(C,BE,P,dQ,vQ,vT);
    
  // Vertical stack, row-major --> horizontal stack of transposes, col-major
  Matrix<subdivision_evaluator_real_t,Dynamic,Dynamic> TT(4,4*BE.rows());
  // draw_skeleton needs a different stack
  MatrixXd T(BE.rows()*(3+1),3);
  for(int w = 0;w<BE.rows();w++)
  {
    Affine3d a = Affine3d::Identity();
    a.translate(vT[w]);
    a.rotate(vQ[w]);
    TT.block(0,w*4,4,4) = 
      a.matrix().transpose().cast<subdivision_evaluator_real_t>();
    T.block(w*(3+1),0,3+1,3) = a.matrix().transpose().block(0,0,3+1,3);
  }
  /////////////////////////////////////////////////////////////////////////
  // Tell subdiv lib about transformations and get modified cage, refined
  /////////////////////////////////////////////////////////////////////////
  if(c.stale || r.stale)
  {
    const subdivision_evaluator_real_t* transforms = TT.data();
    subdivision_evaluator_real_t* vs_modified = new subdivision_evaluator_real_t[ c.V.size() ];
    // Retrieve control mesh
    compute_control_mesh_vertices_given_transforms_for_subdivision_skinning_engine( engine, (subdivision_evaluator_real_t*)transforms, vs_modified );
    const auto & CUT = Map<Matrix<subdivision_evaluator_real_t,Dynamic,Dynamic> >
      (vs_modified,3,c.V.rows());
    c.U = CUT.cast<double>().transpose();
    if((!s.mouse.is_widget_down() && r.must_update && show_cage_on_drag) || 
      r.force_visible)
    {
      if(r.stale)
      {
        update_refined(TT,r.Q,vs_modified,r.U,r.N);
      }
      r.must_update = false;
    }
    delete[] vs_modified;
  }

  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);
  glEnable(GL_NORMALIZE);
  lights();
  push_scene();

  double y_boost = 0;
  // Draw a nice floor
  if(floor_visible)
  {
    if(c.U.rows() > 0)
    { 
      y_boost = 2.*(c.U.col(1).minCoeff()-c.V.col(1).minCoeff())/bbd;
    }
    if(panim.is_animating && render_to_png_on_anim)
    {
      glTranslated(0,-y_boost,0);
    }
    glEnable(GL_DEPTH_TEST);
    glPushMatrix();
    double floor_offset = -2./bbd*(c.V.col(1).maxCoeff()-Vmid(1));
    floor_offset += y_boost;
    glTranslated(0,floor_offset,0);
    const float GREY[4] = {0.5,0.5,0.6,1.0};
    const float DARK_GREY[4] = {0.2,0.2,0.3,1.0};
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    static GLuint floor_dl = 0;
    if(!glIsList(floor_dl))
    {
      floor_dl = glGenLists(1);
      glNewList(floor_dl,GL_COMPILE_AND_EXECUTE);
      draw_floor(GREY,DARK_GREY);
      glEndList();
    }else
    {
      glCallList(floor_dl);
    }
    glDisable(GL_CULL_FACE);
    glPopMatrix();
  }

  push_object();

  const auto & draw_skeleton = [](const MatrixXd & T)
  {
    switch(skel_style)
    {
      default:
      case SKEL_STYLE_TYPE_3D:
      {
        draw_skeleton_3d(C,BE,T,s.colors);
        break;
      }
      case SKEL_STYLE_TYPE_VECTOR_GRAPHICS:
        draw_skeleton_vector_graphics(C,BE,T);
        break;
    }
  };
  // Set material properties
  glDisable(GL_COLOR_MATERIAL);
  glMaterialfv(GL_FRONT, GL_AMBIENT,GOLD_AMBIENT);
  glMaterialfv(GL_FRONT, GL_DIFFUSE,GOLD_DIFFUSE);
  glMaterialfv(GL_FRONT, GL_SPECULAR,GOLD_SPECULAR);
  glMaterialf (GL_FRONT, GL_SHININESS, 128);
  glMaterialfv(GL_BACK, GL_AMBIENT,SILVER_AMBIENT);
  glMaterialfv(GL_BACK, GL_DIFFUSE,FAST_GREEN_DIFFUSE);
  glMaterialfv(GL_BACK, GL_SPECULAR,SILVER_SPECULAR);
  glMaterialf (GL_BACK, GL_SHININESS, 128);
  if(wireframe)
  {
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
  }
  glLineWidth(1.0);

  if((r.must_update && show_cage_on_drag) || c.force_visible)
  {
    if(c.stale)
    {
      c.N.conservativeResize(c.Q.rows(),c.N.cols());
    }
    c.draw_and_cache();
  }

  if((!r.must_update && show_cage_on_drag) || r.force_visible)
  {
    r.draw_and_cache();
  }

  if(l.force_visible)
  {
    // Set material properties
    glDisable(GL_COLOR_MATERIAL);
    glMaterialfv(GL_FRONT, GL_AMBIENT,SILVER_AMBIENT);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,FAST_RED_DIFFUSE);
    glMaterialfv(GL_FRONT, GL_SPECULAR,SILVER_SPECULAR);
    glMaterialf (GL_FRONT, GL_SHININESS, 128);
    glMaterialfv(GL_BACK, GL_AMBIENT,SILVER_AMBIENT);
    glMaterialfv(GL_BACK, GL_DIFFUSE,SILVER_DIFFUSE);
    glMaterialfv(GL_BACK, GL_SPECULAR,SILVER_SPECULAR);
    glMaterialf (GL_BACK, GL_SHININESS, 128);
    if(l.stale)
    {
      l.U = M*T;
    }
    l.draw_and_cache();
  }

  glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
  glDisable(GL_DEPTH_TEST);
  if(show_skeleton)
  {
    draw_skeleton(T);
  }

  if(centroid_is_visible && c.U.rows() > 0)
  {
    Vector3d cen;
    centroid(c.U,c.F,cen);
    glEnable(GL_DEPTH_TEST);
    glPushMatrix();
    glTranslated(cen(0),cen(1),cen(2));
    glScaled(bbd/2.,bbd/2.,bbd/2.);
    glScaled(0.1,0.1,0.1);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(0,-100000);
    draw_beach_ball();
    glDisable(GL_POLYGON_OFFSET_FILL);
    glPopMatrix();
  }

  // Mouse is always on top
  glDisable(GL_DEPTH_TEST);
  if(!panim.is_animating)
  {
    s.mouse.draw();
  }
  pop_object();
  pop_scene();

  report_gl_error();

  bool was_rendering_to_png = false;
  if(render_to_png_on_next || (render_to_png_on_anim && panim.is_animating)) 
  {
    was_rendering_to_png = true;
    string prefix = STR(getenv("HOME")<<"/Desktop/"<< "subdivgui-");
    string filename;
    next_filename(prefix,4,".png",filename);
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT,viewport);
    render_to_png_async(filename.c_str(),viewport[2],viewport[3],true,false);
    if(render_to_png_on_next)
    {
      glClearColor(1,1,1,1);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      glFlush();
      usleep(1000*1);
      render_to_png_on_next = false;
    }
  }

  TwDraw();
  glutSwapBuffers();

  if(canim.is_animating || 
    panim.is_animating || 
    force_anim || 
    was_rendering_to_png)
  {
    r.must_update = true;
    c.stale = l.stale = r.stale = true;
    glutPostRedisplay();
  }
  {
    static int count = 0;
    static double start = get_seconds();
    if(++count > 10)
    {
      double now = get_seconds();
      fps = double(count)/(now-start);
      count = 0;
      start = now;
    }
  }
}

void mouse_wheel(int wheel, int direction, int mouse_x, int mouse_y)
{
  using namespace std;
  using namespace igl;
  using namespace Eigen;
  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT,viewport);
  if(wheel == 0 && TwMouseMotion(mouse_x, viewport[3] - mouse_y))
  {
    static double mouse_scroll_y = 0;
    const double delta_y = 0.125*direction;
    mouse_scroll_y += delta_y;
    TwMouseWheel(mouse_scroll_y);
    return;
  }
  push_undo();

  if(wheel==0)
  {
    // factor of zoom change
    double s = (1.-0.01*direction);
    //// FOV zoom: just widen angle. This is hardly ever appropriate.
    //camera.m_angle *= s;
    //camera.m_angle = min(max(camera.m_angle,1),89);
    camera.push_away(s);
  }else
  {
    // Dolly zoom:
    camera.dolly_zoom((double)direction*1.0);
  }
  glutPostRedisplay();
}

void mouse(int glutButton, int glutState, int mouse_x, int mouse_y)
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;
  bool tw_using = TwEventMouseButtonGLUT(glutButton,glutState,mouse_x,mouse_y);
  const int mod = (glutButton <=2 ? glutGetModifiers() : 0);
  const bool option_down = mod & GLUT_ACTIVE_ALT;
  switch(glutButton)
  {
    case GLUT_RIGHT_BUTTON:
    case GLUT_LEFT_BUTTON:
    {
      push_scene();
      push_object();
      switch(glutState)
      {
        case 1:
        {
          // up
          const bool mouse_was_selecting = s.mouse.is_selecting();
          is_rotating = false;
          s.mouse.up(mouse_x,mouse_y);
          glutSetCursor(GLUT_CURSOR_INHERIT);
          if(mouse_was_selecting)
          {
            s.mouse.set_selection_from_last_drag(C,BE,P,RP);
            MouseController::VectorXb S;
            MouseController::propogate_to_descendants_if(
              s.mouse.selection(),P,S);
            vector<size_t> selected_weights;
            for(size_t s = 0;s<S.rows();s++)
            {
              if(S(s))
              {
                selected_weights.push_back(s);
              }
            }
            if(selected_weights.size() > 0)
            {
              r.selected_weights = selected_weights;
              if(r.show_weights)
              {
                r.stale = true;
              }
            }
            MouseController::color_if(S,MAYA_SEA_GREEN,MAYA_VIOLET,s.colors);
          }
          break;
        }
        case 0:
          if(!tw_using)
          {
            down_x = mouse_x;
            down_y = mouse_y;
            if(option_down || glutButton==GLUT_RIGHT_BUTTON)
            {
              glutSetCursor(GLUT_CURSOR_CYCLE);
              // collect information for trackball
              is_rotating = true;
              down_camera = camera;
            }else
            {
              push_undo();
              s.mouse.down(mouse_x,mouse_y);
              if(s.mouse.is_widget_down())
              {
                r.must_update = true;
                c.stale = l.stale = r.stale = true;
              }
            }
          }
        break;
      }
      pop_object();
      pop_scene();
      break;
    }
    // Scroll down
    case 3:
    {
      mouse_wheel(0,-1,mouse_x,mouse_y);
      break;
    }
    // Scroll up
    case 4:
    {
      mouse_wheel(0,1,mouse_x,mouse_y);
      break;
    }
    // Scroll left
    case 5:
    {
      mouse_wheel(1,-1,mouse_x,mouse_y);
      break;
    }
    // Scroll right
    case 6:
    {
      mouse_wheel(1,1,mouse_x,mouse_y);
      break;
    }
  }
  glutPostRedisplay();
}

void mouse_drag(int mouse_x, int mouse_y)
{
  using namespace igl;
  using namespace std;
  using namespace Eigen;

  push_scene();
  push_object();
  if(is_rotating)
  {
    glutSetCursor(GLUT_CURSOR_CYCLE);
    Quaterniond q;
    switch(rotation_type)
    {
      case ROTATION_TYPE_IGL_TRACKBALL:
      {
        // Rotate according to trackball
        igl::trackball(
          width,
          height,
          2.0,
          down_camera.m_rotation_conj,
          down_x,
          down_y,
          mouse_x,
          mouse_y,
          q);
          break;
      }
      case ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP:
      {
        // Rotate according to two axis valuator with fixed up vector
        two_axis_valuator_fixed_up(
          width, height,
          2.0,
          down_camera.m_rotation_conj,
          down_x, down_y, mouse_x, mouse_y,
          q);
        break;
      }
      default:
        break;
    }
    camera.orbit(q.conjugate());
  }else if(s.mouse.drag(mouse_x,mouse_y))
  {
    if(s.mouse.is_widget_down())
    {
      c.stale = l.stale = r.stale = true;
    }
  }
  pop_object();
  pop_scene();
  glutPostRedisplay();
}

void init_relative()
{
  using namespace Eigen;
  using namespace igl;
  using namespace std;
  const auto Vmax = c.V.colwise().maxCoeff();
  const auto Vmin = c.V.colwise().minCoeff();
  Vmid = 0.5*(Vmax + Vmin);
  bbd = (Vmax-Vmin).norm();
  camera.push_away(2);
}

void undo()
{
  using namespace std;
  if(!undo_stack.empty())
  {
    redo_stack.push(s);
    s = undo_stack.top();
    undo_stack.pop();
    s.mouse.reshape(width,height);
  }
}

void redo()
{
  using namespace std;
  if(!redo_stack.empty())
  {
    undo_stack.push(s);
    s = redo_stack.top();
    redo_stack.pop();
    s.mouse.reshape(width,height);
  }
}

bool save()
{
  using namespace std;
  using namespace igl;
  using namespace Eigen;
  string output_filename;
  next_filename(output_prefix,4,".dmat",output_filename);
  MatrixXd T;
  forward_kinematics(C,BE,P,s.mouse.rotations(),T);
  if(writeDMAT(output_filename,T))
  {
    cout<<GREENGIN("Current pose written to "+output_filename+".")<<endl;
    return true;
  }else
  {
    cout<<REDRUM("Writing to "+output_filename+" failed.")<<endl;
    return false;
  }
}

void key(unsigned char key, int mouse_x, int mouse_y)
{
  using namespace std;
  using namespace igl;
  using namespace Eigen;
  r.must_update = true;
  c.stale = l.stale = r.stale = true;
  int mod = glutGetModifiers();
  const bool command_down = GLUT_ACTIVE_COMMAND & mod;
  const bool shift_down = GLUT_ACTIVE_SHIFT & mod;
  switch(key)
  {
    // ESC
    case char(27):
      rebar.save(REBAR_NAME);
    // ^C
    case char(3):
      exit(EXIT_SUCCESS);
    case 'a':
    {
      panim.is_animating = !panim.is_animating;
      panim.pose = s.mouse.rotations();
      panim.start_time = get_seconds();
      panim.count = 0;
      break;
    }
    case 'D':
    case 'd':
    {
      push_undo();
      s.mouse.clear_selection();
      break;
    }
    case 'P':
    case 'p':
    {
      render_to_png_on_next = true;
      break;
    }
    case 'R':
    {
      push_undo();
      s.mouse.reset_selected_rotations();
      break;
    }
    case 'r':
    {
      push_undo();
      s.mouse.reset_rotations();
      break;
    }
    //case 'S':
    //case 's':
    //{
    //  save();
    //  break;
    //}
    case 'z':
    case 'Z':
      is_rotating = false;
      if(command_down)
      {
        if(shift_down)
        {
          redo();
        }else
        {
          undo();
        }
        break;
      }else
      {
        push_undo();
        Quaterniond q;
        snap_to_canonical_view_quat(camera.m_rotation_conj,1.0,q);
        camera.orbit(q.conjugate());
      }
      break;
    default:
      if(!TwEventKeyboardGLUT(key,mouse_x,mouse_y))
      {
        cout<<"Unknown key command: "<<key<<" "<<int(key)<<endl;
      }
  }

  glutPostRedisplay();
}



int main(int argc, char * argv[])
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;
  string filename = "torus.obj";
  string skel_filename = "torus.tgf";
  string weights_filename = "";
  output_prefix = "";
  int comp_level = 1;
  int eval_level = 3;
  switch(argc)
  {
    case 7:
      output_prefix = argv[6];
      //fall through
    case 6:
      weights_filename = argv[5];
    case 5:
      eval_level = atoi(argv[4]);
    case 4:
      comp_level = atoi(argv[3]);
      //fall through
    case 3:
      skel_filename = argv[2];
      // Read and prepare mesh
      filename = argv[1];
      break;
    default:
      cerr<<"Usage:"<<endl<<"    ./example model.obj skeleton.tgf "<<
        "[comp_level] [eval_level] [weights.dmat] [pose-prefix]"<<endl;
      cout<<endl<<"Opening default rig..."<<endl;
      weights_filename = "torus-weights.dmat";
  }
  if(comp_level>eval_level)
  {
    cerr<<"comp_level ("<<comp_level<<") must be <= eval_level ("<<eval_level
      <<")"<<endl;
    return EXIT_FAILURE;
  }

  // print key commands
  cout<<"[Click] and [drag]     Select a bone/Use onscreen widget to rotate bone."<<endl;
  cout<<"⌥ +[Click] and [drag]  Rotate secene."<<endl;
  cout<<"[,]                    Change scene rotation mouse UI method."<<endl;
  cout<<"A                      Force tight animation rendering loop."<<endl;
  cout<<"a                      Animate from current pose to rest pose."<<endl;
  cout<<"C,c                    Visualize centroid."<<endl;
  cout<<"D,d                    Deselect all."<<endl;
  cout<<"F,f                    Show floor."<<endl;
  cout<<"G                      Force show coarse cage."<<endl;
  cout<<"g                      Force show refined mesh."<<endl;
  cout<<"H,h                    Show coarse cage while dragging."<<endl;
  cout<<"l                      Show wireframe on all."<<endl;
  cout<<"L                      Show wireframe on cage."<<endl;
  cout<<"P,p                    Render next frame to png file."<<endl;
  cout<<"R                      Reset selected rotation(s)."<<endl;
  cout<<"r                      Reset all rotations."<<endl;
  cout<<"S,s                    Save current pose."<<endl;
  cout<<"V,v                    Render directly skinned \"limit\" surface."<<endl;
  cout<<"Z,z                    Snap to canonical view."<<endl;
  cout<<"⌘ Z                    Undo."<<endl;
  cout<<"⇧ ⌘ Z                  Redo."<<endl;
  cout<<"^C,ESC                 Exit (without saving)."<<endl;

  // Load subdivision cage
  string dir,_1,_2,name;
  pathinfo(filename,dir,_1,_2,name);
  {
    vector<vector<int > > vCQ,vCF;
    vector<vector<double> > vCV;
    if(!readOBJ(filename,vCV,vCQ))
    {
      return EXIT_FAILURE;
    }
    list_to_matrix(vCV,c.V);
    polygon_mesh_to_triangle_mesh(vCQ,c.F);
    if(!list_to_matrix(vCQ,4,-1,c.Q))
    {
      cerr<<"Error: "<<filename<<" contains high valence facet."<<endl;
      return EXIT_FAILURE;
    }
  }

  if(output_prefix.size() == 0)
  {
    output_prefix = dir+"/"+name+"-pose-";
  }

  {
    string output_filename;
    next_filename(output_prefix,4,".dmat",output_filename);
    cout<<BLUEGIN("Output set to start with "<<output_filename)<<endl;
  }

  /////////////////////////////////////////////////////////////////////////
  // Tell subdiv lib about cage
  /////////////////////////////////////////////////////////////////////////
  const Matrix<subdivision_evaluator_real_t,Dynamic,Dynamic> CVT = 
    c.V.cast<subdivision_evaluator_real_t>().transpose();
  const subdivision_evaluator_real_t * vs = CVT.data();
  const int num_vs = c.V.rows();
  const MatrixXi CQT = c.Q.cast<int>().transpose();
  const int * faces = CQT.data();
  const int num_faces = c.Q.rows();
  eval = new_subdivision_evaluator( 
   num_vs, (subdivision_evaluator_real_t*)vs, num_faces, (int*)faces, eval_level);
    
  /////////////////////////////////////////////////////////////////////////
  // Ask subdiv lib for refined mesh
  /////////////////////////////////////////////////////////////////////////
  const int num_refined_faces = num_refined_quad_faces_of_subdivision_evaluator( eval );
  int* refined_faces = new int[ num_refined_faces*4 ];
  get_refined_quad_faces_of_subdivision_evaluator( eval, refined_faces );
  const auto & QT = Map<MatrixXi>(refined_faces,4,num_refined_faces);
  r.Q = QT.transpose();
  polygon_mesh_to_triangle_mesh(r.Q,r.F);
  delete[] refined_faces;
    
  const int num_refined_vs = num_refined_vertices_of_subdivision_evaluator( eval );
  subdivision_evaluator_real_t* refined_vs = new subdivision_evaluator_real_t[ num_refined_vs*3 ];
  get_refined_vertices_of_subdivision_evaluator( eval, refined_vs );
  const auto & VT = Map<Matrix<subdivision_evaluator_real_t,Dynamic,Dynamic> >
    (refined_vs,3,num_refined_vs);
  r.V = VT.cast<double>().transpose();
  delete[] refined_vs;

  /////////////////////////////////////////////////////////////////////////
  // Load skeleton and load/compute weights
  /////////////////////////////////////////////////////////////////////////
  // Read in skeleton and precompute hierarchy
  readTGF(skel_filename,C,BE);
  // initialize mouse interface
  s.mouse.set_size(BE.rows());
  // Rigid parts (not used)
  colon<int>(0,BE.rows()-1,RP);
  assert(RP.size() == BE.rows());
  // Bone parents
  bone_parents(BE,P);
  if(weights_filename.size() == 0 || !file_exists(weights_filename.c_str()))
  {
    cout<<YELLOWGIN("Computing weights...")<<endl;
    //robust_bbw(V,F,C,BE,r.W);
    subdiv_weights(c.V,c.Q,C,BE,comp_level,eval_level,WEIGHTS_TYPE_BBW,r.W);
    if(weights_filename.size() > 0)
    {
      if(writeDMAT(weights_filename,r.W))
      {
        cout<<GREENGIN("Saved weights to "<<weights_filename)<<endl;
      }else
      {
        cerr<<REDRUM("Failed to save weights to "<<weights_filename)<<endl;
      }
    }
  }else
  {
    // Read in weights and precompute LBS matrix
    readDMAT(weights_filename,r.W);
  }

  assert(r.W.rows() == r.V.rows() && "#W must match #r.V");
  // Use refined mesh for lbs
  l = r;
  lbs_matrix(l.V,r.W,M);
  if(r.W.cols() != BE.rows())
  {
    cerr<<REDRUM("# cols in weights ("<<r.W.cols()<<
      ") != # bones ("<<BE.rows()<<")")<<endl;
    return EXIT_FAILURE;
  }

  if(num_refined_vs != r.W.rows())
  {
    cerr<<REDRUM("# rows in weights ("<<r.W.rows()<<
      ") != # refined vertices ("<<
      num_refined_vs<<") at level ("<<
      eval_level<<")")<<endl;
    return EXIT_FAILURE;
  }
    
  /////////////////////////////////////////////////////////////////////////
  // Tell subdiv lib about weights
  /////////////////////////////////////////////////////////////////////////
  const Matrix<subdivision_evaluator_real_t,Dynamic,Dynamic> WT = 
    r.W.cast<subdivision_evaluator_real_t>().transpose();
  const subdivision_evaluator_real_t * weights = WT.data();
  engine = new_subdivision_skinning_engine( eval, r.W.cols(), weights );


  // Handle degenerate quads before draw (these only occur in coarse)
  for(int i = 0;i<c.Q.rows();i++)
  {
    if(c.Q(i,3) == -1)
    {
      c.Q(i,3) = c.Q(i,2);
    }
  }

  // initialize some camera settings
  init_relative();


  // Init glut
  glutInit(&argc,argv);
  if( !TwInit(TW_OPENGL, NULL) )
  {
    // A fatal error occured
    fprintf(stderr, "AntTweakBar initialization failed: %s\n", TwGetLastError());
    return EXIT_FAILURE;
  }
  // Create a tweak bar
  rebar.TwNewBar("TweakBar");
  rebar.TwAddVarRO("fps", TW_TYPE_DOUBLE,&fps,"");
  rebar.TwAddVarRW("force_anim", TW_TYPE_BOOLCPP,&force_anim,"key=A");
  rebar.TwAddVarRW("camera_rotation", TW_TYPE_QUAT4D,
    camera.m_rotation_conj.coeffs().data(), "open readonly=true");
  TwType RotationTypeTW = ReTwDefineEnumFromString("RotationType",
    "igl_trackball,two-a...-fixed-up");
  rebar.TwAddVarCB( "rotation_type", RotationTypeTW,
    set_rotation_type,get_rotation_type,NULL,"keyIncr=] keyDecr=[");
  rebar.TwAddVarRW("wireframe", TW_TYPE_BOOLCPP,&wireframe,"key=l");
  rebar.TwAddVarRW("show_cage_on_drag", TW_TYPE_BOOLCPP,&show_cage_on_drag,
    "keyIncr=H keyDecr=h");
  rebar.TwAddVarRW("c_force_visible", TW_TYPE_BOOLCPP,&c.force_visible,"key=G");
  rebar.TwAddVarRW("c_wire_frame", TW_TYPE_BOOLCPP,&c.wireframe,"key=L");
  rebar.TwAddVarRW("r_force_visible", TW_TYPE_BOOLCPP,&r.force_visible,"key=g");
  rebar.TwAddVarRW("l_force_visible", TW_TYPE_BOOLCPP,&l.force_visible,
    "keyIncr=V keyDecr=v");
  rebar.TwAddVarRW("floor_visible", TW_TYPE_BOOLCPP,&floor_visible,"keyIncr=f keyDecr=F");
  rebar.TwAddVarRW("render_to_png_on_anim", TW_TYPE_BOOLCPP,&render_to_png_on_anim,"key='CTRL+a'",false);
  rebar.TwAddVarRW("centroid_is_visible", TW_TYPE_BOOLCPP,&centroid_is_visible,
    "keyIncr=C keyDecr=c label='centroid visible?'");
  TwType SkelStyleTypeTW = ReTwDefineEnumFromString("SkelStyleType",
    "3d,vector-graphics");
  rebar.TwAddVarRW("style",SkelStyleTypeTW,&skel_style,"");
  rebar.TwAddVarRW("r_show_weights", TW_TYPE_BOOLCPP,&r.show_weights,"key=w");
  rebar.TwAddVarRW("show_skeleton", TW_TYPE_BOOLCPP,&show_skeleton,"keyIncr=s keyDecr=S");
  rebar.load(REBAR_NAME);

  // Init antweakbar
  glutInitDisplayString( "rgba depth double samples>=8 ");
  glutInitWindowSize(glutGet(GLUT_SCREEN_WIDTH)/2.0,glutGet(GLUT_SCREEN_HEIGHT)/2.0);
  glutCreateWindow("skeleton-poser");
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(key);
  glutMouseFunc(mouse);
  glutMotionFunc(mouse_drag);
  glutPassiveMotionFunc((GLUTmousemotionfun)TwEventMouseMotionGLUT);
  glutMainLoop();

  // Clean up
  delete_subdivision_skinning_engine( engine );
  delete_subdivision_evaluator( eval );
  return EXIT_SUCCESS;
}
