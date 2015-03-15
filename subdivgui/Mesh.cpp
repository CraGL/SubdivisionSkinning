#include "Mesh.h"
#include <igl/per_vertex_normals.h>
#include <igl/polygon_mesh_to_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/draw_mesh.h>
#include <Eigen/Geometry>

void Mesh::draw_and_cache()
{
  using namespace Eigen;
  using namespace igl;
  if(wireframe || show_weights)
  {
    glPushAttrib(GL_POLYGON_BIT);
    glPushAttrib(GL_ENABLE_BIT);
    glPushAttrib(GL_LIGHTING_BIT);
    glPushAttrib(GL_LINE_BIT);
    if(wireframe)
    {
      glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    }
    glEnable(GL_COLOR_MATERIAL);
    if(show_weights)
    {
      glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT);
      glColor3f(0.1,0.1,0.1);
      glColorMaterial(GL_FRONT_AND_BACK,GL_SPECULAR);
      glColor3f(0.4,0.4,0.4);
    }else
    {
      glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT);
      glColor3f(0.,0,0);
      glColorMaterial(GL_FRONT_AND_BACK,GL_SPECULAR);
      glColor3f(0.,0,0);
    }
    glColorMaterial(GL_FRONT_AND_BACK,GL_DIFFUSE);
    if(wireframe)
    {
      glColor3f(0,0,0);
    }
    if(!glIsList(dl) || stale)
    {
      assert(V.rows() == U.rows());
      if(show_weights && W.cols()>0)
      {
        per_vertex_normals(U,F,N);
        C.resize(V.rows(),3);
        for(size_t v = 0;v<V.rows();v++)
        {
          C(v,0) = 1;
          double w = 0;
          for(size_t selected_weight : selected_weights)
          {
            w+=W(v,selected_weight);
          }
          C(v,1) = 1.-w;
          C(v,2) = 1.-w;
        }
      }
      if(glIsList(dl))
      {
        glDeleteLists(dl,1);
      }
      dl = glGenLists(1);
      glNewList(dl,GL_COMPILE_AND_EXECUTE);
      if(show_weights)
      {
        draw_mesh(U,F,N,C);
      }else
      { 
        draw_mesh(U,Q,N);
      }
      glEndList();
    }else
    {
      glCallList(dl);
    }
    glPopAttrib();
    glPopAttrib();
    glPopAttrib();
    glPopAttrib();
  } else
  {

    if(!glIsBuffer(ibo))
    {
      T.resize(Q.rows()*2,3);
      Matrix<GLuint,Dynamic,Dynamic,RowMajor> TCT(T.rows(),3);
      for(int q=0;q<Q.rows();q++)
      {
        T(q*2+0,0) = Q(q,0);
        T(q*2+0,1) = Q(q,1);
        T(q*2+0,2) = Q(q,2);
        T(q*2+1,0) = Q(q,0);
        T(q*2+1,1) = Q(q,2);
        T(q*2+1,2) = Q(q,3)==-1?Q(q,2):Q(q,3);
      }
      T1 = Q.leftCols(3).eval();
      for(int f=0;f<T.rows();f++)
      {
        for(int c = 0;c<3;c++)
        {
          TCT(f,c) = f*3+c;
        }
      }
      UCT.resize(3*T.rows(),3);
      NCT.resize(3*T.rows(),3);
      glGenBuffers(1,&ibo);
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,ibo);
      glBufferData(GL_ELEMENT_ARRAY_BUFFER,sizeof(GLuint)*TCT.size(),TCT.data(),GL_STATIC_DRAW);
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
      glGenBuffers(1,&vbo);
      glGenBuffers(1,&nbo);
    }

    const int Trows = T.rows();
#   pragma omp parallel for if (Trows>10000)
    for(int f=0;f<Trows;f+=2)
    {
      for(int c = 0;c<3;c++)
      {
        UCT.row(f*3+c) = U.row(T(f,c)).cast<float>();
        UCT.row((f+1)*3+c) = U.row(T(f+1,c)).cast<float>();
      }
      Matrix<float,1,3,RowMajor> v1 = UCT.row(f*3+1)-UCT.row(f*3+0);
      Matrix<float,1,3,RowMajor> v2 = UCT.row(f*3+2)-UCT.row(f*3+0);
      Matrix<float,1,3,RowMajor> n = v1.cross(v2);
      for(int c = 0;c<3;c++)
      {
        for(int d = 0;d<3;d++)
        {
          const float nd = n(d);
          NCT(f*3+c,d) = nd;
          NCT((f+1)*3+c,d) = nd;
        }
      }
    }

      glEnableClientState(GL_VERTEX_ARRAY);
      glBindBuffer(GL_ARRAY_BUFFER,vbo);
      glBufferData(GL_ARRAY_BUFFER,sizeof(float)*UCT.size(),UCT.data(),GL_DYNAMIC_DRAW);
      glVertexPointer(3,GL_FLOAT,0,0);
      glEnableClientState(GL_NORMAL_ARRAY);
      glBindBuffer(GL_ARRAY_BUFFER,nbo);
      glBufferData(GL_ARRAY_BUFFER,sizeof(float)*NCT.size(),NCT.data(),GL_DYNAMIC_DRAW);
      glNormalPointer(GL_FLOAT,0,0);
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,ibo);
      glDrawElements(GL_TRIANGLES,T.size(),GL_UNSIGNED_INT,0);
      glBindBuffer(GL_ARRAY_BUFFER,0);
    }

  stale = false;
}
