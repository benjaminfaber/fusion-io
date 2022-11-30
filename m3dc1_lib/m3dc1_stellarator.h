#ifndef M3DC1_STELLARATOR_H
#define M3DC1_STELLARATOR_H

#include "m3dc1_mesh.h"
#include "m3dc1_field.h"


class m3dc1_stellarator_mesh : public m3dc1_3d_mesh {

   protected:

    virtual void logical_to_global(const double xl,
                     const double Phi, const double zl,
                     double* R, double* Z) const;
    
    virtual void logical_to_local(int i, const double xl,
                     const double Phi, const double zl,
                     double* xi, double* zi, double* eta) const;

    virtual void global_to_local(const double Rst,
                     const double Phi, const double Zst,
                     double* xi, double* zi, double* eta);
  public:
    enum eq_type_enum {
        VMEC = 1,
    };

    double* rst;
    double* zst;
    m3dc1_stellarator_mesh::eq_type_enum eq_type;
    std::string eq_file;
    
    virtual void global_to_logical(const double Rst,
                     const double Phi, const double Zst,
                     double* xl, double* zl) const;


    virtual int in_element(double Rst, double Phi, double Zst,
                double* xi_out=0, double* zi_out=0, double* eta_out=0,
                int guess=-1)
    {
        double xl, zl;
        m3dc1_stellarator_mesh::global_to_logical(Rst, Phi, Zst, &xl, &zl);
        return m3dc1_3d_mesh::in_element(xl, Phi, zl, xi_out, zi_out, eta_out, guess);
    }

                     
    m3dc1_stellarator_mesh(int n, m3dc1_stellarator_mesh::eq_type_enum eq_type, std::string eq_file);
    ~m3dc1_stellarator_mesh();
};


class m3dc1_stellarator_field : public m3dc1_3d_field {
  public:
    m3dc1_3d_field* rst;
    m3dc1_3d_field* zst;
    
  public:
    m3dc1_stellarator_mesh* stell_mesh;
    m3dc1_stellarator_field(m3dc1_stellarator_mesh* m);
    virtual ~m3dc1_stellarator_field();

    virtual bool eval(const double r, const double phi, const double z,
              const m3dc1_field::m3dc1_get_op op, double* val,
              int* element=0);
};

#endif