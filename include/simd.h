#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>

struct LigandSIMD {
    int atoms;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
};

LigandSIMD read_xyz_simd(std::string file_name) {
    std::ifstream infile(file_name);
    
    LigandSIMD lig;
    if (!infile) {
        std::cerr << "Error: could not open file, atoms will be empty\n";
        return lig;
    }

    // First line is number of atoms
    infile >> lig.atoms;

    double x, y, z;
    while(infile >> x >> y >> z) {
        lig.x.push_back(x);
        lig.y.push_back(y);
        lig.z.push_back(z);
    }
    return lig;
}

LigandSIMD transform_ligand_simd(const LigandSIMD& ligand, int pose_index)
{
    LigandSIMD out;
    int n = ligand.atoms;
    out.x.resize(n);
    out.y.resize(n);
    out.z.resize(n);

    // -----------------------------------------------------------
    // 1) Compute ligand centroid
    // -----------------------------------------------------------
    double cx = 0, cy = 0, cz = 0;
    for (int i = 0; i < n; i++) {
        cx += ligand.x[i];
        cy += ligand.y[i];
        cz += ligand.z[i];
    }
    cx /= n;
    cy /= n;
    cz /= n;

    // -----------------------------------------------------------
    // 2) Generate Euler rotation matrix
    // -----------------------------------------------------------
    double alpha = (pose_index * 0.37);
    double beta  = (pose_index * 0.51);
    double gamma = (pose_index * 0.29);

    double R[3][3];
    euler_to_matrix(alpha, beta, gamma, R);

    // -----------------------------------------------------------
    // 3) Translation
    // -----------------------------------------------------------
    double tx = 5.0 * sin(pose_index * 0.21);
    double ty = 5.0 * cos(pose_index * 0.13);
    double tz = 5.0 * sin(pose_index * 0.17 + 0.5);

    // -----------------------------------------------------------
    // 4) Apply:   (p - centroid) → rotate → + centroid → + translation
    // -----------------------------------------------------------
    for (int i = 0; i < n; i++) {
        double x = ligand.x[i] - cx;
        double y = ligand.y[i] - cy;
        double z = ligand.z[i] - cz;

        double rx = R[0][0]*x + R[0][1]*y + R[0][2]*z;
        double ry = R[1][0]*x + R[1][1]*y + R[1][2]*z;
        double rz = R[2][0]*x + R[2][1]*y + R[2][2]*z;

        // Re-center and translate
        out.x[i] = rx + cx + tx;
        out.y[i] = ry + cy + ty;
        out.z[i] = rz + cz + tz;
    }

    return out;
}