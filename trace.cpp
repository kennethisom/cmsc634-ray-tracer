#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <string>
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

#define IMG_WIDTH 0
#define IMG_HEIGHT 1
#define COLOR_R 0
#define COLOR_G 1
#define COLOR_B 2
#define PI 3.14159265358979323846
#define BIAS 0.001
        
void tokenize(string& line, vector<string>& tokens) {
    string token;
    size_t start_pos, space_pos, num_chars;
    
    start_pos = 0;
    space_pos = 0;
        
    while ((space_pos = line.find(" ", start_pos)) != string::npos) {
        if (start_pos == space_pos) { // Character at start_pos is a space. Skip.
            start_pos++;
        } else {
            num_chars = space_pos - start_pos;
            token = line.substr(start_pos, num_chars);
            if (token != "" && token.find(" ") == string::npos) { // Skips tokens that are empty or all spaces.
                tokens.push_back(token);
            }
            start_pos = space_pos + 1;
        }
    }
    
    // If line doesn't end with a space there will be one last token to grab
    if (start_pos < line.size()) {
        token = line.substr(start_pos, string::npos);
        tokens.push_back(token);
    }
}

class Vector3 {

    private:
        float x;
        float y;
        float z;
        
    public:
        Vector3() {
            x = 0.0;
            y = 0.0;
            z = 0.0;
        }
        
        Vector3(float x_val, float y_val, float z_val) {
            x = x_val;
            y = y_val;
            z = z_val;
        }
        
        void set(float x_val, float y_val, float z_val) {
            x = x_val;
            y = y_val;
            z = z_val;
        }
        
        float get_x() {
            return x;
        }
        
        float get_y() {
            return y;
        }
        
        float get_z() {
            return z;
        }
        
        float get_index(int index) {
            float value = 0.0;
            switch(index) {
                case 0:
                    value = x;
                    break;
                case 1:
                    value = y;
                    break;
                case 2:
                    value = z;
                    break;
                default:
                    //exit with error
                    break;
            }
            
            return value;
        }
        
        Vector3* plus_vector(Vector3* term2) {
            float result_x = x + term2->get_x();
            float result_y = y + term2->get_y();
            float result_z = z + term2->get_z();
            return new Vector3(result_x, result_y, result_z);
        }
        
        Vector3* minus_vector(Vector3* term2) {
            float result_x = x - term2->get_x();
            float result_y = y - term2->get_y();
            float result_z = z - term2->get_z();
            return new Vector3(result_x, result_y, result_z);
        }
        
        Vector3* scale(float scalar) {
            return new Vector3(scalar * x, scalar * y, scalar * z);
        }
        
        Vector3* cross_product(Vector3* term2) {
            float result_x = (y * term2->get_z()) - (z * term2->get_y());
            float result_y = (z * term2->get_x()) - (x * term2->get_z());
            float result_z = (x * term2->get_y()) - (y * term2->get_x());
            return new Vector3(result_x, result_y, result_z);
        }
        
        float dot_product(Vector3* term2) {
            float dp = 0;
            dp += x * term2->get_x();
            dp += y * term2->get_y();
            dp += z * term2->get_z();
            
            return dp;
        }
        
        float magnitude() {
            return sqrt((x*x) + (y*y) + (z*z));
        }
        
        Vector3* normalize() {
            float m = magnitude();
            return new Vector3(x / m, y / m, z / m);
        }
        
        Vector3* operator=(Vector3* rhs) {
            set(rhs->get_x(), rhs->get_y(), rhs->get_z());
            return this;
        }
};

class Light {

    private:
        Vector3* position; // Px Py Pz
        float color[3]; // Cx Cy Cz

    public:
        Light() {
            position = new Vector3(0.0, 0.0, 0.0);
            color[COLOR_R] = 1.0;
            color[COLOR_G] = 1.0;
            color[COLOR_B] = 1.0;
        }
        
        ~Light() {
            delete position;
        }
        
        void set_position(float x, float y, float z) {
            position->set(x, y, z);
        }
        
        void set_color(float r, float g, float b) {
            color[COLOR_R] = r;
            color[COLOR_G] = g;
            color[COLOR_B] = b;
        }
        
        void scale_intensity(int num_lights) {
            color[COLOR_R] *= (float) (1 / sqrt(num_lights));
            color[COLOR_G] *= (float) (1 / sqrt(num_lights));
            color[COLOR_B] *= (float) (1 / sqrt(num_lights));
        }
        
        Vector3* get_position() {
            return new Vector3(position->get_x(), position->get_y(), position->get_z());
        }
        
        float get_color(int index) {
            return color[index];
        }
        
        void print() {
            printf("l %g %g %g %g %g %g\n", position->get_x(), position->get_y(), position->get_z(), color[COLOR_R], color[COLOR_G], color[COLOR_B]);
        }
};

class Fill {

    private:
        float color[3]; // Cx Cy Cz
        float diffuse;
        float specular;
        float shine;
        float transmittance;
        float refraction;
        
    public:
        Fill() {
            color[COLOR_R] = 0.0;
            color[COLOR_G] = 0.0;
            color[COLOR_B] = 0.0;
            diffuse = 0.0;
            specular = 0.0;
            shine = 0.0;
            transmittance = 0.0;
            refraction = 0.0;
        }
        
        void set_fill(float r, float g, float b, float kd, float ks, float s, float kt, float ir) {
            color[COLOR_R] = r;
            color[COLOR_G] = g;
            color[COLOR_B] = b;
            diffuse = kd;
            specular = ks;
            shine = s;
            transmittance = kt;
            refraction = ir;
        }
        
        void get_fill(float& r, float& g, float& b, float& kd, float& ks, float& s, float& kt, float& ir) {
            r = color[COLOR_R];
            g = color[COLOR_G];
            b = color[COLOR_B];
            kd = diffuse;
            ks = specular;
            s = shine;
            kt = transmittance;
            ir = refraction;
        }
        
        void copy_color(float (&color_copy)[3]) {
            color_copy[COLOR_R] = color[COLOR_R];
            color_copy[COLOR_G] = color[COLOR_G];
            color_copy[COLOR_B] = color[COLOR_B];
        }
        
        void print() {
            printf("f %g %g %g %g %g %g %g %g\n", color[COLOR_R], color[COLOR_G], color[COLOR_B], diffuse, specular, shine, transmittance, refraction);
        }

};

class Primitive {

    protected:
        Fill* fill;
    
    public:
        ~Primitive() {
            //delete fill;
        }
        
        virtual void set_primitive(vector<string> tokens, Fill* f) {
            fill = f;
            return;
        };
        
        virtual float calculate_t(Vector3* v_e, Vector3* v_d, float t_min, float t_max) {
            return -1.0;
        };
        
        void calculate_color(Vector3* v_e, Vector3* v_d, float t, Vector3* nv_n, int iteration, float current_ir, float (&background)[3], vector<Primitive*>& primitives, vector<Light*>& lights, float (&color)[3]) {
            float r, g, b, kd, ks, s, kt, ir, diffuse, specular, lt, ltmin, ltmax, rt, rt_closest, rtmin;
            float reflection_color[3], refraction_color[3];
            Vector3* v_p;
            
            fill->get_fill(r, g, b, kd, ks, s, kt, ir);
            
            Vector3* v_td;
            v_td = v_d->scale(t);
            v_p = v_e->plus_vector(v_td);
            delete v_td;
            
            color[COLOR_R] = 0.0;
            color[COLOR_G] = 0.0;
            color[COLOR_B] = 0.0;
            
            for (int i=0; i<(int)lights.size(); i++) {
                Vector3* v_light;
                v_light = lights[i]->get_position();
                Vector3* v_l;
                v_l = v_light->minus_vector(v_p);
                delete v_light;
                Vector3* nv_l;
                nv_l = v_l->normalize();
                
                // Is this light visible?  If not continue
                ltmin = BIAS;
                ltmax = 1.0;
                bool visible = true;
                for (int j=0; j<(int)primitives.size(); j++) {
                    lt = primitives[j]->calculate_t(v_p, v_l, ltmin, ltmax);
                    if (lt >= ltmin && lt <= ltmax) {
                        visible = false;
                        break;
                    }
                }
                delete v_l;
                
                if (!visible) { continue; }
                
                Vector3* v_v;
                v_v = v_e->minus_vector(v_p);
                Vector3* nv_v;
                nv_v = v_v->normalize();
                delete v_v;
                Vector3* v_h;
                v_h = nv_l->plus_vector(nv_v);
                delete nv_v;
                Vector3* nv_h;
                nv_h = v_h->normalize();
                delete v_h;
                
                diffuse = max((float) 0.0, nv_n->dot_product(nv_l));
                delete nv_l;
                specular = pow(max((float) 0.0, nv_n->dot_product(nv_h)), s);
                delete nv_h;
                
                color[COLOR_R] += ((kd * r) + (ks * specular)) * diffuse * lights[i]->get_color(COLOR_R);
                color[COLOR_G] += ((kd * g )+ (ks * specular)) * diffuse * lights[i]->get_color(COLOR_G);
                color[COLOR_B] += ((kd * b) + (ks * specular)) * diffuse * lights[i]->get_color(COLOR_B);
            }
            
            // Calculate Color that Results from Reflection
            if (iteration < 5 && ks > 0.0) {
                reflection_color[COLOR_R] = 0.0;
                reflection_color[COLOR_G] = 0.0;
                reflection_color[COLOR_B] = 0.0;
                Vector3* nv_d;
                nv_d = v_d->normalize();
                Vector3* v_r;
                v_r = nv_d->minus_vector(nv_n->scale(2 * nv_d->dot_product(nv_n)));
                delete nv_d;
                Vector3* nv_r;
                nv_r = v_r->normalize();
                delete v_r;
                
                rt = -1;
                rtmin = BIAS;
                rt_closest = INFINITY;
                int primitive_i = -1;
                for (int i=0; i<(int)primitives.size(); i++) {
                    rt = primitives[i]->calculate_t(v_p, nv_r, rtmin, INFINITY);
                    if (rt > rtmin && rt < rt_closest) {
                        rt_closest = rt;
                        primitive_i = i;
                    }
                }
                    
                if (primitive_i != -1) {
                    primitives[primitive_i]->get_color(v_p, nv_r, rt_closest, iteration+1, current_ir, background, primitives, lights, reflection_color);
                    color[COLOR_R] += ks * reflection_color[COLOR_R];
                    color[COLOR_G] += ks * reflection_color[COLOR_G];
                    color[COLOR_B] += ks * reflection_color[COLOR_B];
                } else {
                    color[COLOR_R] += ks * background[COLOR_R];
                    color[COLOR_G] += ks * background[COLOR_G];
                    color[COLOR_B] += ks * background[COLOR_B];
                }
                delete nv_r;
            }
            
            // Calculate Color that Results from Refraction
            if (iteration < 5 && kt > 0.0) {
                refraction_color[COLOR_R] = 0.0;
                refraction_color[COLOR_G] = 0.0;
                refraction_color[COLOR_B] = 0.0;
                
                bool refracted = true;
                
                Vector3* nv_d;
                nv_d = v_d->normalize();
                float dp_dn = nv_d->dot_product(nv_n);
                
                Vector3* v_t = new Vector3();
                if (current_ir == ir) {
                    delete v_t;
                    v_t = new Vector3(v_d->get_x(), v_d->get_y(), v_d->get_z());
                } else {
                    float cosphi = 1 - (((current_ir*current_ir)/(ir*ir)) * (1 - (dp_dn * dp_dn)));
                    if (cosphi >= 0.0) {
                        Vector3* v_ncostheta;
                        v_ncostheta = nv_n->scale(dp_dn);
                        Vector3* v_ncosthetad;
                        v_ncosthetad = nv_d->minus_vector(v_ncostheta);
                        delete v_ncostheta;
                        float ndnt = current_ir / ir;
                        Vector3* v_bsintheta;
                        v_bsintheta = v_ncosthetad->scale(ndnt);
                        delete v_ncosthetad;
                        cosphi = sqrt(cosphi);
                        Vector3* v_ncosphi;
                        v_ncosphi = nv_n->scale(cosphi);
                        delete v_t;
                        v_t = v_bsintheta->minus_vector(v_ncosphi);
                        delete v_ncosphi;
                        delete v_bsintheta;
                    } else {
                        refracted = false;
                    }
                }
                
                if (refracted) {
                    Vector3* nv_t;
                    nv_t = v_t->normalize();
                    
                    rt = -1;
                    rtmin = BIAS;
                    rt_closest = INFINITY;
                    int primitive_i = -1;
                    for (int i=0; i<(int)primitives.size(); i++) {
                        rt = primitives[i]->calculate_t(v_p, nv_t, rtmin, INFINITY);
                        if (rt > rtmin && rt < rt_closest) {
                            rt_closest = rt;
                            primitive_i = i;
                        }
                    }
                        
                    if (primitive_i != -1) {
                        primitives[primitive_i]->get_color(v_p, nv_t, rt_closest, iteration+1, ir, background, primitives, lights, refraction_color);
                        color[COLOR_R] += kt * refraction_color[COLOR_R];
                        color[COLOR_G] += kt * refraction_color[COLOR_G];
                        color[COLOR_B] += kt * refraction_color[COLOR_B];
                    } else {
                        color[COLOR_R] += kt * background[COLOR_R];
                        color[COLOR_G] += kt * background[COLOR_G];
                        color[COLOR_B] += kt * background[COLOR_B];
                    }
                    delete nv_t;
                }
                delete nv_d;
                delete v_t;
            }
            
            delete v_p;
        }
        
        virtual void get_color(Vector3* v_e, Vector3* v_d, float t, int iteration, float current_ir, float (&background)[3], vector<Primitive*>& primitives, vector<Light*>& lights, float (&color)[3]) {
            fill->copy_color(color);
        }
        
        virtual void print() {
            return;
        };

};

class Cone : public Primitive {

    private:
        Vector3* base_position; // BPx BPy BPz
        Vector3* apex_position; // APx APy APz
        float base_radius;
        float apex_radius;
    
    public:
        ~Cone() {
            delete base_position;
            delete apex_position;
        }
        
        void set_primitive(vector<string> tokens, Fill* f) {
            fill = f;
            
            if (tokens.size() > 1) {
                base_position = new Vector3((float) atof(tokens[1].c_str()), (float) atof(tokens[2].c_str()), (float) atof(tokens[3].c_str()));
                base_radius = (float) atof(tokens[4].c_str());
                
                apex_position = new Vector3((float) atof(tokens[5].c_str()), (float) atof(tokens[6].c_str()), (float) atof(tokens[7].c_str()));
                apex_radius = (float) atof(tokens[8].c_str());
            } else {
                string line;
                
                getline(cin, line);
                tokens.clear();
                tokenize(line, tokens);
                
                base_position = new Vector3((float) atof(tokens[0].c_str()), (float) atof(tokens[1].c_str()), (float) atof(tokens[2].c_str()));
                base_radius = (float) atof(tokens[3].c_str());
                
                getline(cin, line);
                tokens.clear();
                tokenize(line, tokens);
                
                apex_position = new Vector3((float) atof(tokens[0].c_str()), (float) atof(tokens[1].c_str()), (float) atof(tokens[2].c_str()));
                apex_radius = (float) atof(tokens[3].c_str());
            }
        };
        
        float calculate_t(Vector3* v_e, Vector3* v_d, float t_min, float t_max) {
            return -1.0;
        };
        
        void get_color(Vector3* v_e, Vector3* v_d, float t, int iteration, float current_ir, float (&background)[3], vector<Primitive*>& primitives, vector<Light*>& lights, float (&color)[3]) {
            fill->copy_color(color);
        }
        
        void print() {
            fill->print();
            printf("c\n");
            printf("%g %g %g %g\n", base_position->get_x(), base_position->get_y(), base_position->get_z(), base_radius);
            printf("%g %g %g %g\n", apex_position->get_x(), apex_position->get_y(), apex_position->get_z(), apex_radius);
            
        };
        
};

class Sphere : public Primitive {

    private:
        Vector3* center; // BPx BPy BPz
        float radius;
    
    public:
        ~Sphere() {
            delete center;
        }
        
        void set_primitive(vector<string> tokens, Fill* f) {
            fill = f;
            center = new Vector3((float) atof(tokens[1].c_str()), (float) atof(tokens[2].c_str()), (float) atof(tokens[3].c_str()));
            radius = (float) atof(tokens[4].c_str());
        };
        
        float calculate_t(Vector3* v_e, Vector3* v_d, float t_min, float t_max) {
            float t, dp_dce, dp_dd, dp_cece;
            float discriminant;
            Vector3* v_ce;
            
            v_ce = v_e->minus_vector(center);
            
            dp_dce = v_d->dot_product(v_ce);
            dp_dd = v_d->dot_product(v_d);
            dp_cece = v_ce->dot_product(v_ce);
            discriminant = (dp_dce * dp_dce) - (dp_dd * (dp_cece - (radius * radius)));
            
            if (discriminant < 0) {
                t = -1.0;
            } else {
                t = (-dp_dce - sqrt(discriminant)) / dp_dd;
                if (t < t_min) {
                    t = (-dp_dce + sqrt(discriminant)) / dp_dd;
                }
            }
            
            delete v_ce;
            
            return t;
        };
        
        void get_color(Vector3* v_e, Vector3* v_d, float t, int iteration, float current_ir, float (&background)[3], vector<Primitive*>& primitives, vector<Light*>& lights, float (&color)[3]) {
            Vector3* v_td;
            v_td = v_d->scale(t);
            Vector3* v_p;
            v_p = v_e->plus_vector(v_td);
            delete v_td;
            Vector3* v_n;
            v_n = v_p->minus_vector(center);
            delete v_p;
            Vector3* nv_n;
            nv_n = v_n->normalize();
            delete v_n;
            
            if(in_sphere(v_e)) {
                Vector3* inv_n;
                inv_n = nv_n->scale(-1.0);
                delete nv_n;
                nv_n = new Vector3(inv_n->get_x(), inv_n->get_y(), inv_n->get_z());
                delete inv_n;
            }
            
            calculate_color(v_e, v_d, t, nv_n, iteration, current_ir, background, primitives, lights, color);
            delete nv_n;
        }
        
        bool in_sphere(Vector3* v_e) {
            float t1, t2, dp_dce, dp_dd, dp_cece;
            float discriminant;
            Vector3* v_d;
            Vector3* v_ce;
            
            bool in = false;
            
            v_ce = v_e->minus_vector(center);
            v_d = new Vector3(1, 0, 0);
            
            dp_dce = v_d->dot_product(v_ce);
            dp_dd = v_d->dot_product(v_d);
            dp_cece = v_ce->dot_product(v_ce);
            discriminant = (dp_dce * dp_dce) - (dp_dd * (dp_cece - (radius * radius)));
            
            if (discriminant > 0) {
                t1 = (-dp_dce - sqrt(discriminant)) / dp_dd;
                t2 = (-dp_dce + sqrt(discriminant)) / dp_dd;
                if (t1 < t2) {
                    if (t1 < -BIAS && t2 > BIAS) {
                        in = true;
                    }
                } else {
                    if (t2 < -BIAS && t1 > BIAS) {
                        in = true;
                    }
                }
                
            }
            
            delete v_ce;
            
            return in;
        }
        
        void print() {
            fill->print();
            printf("s %g %g %g %g\n", center->get_x(), center->get_y(), center->get_z(), radius);
        };
        
};

class Polygon : public Primitive {

    private:
        vector<Vector3*> vertices; // Vx Vy Vz
        Vector3* normal; // Nx Ny Nz
    
    public:
        ~Polygon() {
            for(int i=0; i<(int)vertices.size(); i++) {
                delete vertices[i];
            }
            delete normal;
        }
        
        void set_primitive(vector<string> tokens, Fill* f) {
            fill = f;
            string line;
            
            int num_vertices = atoi(tokens[1].c_str());
            for (int i=0; i<num_vertices; i++) {
                getline(cin, line);
                tokens.clear();
                tokenize(line, tokens);
                
                vertices.push_back(new Vector3((float) atof(tokens[0].c_str()), (float) atof(tokens[1].c_str()), (float) atof(tokens[2].c_str())));
            }
            
            
            Vector3* v_p0p1;
            Vector3* v_p0p2;
            v_p0p1 = vertices[1]->minus_vector(vertices[0]);
            v_p0p2 = vertices[2]->minus_vector(vertices[0]);
            Vector3* v_n;
            v_n = v_p0p1->cross_product(v_p0p2);
            delete v_p0p1;
            delete v_p0p2;
            normal = v_n->normalize();
            delete v_n;
            
        };
        
        float calculate_t(Vector3* v_e, Vector3* v_d, float t_min, float t_max) {
            float dp_ep0n, dp_dn;
            float t;
            
            Vector3* v_ep0;
            v_ep0 = vertices[0]->minus_vector(v_e);
            
            dp_ep0n = v_ep0->dot_product(normal);
            delete v_ep0;
            dp_dn = v_d->dot_product(normal);
            
            if (dp_dn == 0) {
                t = -1.0;
            } else {
                t = dp_ep0n / dp_dn;
            }
            
            Vector3* v_td;
            v_td = v_d->scale(t);
            Vector3* v_p;
            v_p = v_e->plus_vector(v_td);
            delete v_td;
            
            // Compute p1 - v_e as v_ep1
            // Compute dp_ep1n
            // Compute dp_dn
            // Compute t = dp_ep1n / dp_dn
            // Compute v_p = v_e + t * v_d
            
            // Determine if v_p is inside or outside the polygon
            // CrossedEdges = 0
            // For all Consecutive Vertice Pairs:
            //  E = Vi+1 - Vi
            //  R = [1, 0]
            //  Find j - the scalar along Vi + jE where P + kR will intersect it.
            //  Find k - the scalar along P + kR where Vi + jE will intersect it.
            //  If k > 0 and k < INFINITY and j > 0 and j < 1 then intersection occured, increase CrossedEdges
            
            // Given P + kR = Vi + jE => jE - kR = P - Vi
            //
            // And |a b| |j| = |c|
            //     |d e| |k|   |f|
            //
            // Then
            // j = (ce - bf)/(ea - bd)
            // k = (af - dc)/(ae - db)
            //
            // Since b = -1 and e = 0
            // j = f/d
            // k = (af - dc)/d
            
            float a, b, c, d, e, f, j, k;
            float v_e2[2], nv_r[2], v_v0[2], v_v1[2];
            int index0, index1;
            
            if (abs(normal->get_z()) > abs(normal->get_x()) && abs(normal->get_z()) > abs(normal->get_y())) {
                index0 = 0;
                index1 = 1;
            } else if (abs(normal->get_y()) > abs(normal->get_x())) {
                index0 = 0;
                index1 = 2;
            } else {
                index0 = 1;
                index1 = 2;
            }
            
            nv_r[0] = 1.0;
            nv_r[1] = 0.0;
            
            int CrossedEdges = 0;
            for (int i=0; i < (int)vertices.size(); i++) {
                v_v0[0] = vertices[i]->get_index(index0);
                v_v0[1] = vertices[i]->get_index(index1);
                if (i+1 == (int)vertices.size()) {
                    v_v1[0] = vertices[0]->get_index(index0);
                    v_v1[1] = vertices[0]->get_index(index1);
                } else {
                    v_v1[0] = vertices[i+1]->get_index(index0);
                    v_v1[1] = vertices[i+1]->get_index(index1);
                }
                v_e2[0] = v_v1[0] - v_v0[0];
                v_e2[1] = v_v1[1] - v_v0[1];
                
                a = v_e2[0];
                b = -nv_r[0];
                c = v_p->get_index(index0) - v_v0[0];
                d = v_e2[1];
                e = -nv_r[1];
                f = v_p->get_index(index1) - v_v0[1];
                
                j = f / d;
                k = ((a * f) - (d * c)) / d;
                
                if (k > 0 && k < INFINITY && j > 0 && j < 1) {
                    CrossedEdges++;
                }
            }
            delete v_p;
            
            if (CrossedEdges % 2 == 0) {
                t = -1.0;
            }
            
            return t;
        };
        
        void get_color(Vector3* v_e, Vector3* v_d, float t, int iteration, float current_ir, float (&background)[3], vector<Primitive*>& primitives, vector<Light*>& lights, float (&color)[3]) {
            calculate_color(v_e, v_d, t, normal, iteration, current_ir, background, primitives, lights, color);
        };
        
        void print() {
            fill->print();
            printf("p %d\n", (int) vertices.size());
            for (int i=0; i<(int)vertices.size(); i++) {
                printf("%g %g %g\n", vertices[i]->get_x(), vertices[i]->get_y(), vertices[i]->get_z());
            }
        };

};

class Patch : public Primitive {

    private:
        vector<Vector3*> vertices; // Vx Vy Vz
        vector<Vector3*> normals; // Nx Ny Nz
    
    public:
        ~Patch() {
            for(int i=0; i<(int)vertices.size(); i++) {
                delete vertices[i];
                delete normals[i];
            }
            
            //vertices.clear();
            //normals.clear();
        }
        
        void set_primitive(vector<string> tokens, Fill* f) {
            fill = f;
            string line;
            
            int num_vertices = atoi(tokens[1].c_str());
            for (int i=0; i<num_vertices; i++) {
                getline(cin, line);
                tokens.clear();
                tokenize(line, tokens);
                
                vertices.push_back(new Vector3((float) atof(tokens[0].c_str()), (float) atof(tokens[1].c_str()), (float) atof(tokens[2].c_str())));
                normals.push_back(new Vector3((float) atof(tokens[3].c_str()), (float) atof(tokens[4].c_str()), (float) atof(tokens[5].c_str())));
            }
        };
        
        float calculate_t(Vector3* v_e, Vector3* v_d, float t_min, float t_max) {
            return -1.0;
        };
        
        void get_color(Vector3* v_e, Vector3* v_d, float t, int iteration, float current_ir, float (&background)[3], vector<Primitive*>& primitives, vector<Light*>& lights, float (&color)[3]) {
            fill->copy_color(color);
        };
        
        void print() {
            fill->print();
            printf("pp %d\n", (int) vertices.size());
            for (int i=0; i<(int)vertices.size(); i++) {
                printf("%g %g %g %g %g %g\n", vertices[i]->get_x(), vertices[i]->get_y(), vertices[i]->get_z(), normals[i]->get_x(), normals[i]->get_y(), normals[i]->get_z());
            }
        };

};

class Scene {

    private:
        Vector3* from; // Fx Fy Fz
        Vector3* at; // Ax Ay Az
        Vector3* up; // Ux Uy Uz
        float angle;
        float hither;
        int resolution[2]; //Rx Ry
        float background[3]; // Br Bg Bb
        vector<Light*> lights;
        vector<Primitive*> primitives;

    public:
        Scene() {
            from = new Vector3();
            at = new Vector3();
            up = new Vector3();
            angle = 0.0;
            hither = 0.0;
            resolution[IMG_WIDTH] = 0.0;
            resolution[IMG_HEIGHT] = 0.0;
            background[COLOR_R] = 0.0;
            background[COLOR_G] = 0.0;
            background[COLOR_B] = 0.0;
            lights.clear();
            primitives.clear();
        }
        
        ~Scene() {
            delete from;
            delete at;
            delete up;
            
            for(int i=0; i<(int)lights.size(); i++) {
                delete lights[i];
            }
            for(int i=0; i<(int)primitives.size(); i++) {
                delete primitives[i];
            }
            
            //lights.clear();
            //primitives.clear();
        }
        
        void set_background(float r, float g, float b) {
            background[COLOR_R] = r;
            background[COLOR_G] = g;
            background[COLOR_B] = b;
        }
        
        void set_view() {
            string line;
            vector<string> tokens;
            
            for (int i=0; i<6; i++) {
                getline(cin, line);
                tokens.clear();
                tokenize(line, tokens);
                
                if (tokens[0] == "from") {
                    from->set((float) atof(tokens[1].c_str()), (float) atof(tokens[2].c_str()), (float) atof(tokens[3].c_str()));
                } else if (tokens[0] == "at") {
                    at->set((float) atof(tokens[1].c_str()), (float) atof(tokens[2].c_str()), (float) atof(tokens[3].c_str()));
                } else if (tokens[0] == "up") {
                    up->set((float) atof(tokens[1].c_str()), (float) atof(tokens[2].c_str()), (float) atof(tokens[3].c_str()));
                } else if (tokens[0] == "angle") {
                    angle = (float) atof(tokens[1].c_str());
                } else if (tokens[0] == "hither") {
                    hither = (float) atof(tokens[1].c_str());
                } else if (tokens[0] == "resolution") {
                    resolution[IMG_WIDTH] = atoi(tokens[1].c_str());
                    resolution[IMG_HEIGHT] = atoi(tokens[2].c_str());
                } else {
                    // Error
                }
            }
        }

        void parse_input() {
            Fill* fill = NULL;
            Light* light = NULL;
            Primitive* primitive;
            string line;
            vector<string> tokens;
            
            while (getline(cin, line)) {
                if ((int)line.size() == 0) {
                    continue;
                }
                
                // Tokenize String
                tokens.clear();
                tokenize(line, tokens);
                
                // Process Tokens
                switch(tokens[0][0]) {
                    case 'b':
                        set_background((float) atof(tokens[1].c_str()), (float) atof(tokens[2].c_str()), (float) atof(tokens[3].c_str()));
                        break;
                    case 'v':
                        set_view();
                        break;
                    case 'l':
                        light = new Light();
                        light->set_position((float) atof(tokens[1].c_str()), (float) atof(tokens[2].c_str()), (float) atof(tokens[3].c_str()));
                        if (tokens.size() == 7) {
                            light->set_color((float) atof(tokens[4].c_str()), (float) atof(tokens[5].c_str()), (float) atof(tokens[6].c_str()));
                        }
                        lights.push_back(light);
                        break;
                    case 'f':
                        fill = new Fill();
                        fill->set_fill((float) atof(tokens[1].c_str()), (float) atof(tokens[2].c_str()), (float) atof(tokens[3].c_str()),
                                       (float) atof(tokens[4].c_str()), (float) atof(tokens[5].c_str()), (float) atof(tokens[6].c_str()),
                                       (float) atof(tokens[7].c_str()), (float) atof(tokens[8].c_str()));
                        break;
                    case 'c':
                        primitive = new Cone;
                        primitive->set_primitive(tokens, fill);
                        primitives.push_back(primitive);
                        break;
                    case 's':
                        primitive = new Sphere;
                        primitive->set_primitive(tokens, fill);
                        primitives.push_back(primitive);
                        break;
                    case 'p':
                        //if (tokens[0] == "pp") {
                        //    primitive = new Patch;
                        //    primitive->set_primitive(tokens, fill);
                        //    primitives.push_back(primitive);
                        //} else {
                            primitive = new Polygon;
                            primitive->set_primitive(tokens, fill);
                            primitives.push_back(primitive);
                        //}
                        break;
                    case '#':
                        // Do nothing
                        break;
                    default:
                        // Exit with error
                        break;
                }
            }
            
            // Scale Light Intensity
            if ((int)lights.size() > 0) {
                for (int i=0; i<(int)lights.size(); i++) {
                    lights[i]->scale_intensity((int)lights.size());
                }
            }
            
            //delete fill;
            //delete light;
            //delete primitive;
            tokens.clear();
        }
        
        void render_image(float aperture_size, int num_rays) {
            int width, height;
            int primitive_i;
            float color[3];
            float pixel_color[3];
            float angle_rad;
            float t, closest;
            float ic_b, ic_u, ic_l, ic_r;
            float vs, us;
            float d;
            
            //Calculate vE, nvW, nvU, nvV
            Vector3* v_w;
            v_w = from->minus_vector(at);
            d = v_w->magnitude();
            Vector3* nv_w;
            nv_w = v_w->normalize();
            delete v_w;
            
            Vector3* v_u;
            v_u = up->cross_product(nv_w);
            Vector3* nv_u;
            nv_u = v_u->normalize();
            delete v_u;
            
            Vector3* v_v;
            v_v = nv_w->cross_product(nv_u);
            Vector3* nv_v;
            nv_v = v_v->normalize();
            delete v_v;
            
            angle_rad = angle * PI / 180;
            width = resolution[IMG_WIDTH];
            height = resolution[IMG_HEIGHT];
            
            unsigned char pixels[height][width][3];
            
            //Calculate u, b, r, l
            ic_u = tan(angle_rad/2) * d;
            ic_r = ic_u;
            ic_l = 0 - ic_r;
            ic_b = 0 - ic_u;
            
            int j_max = height - 1;
            for (int j=0; j<height; j++) {
                for (int i=0; i<width; i++) {
                    us = ic_l + (ic_r - ic_l) * (i + 0.5) / width;
                    vs = ic_b + (ic_u - ic_b) * (j + 0.5) / height;
                    
                    Vector3* v_usu;
                    v_usu = nv_u->scale(us);
                    Vector3* v_vsv;
                    v_vsv = nv_v->scale(vs);
                    
                    pixel_color[COLOR_R] = 0.0;
                    pixel_color[COLOR_G] = 0.0;
                    pixel_color[COLOR_B] = 0.0;
                    
                    for (int k=0; k<num_rays; k++) {
                        closest = INFINITY;
                        primitive_i = -1;
                    
                        color[COLOR_R] = background[COLOR_R];
                        color[COLOR_G] = background[COLOR_G];
                        color[COLOR_B] = background[COLOR_B];
                        
                        //Pick random vector within from aperture => v_astart
                        Vector3* v_astart;
                        if (aperture_size > 0.0) {
                            // pick random angle
                            // pick random magnitude
                            // aus = magnitude * cos(angle)
                            // avs = magnitude * sin(angle)
                            // v_astart = from + aus*nv_u
                            // v_astart = v_a + avs*nv_v
                            float a_angle, a_angle_rad, a_magnitude, a_us, a_vs;
                            a_magnitude = (float) (aperture_size * rand() / RAND_MAX);
                            a_angle = (float) (360 * rand() / RAND_MAX);
                            a_angle_rad = a_angle * PI / 180;
                            a_us = a_magnitude * cos(a_angle_rad);
                            a_vs = a_magnitude * sin(a_angle_rad);
                            Vector3* v_ausu;
                            v_ausu = nv_u->scale(a_us);
                            Vector3* v_avsv;
                            v_avsv = nv_v->scale(a_vs);
                            v_astart = from->plus_vector(v_ausu);
                            delete v_ausu;
                            v_astart = v_astart->plus_vector(v_avsv);
                            delete v_avsv;
                        } else {
                            v_astart = new Vector3(from->get_x(), from->get_y(), from->get_z());
                        }
                        
                        Vector3* v_dw;
                        v_dw = nv_w->scale(-d);
                        
                        Vector3* v_s;
                        v_s = from->plus_vector(v_usu);
                        v_s = v_s->plus_vector(v_vsv);
                        v_s = v_s->plus_vector(v_dw);
                        delete v_dw;
                        
                        Vector3* v_d;
                        v_d = v_s->minus_vector(v_astart);
                        delete v_s;
                        Vector3* nv_d;
                        nv_d = v_d->normalize();
                        delete v_d;
                        
                        for (int l=0; l<(int)primitives.size(); l++) {
                            t = primitives[l]->calculate_t(v_astart, nv_d, hither, INFINITY);
                            if (t > hither && t < closest) {
                                closest = t;
                                primitive_i = l;
                            }
                        }
                    
                        if (primitive_i != -1) {
                            primitives[primitive_i]->get_color(v_astart, nv_d, closest, 0, 1.0, background, primitives, lights, color);
                        }
                        
                        if (color[COLOR_R] > 1.0) {
                            color[COLOR_R] = 1.0;
                        }
                        if (color[COLOR_G] > 1.0) {
                            color[COLOR_G] = 1.0;
                        }
                        if (color[COLOR_B] > 1.0) {
                            color[COLOR_B] = 1.0;
                        }
                        pixel_color[COLOR_R] += color[COLOR_R];
                        pixel_color[COLOR_G] += color[COLOR_G];
                        pixel_color[COLOR_B] += color[COLOR_B];
                        delete v_astart;
                        delete nv_d;
                    }
                    delete v_usu;
                    delete v_vsv;
                    
                    pixels[j_max - j][i][COLOR_R] = pixel_color[COLOR_R] / num_rays * 255;
                    pixels[j_max - j][i][COLOR_G] = pixel_color[COLOR_G] / num_rays * 255;
                    pixels[j_max - j][i][COLOR_B] = pixel_color[COLOR_B] / num_rays * 255;
                }
                if (j%32 == 31) {
                    printf("Rendered %d Lines\n", j+1);
                }
            }
            
            delete nv_w;
            delete nv_u;
            delete nv_v;
            
            printf("Creating Image File\n");
            FILE *f = fopen("trace.ppm","wb");
            fprintf(f, "P6\n%d %d\n%d\n", width, height, 255);
            fwrite(pixels, 1, height*width*3, f);
            fclose(f);
        }
        
        void print() {
            printf("b %g %g %g\n", background[COLOR_R], background[COLOR_G], background[COLOR_B]);
            printf("v\n");
            printf("from %g %g %g\n", from->get_x(), from->get_y(), from->get_z());
            printf("at %g %g %g\n", at->get_x(), at->get_y(), at->get_z());
            printf("up %g %g %g\n", up->get_x(), up->get_y(), up->get_z());
            printf("angle %g\n", angle);
            printf("hither %g\n", hither);
            printf("resolution %d %d\n", resolution[IMG_WIDTH], resolution[IMG_HEIGHT]);
            for(int i=0; i<(int)lights.size(); i++) {
                lights[i]->print();
            }
            for(int i=0; i<(int)primitives.size(); i++) {
                primitives[i]->print();
            }
        }

};

void print_help() {
    printf("Ray Tracer Command Line Arguments\n\n");
    printf("  -a [value]\tRadius of the aperture. Float Value. Default = 0\n");
    printf("  -r [value]\tNumber of rays for each pixel. Integer Value. Default = 1\n\n");
}

int main(int argc, char *argv[]) {
    int num_rays = 1;
    float aperture_size = 0;
    
    for (int i=1; i<argc; i++) {
        if (argv[i][0] == '-') {
            switch (argv[i][1]) {
                case 'h':
                    print_help();
                    exit(EXIT_SUCCESS);
                    break;
                case 'a':
                    i++;
                    aperture_size = (float) atof(argv[i]);
                    if (aperture_size < 0) {
                        printf("ERROR: Aperture size must be greater than or equal to 0.\n");
                        exit(EXIT_FAILURE);
                    }
                    break;
                case 'r':
                    i++;
                    num_rays = atoi(argv[i]);
                    if (num_rays < 1) {
                        printf("ERROR: Aperture size must be a positive integer greater than 0.\n");
                        exit(EXIT_FAILURE);
                    }
                    break;
                default:
                    printf("ERROR: Incorrect Command Line Arguments Provided!\n\n");
                    print_help();
                    exit(EXIT_FAILURE);
            }
        } else {
            printf("ERROR: Incorrect Command Line Arguments Provided!\n\n");
            print_help();
            exit(EXIT_FAILURE);
        }
    }
    
    Scene* scene = new Scene();
    printf("Parsing Input... ");
    scene->parse_input();
    printf(" Done!\n");
    printf("Rendering Image...\n");
    scene->render_image(aperture_size, num_rays);
    printf("Image Rendering Complete!\n");
    
    delete scene;
    
    return 0;
}
