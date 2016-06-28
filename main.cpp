#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#if defined(__APPLE__)
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
#include <windows.h>
#endif
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

const int MAX_DEPTH = 5; // raytrace recursion depth
const float EPSILON = 0.0001f; // float-comparsion error

// application window
const int screenWidth = 600;
const int screenHeight = 600;

struct Vector {
  float x, y, z;

  Vector() {
    x = y = z = 0;
  }

  Vector(float x0, float y0, float z0 = 0) {
    x = x0; y = y0; z = z0;
  }

  Vector operator+(float a) const {
    return Vector(x + a, y + a, z + a);
  }

  Vector operator*(float a) const {
    return Vector(x * a, y * a, z * a);
  }

  Vector operator/(float a) const {
    return Vector(x / a, y / a, z / a);
  }

  Vector operator+(const Vector& v) const  {
    return Vector(x + v.x, y + v.y, z + v.z);
  }

  Vector operator-(const Vector& v) const {
    return Vector(x - v.x, y - v.y, z - v.z);
  }

  Vector operator-() const {
    return Vector(-x, -y, -z);
  }

  // dot product
  float operator*(const Vector& v) const {
    return (x * v.x + y * v.y + z * v.z);
  }

  // cross product
  Vector operator%(const Vector& v) const {
    return Vector(y*v.z-z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);
  }

  float Length() const {
    return sqrt(x * x + y * y + z * z);
  }

  Vector normalize() const {
    float len = this->Length();
    return Vector(x/len, y/len, z/len);
  }
};

struct Color {
  float r, g, b;

  Color() {
    r = g = b = 0;
  }

  Color(float r0, float g0, float b0) {
    r = r0; g = g0; b = b0;
  }

  Color operator*(float a) {
    return Color(r * a, g * a, b * a);
  }

  Color operator+(float a) {
    return Color(r + a, g + a, b + a);
  }

  Color operator-(float a) {
    return Color(r - a, g - a, b - a);
  }

  Color operator*(const Color& c) {
    return Color(r * c.r, g * c.g, b * c.b);
  }

  Color operator/(float a) {
    return Color(r / a, g / a, b / a);
  }

  Color operator/(const Color& c) {
    return Color(r / c.r, g / c.g, b / c.b);
  }

  Color operator+(const Color& c) {
    return Color(r + c.r, g + c.g, b + c.b);
  }

  Color operator-(const Color& c) {
    return Color(r - c.r, g - c.g, b - c.b);
  }
};

struct Material {
  Color F0, ka, kd, ks, k, n;
  float shininess;
  bool isReflective, isRefractive;

  Material(Color ka, Color kd, Color ks, Color k, Color n,
      float shininess, bool isReflective, bool isRefractive) {
    this->ka = ka; this->kd = kd;
    this->ks = ks; this->n = n;
    this->shininess = shininess; this->isReflective = isReflective;
    this->isRefractive = isRefractive;

    F0 = ( (n - 1) * (n - 1) + k * k )
       / ( (n + 1) * (n + 1) + k * k);
  }

  Color Fresnel(Vector inDir, Vector normal) {
    float cosa = fabs( normal * inDir );
    return F0 + ( Color(1.0, 1.0, 1.0) - F0 ) * pow((1.0f - cosa), 5);
  }

  Vector reflect(Vector inDir, Vector normal) {
    return inDir - normal * ( normal * inDir ) * 2.0f;
  }

  Vector refract(Vector inDir, Vector normal) {
    float ior = n.r;
    float cosa = (normal * inDir) * (-1.0);

    if (cosa < 0) { cosa = -cosa; normal = normal * (-1.0); ior = 1.0 / n.r; }

    float disc = 1 - (1 - cosa * cosa) / ior / ior;
    if (disc < 0) { return reflect(inDir, normal); }

    return inDir/ior + normal * ( cosa/ior - sqrt(disc) );
  }

  Color shade(Vector normal, Vector viewDir, Vector lightDir, Color inRad) {
    float cosTheta = normal * lightDir;
    Color reflRad(0.0, 0.0, 0.0);

    if (cosTheta < 0) { return reflRad; }

    reflRad = inRad * kd * cosTheta;
    Vector halfway = (viewDir + lightDir).normalize();
    float cosDelta = normal * halfway;

    if (cosDelta < 0) { return reflRad; }
    return reflRad + inRad * ks * pow(cosDelta, shininess);
  }
};

Material gold( Color(0.25, 0.2, 0.075), Color(0.75, 0.61, 0.23),
               Color(0.83, 0.75, 0.56), Color(3.1, 2.7, 1.9),
               Color(0.17, 0.35, 1.5), 5.0f, true, false );
Material gold_diff( Color(0.25, 0.2, 0.075), Color(0.45, 0.31, 0.13),
              Color(0.83, 0.75, 0.56), Color(0, 0, 0),
              Color(0, 0, 0), 50.0f, false, false );
Material green( Color(0.1, 0.2, 0.1), Color(0.0, 0.4, 0.1),
                Color(0.1, 0.1, 0.1), Color(0, 0, 0),
                Color(0, 0, 0), 50.0f, false, false );
Material blue( Color(0.0, 0.0, 0.15), Color(0.05, 0.15, 0.42),
               Color(0.25, 0.35, 0.85), Color(0, 0, 0),
               Color(0, 0, 0), 25.0f, false, false );
Material red( Color(0.15, 0.0, 0.0), Color(0.3, 0.05, 0.05),
                Color(0.3, 0.3, 0.5), Color(0, 0, 0),
                Color(0, 0, 0), 25.0f, false, false );

struct Light {
  Vector position;
  Color output;

  Light(Vector position, Color output) {
    this->position = position;
    this->output = output;
  }

  Color getLightAt(Vector position) {
    float dist = (position - this->position).Length();
    return output / ( dist * dist * 4 * M_PI );
  }
};

struct Hit {
  float t;
  Vector position, normal;
  Material* material;
  bool blackHoleHit;

  Hit() { t = -1; blackHoleHit = false; }
};

struct Ray {
  Vector position, direction;

  Ray(Vector position = Vector(0,0,0), Vector direction = Vector(0,0,0)) {
    this->position = position;
    this->direction = direction;
  }
};

struct Intersectable {
  Material* material;

  Intersectable(Material* material) {
    this->material = material;
  }

  virtual Hit intersect(const Ray& ray) = 0;
};

struct Triangle : Intersectable {
  Vector a, b, c;
  Vector normal;

  Triangle(Vector a, Vector b, Vector c, Material* material)
      : Intersectable(material) {
    this->a = a; this->b = b; this->c = c;
    normal = ( (b-a) % (c-a) ).normalize();
  }

  Hit intersect(const Ray& ray) {
    Hit hit;

    float t = (a - ray.position) * normal / (ray.direction * normal);
    if (t < 0) { return hit; }
    Vector point = ray.position + ray.direction * t;

    bool isInside = ( (b-a) % (point-a) ) * normal > 0
            && ( (c-b) % (point-b) ) * normal > 0
            && ( (a-c) % (point-c) ) * normal > 0;

    if (isInside) {
      hit.t = t;
      hit.position = point;
      hit.normal = this->normal;
      hit.material = this->material;
    }

    return hit;
  }
};

struct Sphere : Intersectable {
  Vector position;
  float radius;

  Sphere(Vector position, float radius, Material* material)
      : Intersectable(material) {
    this->position = position;
    this->radius = radius;
  }

  Hit intersect(const Ray& ray) {
    Hit hit;
    Vector dist = ray.position - position;

    float a = (ray.direction * ray.direction);
    float b = (dist * ray.direction) * 2.0;
    float c = ( dist * dist ) - radius * radius;

    float discr = b * b - 4.0 * a * c;
    if (discr < 0) { return hit; }
    float sqrt_discr = sqrt(discr);

    float t1 = (-b + sqrt_discr)/2.0/a;
    float t2 = (-b - sqrt_discr)/2.0/a;

    if (t1 < EPSILON) { t1 = EPSILON * (-1.0f); }
    if (t2 < EPSILON) { t2 = EPSILON * (-1.0f); }

    if (t1 < 0.0f) { return hit; }

    float t = t2 > 0.0f ? t2 : t1;

    hit.t = t;
    hit.position = ray.position + ray.direction * t;
    hit.normal = (hit.position - position).normalize();
    hit.material = this->material;
    return hit;
  }
};

struct Torus : Intersectable {
  const static int U_PRECISION = 6;
  const static int V_PRECISION = 5;

  Vector position;
  Sphere* boundingSphere;
  float R, r;
  Triangle* triangles[ (U_PRECISION+1) * (V_PRECISION+1) * 2 ];
  int triangleCount;

  Torus(Vector position, float R, float r, Material* material)
      : Intersectable(material) {
    this->position = position; this->R = R; this->r = r;
    this->boundingSphere = new Sphere( position, R+r, NULL );
    this->triangleCount = 0;

    for (int i = 0; i <= U_PRECISION; ++i) {
      float u = (i / (float)U_PRECISION) * 2 * M_PI;
      float u_delta = 2 * M_PI / (float)U_PRECISION;

      for (int j = 0; j <= V_PRECISION; j++) {
        float v = (j / (float)V_PRECISION) * 2 * M_PI;
        float v_delta = 2 * M_PI / (float)V_PRECISION;

        Vector p1 = pointFromParams( u, v );
        Vector p2 = pointFromParams( u + u_delta, v );
        Vector p3 = pointFromParams( u, v + v_delta );
        Vector p4 = pointFromParams( u + u_delta, v + v_delta );

        addTriangle( new Triangle(p1, p2, p3, material) );
        addTriangle( new Triangle(p2, p4, p3, material) );
      }
    }
  }

  ~Torus() {
    delete boundingSphere;

    for (int i = 0; i < triangleCount; ++i) {
      delete triangles[i];
    }
  }

  void addTriangle(Triangle* t) {
    triangles[triangleCount] = t;
    triangleCount += 1;
  }

  Vector pointFromParams(float u, float v) {
    Vector p;
    p.x = R * cos(v) + r * cos(u) * cos(v) + position.x;
    p.y = R * sin(v) + r * cos(u) * sin(v) + position.y;
    p.z = r * sin(u) + position.z;
    return p;
  }

  Hit intersect(const Ray& ray) {
    Hit hit;
    hit.t = 99999.9f;
    bool didIntersect = false;

    Hit boundingHit = boundingSphere->intersect(ray);
    if (boundingHit.t < 0) { return hit; }

    for (int i = 0; i < triangleCount; ++i) {
      Hit trHit = triangles[i]->intersect(ray);

      if (trHit.t > EPSILON) {
        if (trHit.t < hit.t) {
          hit = trHit;
          didIntersect = true;
        }
      }
    }

    if (!didIntersect) { hit.t = -1.0f; }
    return hit;
  }
};

struct Plane : Intersectable {
  Vector position, normal;
  Material* material2;

  Plane(Vector position, Vector normal, Material* material, Material* material2)
      : Intersectable(material) {
    this->position = position;
    this->normal = normal;
    this->material2 = material2;
  }

  Hit intersect(const Ray& ray) {
    Hit hit;

    float d = normal * ray.direction;
    if ( d < EPSILON && d > (EPSILON * (-1.0)) ) { return hit; }

    float t = ( normal * (ray.position - position) / d ) * (-1.0);

    if (t > EPSILON) {
      hit.t = t;
      hit.position = ray.position + ray.direction * t;
      hit.normal = this->normal;
      hit.material = getMaterialAt( hit.position );
    }

    return hit;
  }

  Material* getMaterialAt(Vector pos) {
    const float mod = 0.2;
    const float cmp = 0.005;

    float x = fabs(fmod( pos.x, mod ));
    float y = fabs(fmod( pos.y, mod ));
    float z = fabs(fmod( pos.z, mod ));

    bool b = (( x < cmp && x > EPSILON )
           ^ ( y < cmp && y > EPSILON )
           ^ ( z < cmp && z > EPSILON ))
           || ( x < cmp && y < cmp )
           || ( y < cmp && z < cmp )
           || ( x < cmp && z < cmp );

    return ( b ? material : material2  );
  }
};

struct BlackHole {
  const static float earthMass = 0.6 * 10e25;
  const static float c = 3 * 10e8;
  const static float G = 6.67 * 10e-11;

  Vector position;
  float radius;

  BlackHole(Vector position = Vector(0,0,0)) {
    this->position = position;
    radius = 10 * 2 * G * earthMass / c / c;
  }
};

struct Camera {
  Vector eyepos, direction, up, right;
  float XM, YM;

  Camera( Vector eyepos = Vector(0, 0, -1.0),
          Vector direction = Vector(0, 0, 0),
          Vector up = Vector(0, 1, 0),
          Vector right = Vector(1, 0, 0),
          float XM = 600,
          float YM = 600) {
    this->eyepos = eyepos; this->direction = direction;
    this->up = up; this->right = right;
    this->XM = XM; this->YM = YM;
  }

  Ray getRay(int x, int y) {
    Ray ray;
    Vector p = direction + right * ( 2 * (float)x / XM - 1 )
                         + up    * ( 2 * (float)y / YM - 1 )
                         + Vector(0.5 / XM, 0.5 / YM);
    ray.position = eyepos;
    ray.direction = (p - eyepos).normalize();
    return ray;
  }
};

struct Scene {
  Color image[screenWidth*screenHeight];

  Intersectable* objects[15];
  int objCount;

  Light* lights[5];
  int lightCount;

  Color ambientLight;
  Camera camera;
  BlackHole blackHole;

  Scene() {
    objCount = 0;
    lightCount = 0;
    ambientLight = Color(1, 1, 1);
  }

  ~Scene() {
    for (int i = 0; i < objCount; ++i)
      delete objects[i];

    for (int i = 0; i < lightCount; ++i)
      delete lights[i];
  }

  void add(Intersectable* intersectable) {
    objects[ objCount ] = intersectable;
    objCount += 1;
  }

  void add(Light* light) {
    lights[ lightCount ] = light;
    lightCount += 1;
  }

  void build() {
    camera.eyepos = Vector(0.5f, 0.5f, 0);
    camera.direction = Vector(0.5f, 0.5f, 1);

    blackHole.position = Vector(0.3f, 0.3f, 0.6f);

    add( new Plane( Vector(0, 0, 0), Vector(0, 1, 0), &gold_diff, &red ) );
    add( new Plane( Vector(0, 1, 0), Vector(0, -1, 0), &gold_diff, &red ) );
    add( new Plane( Vector(0, 0, 0), Vector(1, 0, 0), &gold_diff, &blue ) );
    add( new Plane( Vector(1, 0, 0), Vector(-1, 0, 0), &gold_diff, &blue ) );
    add( new Plane( Vector(0, 0, 1), Vector(0, 0, -1), &gold_diff, &blue ) );

    add( new Torus( Vector(0.6f, 0.6f, 0.55f), 0.2f, 0.055f, &gold ) );

    add( new Light( Vector(0.4f, 0.6f, 0.45f), Color( 2, 2, 2 ) ) );
    add( new Light( Vector(0.7f, 0.1f, 0.1f), Color( 4, 4, 4 ) ) );
  }

  Hit firstIntersect(Ray ray) {
    Hit bestHit;

    for (int i = 0; i < objCount; ++i) {
      Hit hit = objects[i]->intersect(ray);
      if (hit.t > 0 && (bestHit.t < 0 || hit.t < bestHit.t)) {
        bestHit = hit;
      }
    }

    return bestHit;
  }

  Hit rayMarch(Ray ray) {
    const int MAX_CYCLES = 800;
    const float MARCH_STEP = 0.05;

    Hit hit;
    Ray nRay = ray;

    for (int i = 0; i < MAX_CYCLES; ++i) {
      hit = firstIntersect(nRay);

      if ((nRay.position - blackHole.position).Length() < blackHole.radius) {
        hit.blackHoleHit = true;
        return hit;
      }

      if ((hit.position - nRay.position).Length() < MARCH_STEP) {
        return hit;
      }

      float r = (blackHole.position - nRay.position).Length();
      Vector d = nRay.direction * MARCH_STEP
        + ( (blackHole.position - nRay.position) / r ) * 0.5
        * blackHole.radius * MARCH_STEP * MARCH_STEP / ( r * r );

      nRay.direction = d.normalize();
      nRay.position = nRay.position + nRay.direction * MARCH_STEP;
    }

    return hit;
  }

  Color trace(Ray ray, int depth) {
    if (depth > MAX_DEPTH) { return ambientLight; }
    Hit hit = rayMarch(ray);
    if (hit.t < 0) { return ambientLight; }
    if (hit.blackHoleHit) { return Color(0, 0, 0); }

    Color color = ambientLight * hit.material->ka;

    for (int i = 0; i < lightCount; ++i) {
      Ray shadowRay( hit.position, lights[i]->position - hit.position );
      Hit shadowHit = firstIntersect(shadowRay);

      if (shadowHit.t < 0 || (hit.position - shadowHit.position).Length()
          > (hit.position - lights[i]->position).Length() ) {
        Vector viewDir = ray.direction.normalize() * (-1.0);
        Vector lightDir = shadowRay.direction.normalize();

        color = color + hit.material->shade(hit.normal, lightDir,
               viewDir, lights[i]->getLightAt( hit.position ));
      }
    }

    Color fresnel = hit.material->Fresnel(hit.normal, ray.direction);

    if (hit.material->isReflective) {
      Ray reflected( hit.position,
        hit.material->reflect( ray.direction, hit.normal ) );

      color = color + fresnel * trace( reflected, depth + 1 );
    }

    if (hit.material->isRefractive) {
      Ray refracted( hit.position,
        hit.material->refract( ray.direction, hit.normal ) );

      color = color + fresnel * trace( refracted, depth + 1 );
    }

    return color;
  }

  void render() {
    for (int y = 0; y < screenHeight; y++) {
      for (int x = 0; x < screenWidth; x++) {
        Ray ray = camera.getRay(x, y);
        Color color = trace(ray, 0);
        color = color/(Color(1,1,1) + color);

        image[y * screenWidth + x] = color;
      }
      if (y % 8 == 0) {
        printf("Rendering progress: %.2f%%\n", 100*y/(float)screenHeight);
      }
    }
    glutPostRedisplay();
  }
} scene;

void onInitialization() {
  glViewport(0, 0, screenWidth, screenHeight);
  scene.build();
  scene.render();
}


void onDisplay() {
  glDrawPixels(screenWidth, screenHeight,
               GL_RGB, GL_FLOAT, scene.image);
  glutSwapBuffers();
}

// entry point
int main(int argc, char **argv) {
    glutInit(&argc, argv);
    glutInitWindowSize(600, 600); // 600x600 fixed window
    glutInitWindowPosition(100, 100);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);

    glutCreateWindow("Grafika hazi feladat");

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    onInitialization();
    glutDisplayFunc(onDisplay);

    glutMainLoop();

    return 0;
}
