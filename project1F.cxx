// For this project, you will remove hidden surfaces using the zbuffer algorithm.
// You will also add the ability to interpolate colors within a triangle.

#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <math.h>

#define NORMALS
using std::cerr;
using std::endl;

struct LightingParameters
{
    LightingParameters(void)
    {
         lightDir[0] = -0.6;
         lightDir[1] = 0;
         lightDir[2] = -0.8;
         Ka = 0.3;
         Kd = 0.7;
         Ks = 2.3;
         alpha = 2.5;
    };


    double lightDir[3]; // The direction of the light source
    double Ka;           // The coefficient for ambient lighting.
    double Kd;           // The coefficient for diffuse lighting.
    double Ks;           // The coefficient for specular lighting.
    double alpha;        // The exponent term for specular lighting.
};

LightingParameters lp;

class Matrix
{
  public:
    double          A[4][4];  // A[i][j] means row i, column j

    void            TransformPoint(const double *ptIn, double *ptOut);
    static Matrix   ComposeMatrices(const Matrix &, const Matrix &);
    void            Print(ostream &o);
};

void
Matrix::Print(ostream &o)
{
    for (int i = 0 ; i < 4 ; i++)
    {
        char str[256];
        sprintf(str, "(%.7f %.7f %.7f %.7f)\n", A[i][0], A[i][1], A[i][2], A[i][3]);
        o << str;
    }
}

Matrix
Matrix::ComposeMatrices(const Matrix &M1, const Matrix &M2)
{
    Matrix rv;
    for (int i = 0 ; i < 4 ; i++)
        for (int j = 0 ; j < 4 ; j++)
        {
            rv.A[i][j] = 0;
            for (int k = 0 ; k < 4 ; k++)
                rv.A[i][j] += M1.A[i][k]*M2.A[k][j];
        }

    return rv;
}

void
Matrix::TransformPoint(const double *ptIn, double *ptOut)
{
    ptOut[0] = ptIn[0]*A[0][0]
             + ptIn[1]*A[1][0]
             + ptIn[2]*A[2][0]
             + ptIn[3]*A[3][0];
    ptOut[1] = ptIn[0]*A[0][1]
             + ptIn[1]*A[1][1]
             + ptIn[2]*A[2][1]
             + ptIn[3]*A[3][1];
    ptOut[2] = ptIn[0]*A[0][2]
             + ptIn[1]*A[1][2]
             + ptIn[2]*A[2][2]
             + ptIn[3]*A[3][2];
    ptOut[3] = ptIn[0]*A[0][3]
             + ptIn[1]*A[1][3]
             + ptIn[2]*A[2][3]
             + ptIn[3]*A[3][3];
}

class Camera
{
  public:
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];

    Matrix          ViewTransform(void);
    Matrix          CameraTransform(void);
    Matrix          DeviceTransform(void);
};

// Matrix contains attribute:  double A[4][4] // A[i][j] means row i, column j
Matrix Camera::CameraTransform(void){
    // Calculate Camera Frame:
    //    - camera frame must be a basis: spans space... can get any point through a linear combination of basis vectors
    //    - every member must be linearly independent (no basis vector can be represented via the other basis vectors)
    //                  - linearly independent -> perpendicular vectors
    //  - Consists of (u,v,w,O) -> 4 vectors
    //          - O = origin = camera position (Camera.position[x,y,z])
    //          - w = O - focus = (Camera.position[x,y,z] - Camera.focus[x,y,z]) 
    //          - u = up * w
    //          - v = w * u   up vector = what is up? top of nose to forehead of viewer. = Camera.up[x,y,z]

    double* O = this->position;
    double w[3];
    double u[3]; // CROSS PRODUCT -> UP x (O-FOCUS)
    double v[3]; // CROSS PRODUCT -> (O-FOCUS) X U
    double t[3];
    for( int i = 0; i < 3; i++){
       w[i] = O[i] - this->focus[i];
       t[i] = 0 - O[i];
    }
    //u[i] = Camera::up X w;
    //v[i] = w X U;
    
    //UPxW
    //AXB = (A.y * B.z - A.z*B.y,
    //       B.x * A.z - A.x*B.z,
    //       A.x * B.y - A.y*B.x)
    
    u[0] = this->up[1] * w[2] - (this->up[2] * w[1]);
    u[1] = w[0] * this->up[2] - (this->up[0] * w[2]);
    u[2] = this->up[0] * w[1] - (this->up[1] * w[0]);

    v[0] = w[1] * u[2] - (w[2] * u[1]);
    v[1] = u[0] * w[2] - (w[0] * u[2]);
    v[2] = w[0] * u[1] - (w[1] * u[0]);
    
    double u_norm = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]); 
    double v_norm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]); 
    double w_norm = sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]); 
    for(int i = 0; i < 3; i++){
        u[i] = u[i] / u_norm;
        v[i] = v[i] / v_norm;
        w[i] = w[i] / w_norm;
    }  

    // Camera frame is now created! (u,v,w,O) 
    // We can express any Cartesian vector (source triangle data coordinates from world space) with some linear combination of the three basis vectors that make up the Camera Frame
    Matrix cameraSpace; 
    cameraSpace.A[0][0] = u[0];
    cameraSpace.A[0][1] = v[0];
    cameraSpace.A[0][2] = w[0];
    cameraSpace.A[0][3] = 0;
    
    cameraSpace.A[1][0] = u[1];
    cameraSpace.A[1][1] = v[1];
    cameraSpace.A[1][2] = w[1];
    cameraSpace.A[1][3] = 0;
    
    cameraSpace.A[2][0] = u[2];
    cameraSpace.A[2][1] = v[2];
    cameraSpace.A[2][2] = w[2];
    cameraSpace.A[2][3] = 0;

    cameraSpace.A[3][0] = u[0]*t[0] + u[1]*t[1] + u[2]*t[2]; // dot product of vector u and vector t (where t is (0,0,0) - O) (t is difference between origin of carditian and camera frame)
    cameraSpace.A[3][1] = v[0]*t[0] + v[1]*t[1] + v[2]*t[2];
    cameraSpace.A[3][2] = w[0]*t[0] + w[1]*t[1] + w[2]*t[2];
    cameraSpace.A[3][3] = 1;
    return cameraSpace;
}
Matrix Camera::ViewTransform(void){
    Matrix viewSpace;
    //forming a 3d cube of what is visible between the camera's near plane and far plane  ... with perspective!
    //VIEW TRANSFORM:
    // Input Parameters: (Angle alpha, near-plane, far-plane) -> ANGLES ALWAYS IN RADIANS!
    // Transforms view frustum to image space cube
    //       - view frustum: bounded by viewing pyramid and near-plane / far-plane
    //           - Near-plane: z = -n and Far-plane: z = -f
    //       - Image Space Cube: -1 <= u,v,w <= 1
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            viewSpace.A[i][j] = 0;
        }
    }
    // cot x = 1 / tan(x)
    viewSpace.A[0][0] = 1/tan(this->angle/2);
    viewSpace.A[1][1] = 1/tan(this->angle/2);
    viewSpace.A[2][2] = (this->far + this->near)/(this->far - this->near);
    viewSpace.A[2][3] = -1;
    viewSpace.A[3][2] = (2*(this->far * this->near)) / (this->far - this->near);
    return viewSpace;
}
Matrix Camera::DeviceTransform(void){
    // Start with a set of coordinates that are in the Image Space... (Result of View Transform)
    Matrix deviceSpace;
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            deviceSpace.A[i][j] = 0;
        }
    }
    deviceSpace.A[0][0] = 1000/2;
    deviceSpace.A[1][1] = 1000/2;
    deviceSpace.A[2][2] = 1;
    deviceSpace.A[3][0] = 1000/2;
    deviceSpace.A[3][1] = 1000/2;
    deviceSpace.A[3][3] = 1;

    // Then we will device transform it to the device space 
    // Image Space -> Device Space  (Image Space range: -1 <= x,y,z <= 1 -> Device Space range: 0 <=x<=width, 0<=y<=height z=z)
    // (x,y,z)     -> (x',y',z')
    // x' = width*(x+1)/2
    // y' = heigth*(y+1)/2
    // z' = z
    
    return deviceSpace;
}

double SineParameterize(int curFrame, int nFrames, int ramp)
{
    int nNonRamp = nFrames-2*ramp;
    double height = 1./(nNonRamp + 4*ramp/M_PI);
    if (curFrame < ramp)
    {
        double factor = 2*height*ramp/M_PI;
        double eval = cos(M_PI/2*((double)curFrame)/ramp);
        return (1.-eval)*factor;
    }
    else if (curFrame > nFrames-ramp)
    {
        int amount_left = nFrames-curFrame;
        double factor = 2*height*ramp/M_PI;
        double eval =cos(M_PI/2*((double)amount_left/ramp));
        return 1. - (1-eval)*factor;
    }
    double amount_in_quad = ((double)curFrame-ramp);
    double quad_part = amount_in_quad*height;
    double curve_part = height*(2*ramp)/M_PI;
    return quad_part+curve_part;
}

Camera
GetCamera(int frame, int nframes)
{
    double t = SineParameterize(frame, nframes, nframes/10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 40*sin(2*M_PI*t);
    c.position[1] = 40*cos(2*M_PI*t);
    c.position[2] = 40;
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    return c;
}

class Triangle
{
  public:
      double         X[3];
      double         Y[3];
      double         Z[3];
      // double representation will be converted back to char upon coloring  pixel
      // first index of 2D array represents vertex # of triangle
      // second index of 2D array represents Red, Green, or Blue
      double         colors[3][3];
      double         normals[3][3];
      // normals[][] is indexed by the vertex first and the dimension second normals[VERTEX][DIMENSIONS]
      // int vertexId = 0
      // int x = 0, y = 1, z = 2
      // normals[vertexId][y] = ...; 
      double         shading[3];
              
      bool           isGoingDown;  
      bool           isGoingUp;
      
      int            leftFlatIndex;
      int            rightFlatIndex;
      int            downOrUpIndex;
};

std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1e_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    cerr << "Done reading" << endl;
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();

    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
    double *color_ptr = var->GetPointer(0);
    //vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    //float *color_ptr = var->GetPointer(0);
    vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
    float *normals = n->GetPointer(0);
    std::vector<Triangle> tris(numTris);
    vtkIdType npts;
    vtkIdType *ptIds;
    int idx;
    for (idx = 0, cells->InitTraversal() ; cells->GetNextCell(npts, ptIds) ; idx++)
    {
        if (npts != 3)
        {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        double *pt = NULL;
        pt = pts->GetPoint(ptIds[0]);
        tris[idx].X[0] = pt[0];
        tris[idx].Y[0] = pt[1];
        tris[idx].Z[0] = pt[2];
#ifdef NORMALS
        tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
        tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
        tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
#endif
        pt = pts->GetPoint(ptIds[1]);
        tris[idx].X[1] = pt[0];
        tris[idx].Y[1] = pt[1];
        tris[idx].Z[1] = pt[2];
#ifdef NORMALS
        tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
        tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
        tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
#endif
        pt = pts->GetPoint(ptIds[2]);
        tris[idx].X[2] = pt[0];
        tris[idx].Y[2] = pt[1];
        tris[idx].Z[2] = pt[2];
#ifdef NORMALS
        tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
        tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
        tris[idx].normals[2][2] = normals[3*ptIds[2]+2];
#endif

        // 1->2 interpolate between light blue, dark blue
        // 2->2.5 interpolate between dark blue, cyan
        // 2.5->3 interpolate between cyan, green
        // 3->3.5 interpolate between green, yellow
        // 3.5->4 interpolate between yellow, orange
        // 4->5 interpolate between orange, brick
        // 5->6 interpolate between brick, salmon
        double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
        double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
        unsigned char RGB[8][3] = { { 71, 71, 219 },
                                    { 0, 0, 91 },
                                    { 0, 255, 255 },
                                    { 0, 128, 0 },
                                    { 255, 255, 0 },
                                    { 255, 96, 0 },
                                    { 107, 0, 0 },
                                    { 224, 76, 76 }
                                  };
        for (int j = 0 ; j < 3 ; j++)
        {
            float val = color_ptr[ptIds[j]];
            int r;
            for (r = 0 ; r < 7 ; r++)
            {
                if (mins[r] <= val && val < maxs[r])
                    break;
            }
            if (r == 7)
            {
                cerr << "Could not interpolate color for " << val << endl;
                exit(EXIT_FAILURE);
            }
            double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
            tris[idx].colors[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
            tris[idx].colors[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
            tris[idx].colors[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
        }
    }

    return tris;
}
double ceil_441(double f)
{
    return ceil(f-0.00001);
}
double floor_441(double f)
{
    return floor(f+0.00001);
}

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void
WriteImage(vtkImageData *img, const char *filename)
{
   std::string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
   writer->SetInputData(img);
   writer->SetFileName(full_filename.c_str());
   writer->Write();
   writer->Delete();
}

class Screen
{
  public:
      unsigned char *buffer;
      double *zbuffer;
      int width, height;
      vtkImageData* image;
};

int GetPixelIndex(int r, int c, Screen screen){
   return 3*((r*screen.width) + c);
} 
int GetZBufferIndex(int r, int c, Screen screen){
   return ((r*screen.width) + c);
} 

int GetRowMax(Triangle t){
   double maxY = t.Y[0];
   for(int i = 1; i < 3; i++){
       maxY = fmax(maxY,t.Y[i]);
   }
   return floor_441(maxY);
}
int GetRowMin(Triangle t){
   double minY = t.Y[0];
   for(int i = 1; i < 3; i++){
       minY = fmin(minY,t.Y[i]);
   }
   return ceil_441(minY);
}

void ImageColor(double red, double green, double blue, Screen screen, int r,int c){
   int pixelIndex = GetPixelIndex(r,c,screen);
   if(pixelIndex < screen.width*screen.height*3 && r >= 0 && c >= 0 && r < screen.height && c < screen.width){
      screen.buffer[pixelIndex+0] = ceil_441(red*255);
      screen.buffer[pixelIndex+1] = ceil_441(green*255);
      screen.buffer[pixelIndex+2] = ceil_441(blue*255);
   }
}

int getMaxYIndexOfTriangle(Triangle t){
    double maxY = fmax(fmax(t.Y[0],t.Y[1]),t.Y[2]);
    if(fabs(maxY - t.Y[0]) < 0.0000001){
        return 0;
    }
    else if(fabs(maxY - t.Y[1]) < 0.0000001){
        return 1;
    }
    return 2;
}
int getMinYIndexOfTriangle(Triangle t){
    double minY = fmin(fmin(t.Y[0],t.Y[1]),t.Y[2]);
    if(fabs(minY - t.Y[0]) < 0.0000001){
        return 0;
    }
    else if(fabs(minY - t.Y[1]) < 0.0000001){
        return 1;
    }
    return 2;
}

void RasterizeArbitraryTriangle(Triangle t, Screen screen);

void ApplyPixelDepth(Triangle t, Screen screen, int r, int c, double z){
   int pixelIndex = GetZBufferIndex(r,c,screen);
   if(pixelIndex < screen.width*screen.height && r >= 0 && c >= 0 && r < screen.height && c < screen.width){
      screen.zbuffer[pixelIndex] = z; 
   }
}
double Interpolate(double A, double B, double fA, double fB, double X){
    // General equation to interpolate:
    // F(X) = F(A) + t*(F(B) - F(A))
    // t = (X-A)/(B-A)
    // these values are defined by what they are in respect to (X or Y?)
    // for leftEnd, rightEnd, it will be in terms of y
    double temp; 
    if(A  > B){
        temp = A;
        A = B;
        B = temp; 
        temp = fA;
        fA = fB;
        fB = temp;
    }
    double t = (X - A) / ((B - A)); 
    if((B - A) < 0.0000001){
        return fA;
    }
    return (fA) + (t * (fB - fA));
}

void RasterizeTriangle(Triangle t, Screen screen){
   if(!t.isGoingDown && !t.isGoingUp){
      RasterizeArbitraryTriangle(t,screen);
      return;
   }
   int rowMin = GetRowMin(t), rowMax = GetRowMax(t);
   int rightFlatIndex = t.rightFlatIndex;
   int leftFlatIndex  = t.leftFlatIndex;
   int downOrUp       =  t.downOrUpIndex;

   for(int r = rowMin; r <= rowMax; r++){
      // Interpolate zbuffer(leftEnd) and zbuffer(rightEnd) from original triangle vertices
      double leftEnd = Interpolate(t.Y[leftFlatIndex],t.Y[downOrUp],t.X[leftFlatIndex],t.X[downOrUp],r);
      double rightEnd = Interpolate(t.Y[rightFlatIndex],t.Y[downOrUp],t.X[rightFlatIndex],t.X[downOrUp],r);

      double leftEndZ = Interpolate(t.Y[leftFlatIndex],t.Y[downOrUp],t.Z[leftFlatIndex],t.Z[downOrUp],r);
      double rightEndZ = Interpolate(t.Y[rightFlatIndex],t.Y[downOrUp],t.Z[rightFlatIndex],t.Z[downOrUp],r);

      double leftEndColorRed = Interpolate(t.Y[leftFlatIndex],t.Y[downOrUp],t.colors[leftFlatIndex][0],t.colors[downOrUp][0], r);
      double leftEndColorGreen = Interpolate(t.Y[leftFlatIndex],t.Y[downOrUp],t.colors[leftFlatIndex][1],t.colors[downOrUp][1], r);
      double leftEndColorBlue = Interpolate(t.Y[leftFlatIndex],t.Y[downOrUp],t.colors[leftFlatIndex][2],t.colors[downOrUp][2], r);

      double rightEndColorRed = Interpolate(t.Y[downOrUp],t.Y[rightFlatIndex],t.colors[downOrUp][0],t.colors[rightFlatIndex][0], r);
      double rightEndColorGreen = Interpolate(t.Y[downOrUp],t.Y[rightFlatIndex],t.colors[downOrUp][1],t.colors[rightFlatIndex][1], r);;
      double rightEndColorBlue = Interpolate(t.Y[downOrUp],t.Y[rightFlatIndex],t.colors[downOrUp][2],t.colors[rightFlatIndex][2], r);

      int rightEndIndex = floor_441(rightEnd);

      if(rightEndIndex > screen.width-1){
         rightEndIndex = screen.width-1;
      }
      for(int c = ceil_441(leftEnd); c <= rightEndIndex; c++){

         double r_c_ColorRed = Interpolate(leftEnd,rightEnd,leftEndColorRed,rightEndColorRed, c);
         double r_c_ColorGreen = Interpolate(leftEnd,rightEnd,leftEndColorGreen,rightEndColorGreen, c);
         double r_c_ColorBlue = Interpolate(leftEnd,rightEnd,leftEndColorBlue,rightEndColorBlue, c);
         // leftEnd and rightEnd are in x-value, so we must give operate the Interpolate function in terms of X-value ( the column)
         double depthFieldValue = Interpolate(leftEnd,rightEnd,leftEndZ,rightEndZ,c);
         if(depthFieldValue > screen.zbuffer[GetZBufferIndex(r,c,screen)]){
             ApplyPixelDepth(t,screen,r,c,depthFieldValue);
             ImageColor(r_c_ColorRed, r_c_ColorGreen,r_c_ColorBlue,screen,r,c);
         }
      }
   }
}

int getArrayIndexOfMaxValue(double x1, double x2, int index1, int index2){
    double maxX = fmax(x1,x2);
    if(fabs(maxX - x1) < 0.000001){  //x1 is max
       return index1;
    }
    return index2;
}

int getArrayIndexOfMinValue(double x1, double x2, int index1, int index2){
    double minX = fmin(x1,x2);
    if(fabs(minX - x1) < 0.000001){  //x1 is max
       return index1;
    }
    return index2;
}

double CalculateShading(){
    return 0.5;
}

void FormNewTriangle(Triangle t, Triangle* new_t){
    int topIndex = getMaxYIndexOfTriangle(t);
    int bottomIndex = getMinYIndexOfTriangle(t);
    int middleIndex = 3 - topIndex - bottomIndex;
    
    double newXvalue = Interpolate(t.Y[bottomIndex],t.Y[topIndex],t.X[bottomIndex],t.X[topIndex],t.Y[middleIndex]);
    double newZvalue = Interpolate(t.Y[bottomIndex],t.Y[topIndex],t.Z[bottomIndex],t.Z[topIndex],t.Y[middleIndex]);

    if(new_t->isGoingDown){
       // new goingDown triangle will include: bottom vertex, middle vertex, and new vertex (middleX,middleY)
       new_t->downOrUpIndex = 0;
       new_t->X[0] = t.X[bottomIndex];
       new_t->Y[0] = t.Y[bottomIndex];
       new_t->Z[0] = t.Z[bottomIndex];
       new_t->colors[0][0] = t.colors[bottomIndex][0];
       new_t->colors[0][1] = t.colors[bottomIndex][1];
       new_t->colors[0][2] = t.colors[bottomIndex][2];
    }
    else{
       // new goingUp triangle will include: top vertex, middle vertex, and new vertex (middleX,middleY)
       new_t->downOrUpIndex = 0;
       new_t->X[0] = t.X[topIndex];
       new_t->Y[0] = t.Y[topIndex];
       new_t->Z[0] = t.Z[topIndex];
       new_t->colors[0][0] = t.colors[topIndex][0];
       new_t->colors[0][1] = t.colors[topIndex][1];
       new_t->colors[0][2] = t.colors[topIndex][2];
    }

    // which X value of the original triangle is greatest? middleIndex or newXvalue?
    new_t->rightFlatIndex = 1; 
    new_t->leftFlatIndex = 2;
    int rightFlatIndex = 1;
    int leftFlatIndex = 2; 

    // PROJECT1F - fake shading... 0.5 for each vertex
    new_t->shading[0] = CalculateShading();
    new_t->shading[1] = CalculateShading();
    new_t->shading[2] = CalculateShading();
  
    
    if(t.X[middleIndex] > newXvalue){
        // middleIndex represents the right flat side of new triangle
        new_t->X[rightFlatIndex] = t.X[middleIndex];
        new_t->Y[rightFlatIndex] = t.Y[middleIndex];
        new_t->Z[rightFlatIndex] = t.Z[middleIndex];
	
        new_t->X[leftFlatIndex] = newXvalue;
        new_t->Y[leftFlatIndex] = t.Y[middleIndex];
        new_t->Z[leftFlatIndex] = newZvalue;

        new_t->colors[rightFlatIndex][0] = t.colors[middleIndex][0];
        new_t->colors[rightFlatIndex][1] = t.colors[middleIndex][1];
        new_t->colors[rightFlatIndex][2] = t.colors[middleIndex][2];

        new_t->colors[leftFlatIndex][0] = Interpolate(t.Y[bottomIndex],t.Y[topIndex],t.colors[bottomIndex][0],t.colors[topIndex][0],t.Y[middleIndex]);
        new_t->colors[leftFlatIndex][1] = Interpolate(t.Y[bottomIndex],t.Y[topIndex],t.colors[bottomIndex][1],t.colors[topIndex][1],t.Y[middleIndex]);
        new_t->colors[leftFlatIndex][2] = Interpolate(t.Y[bottomIndex],t.Y[topIndex],t.colors[bottomIndex][2],t.colors[topIndex][2],t.Y[middleIndex]);
    }
    else{ 
        // middle index represents the left flat side of new triangle  
        new_t->X[leftFlatIndex] = t.X[middleIndex];
        new_t->Y[leftFlatIndex] = t.Y[middleIndex];
        new_t->Z[leftFlatIndex] = t.Z[middleIndex];
	
        new_t->X[rightFlatIndex] = newXvalue;
        new_t->Y[rightFlatIndex] = t.Y[middleIndex];
        new_t->Z[rightFlatIndex] = newZvalue;

        new_t->colors[leftFlatIndex][0] = t.colors[middleIndex][0];
        new_t->colors[leftFlatIndex][1] = t.colors[middleIndex][1];
        new_t->colors[leftFlatIndex][2] = t.colors[middleIndex][2];

        new_t->colors[rightFlatIndex][0] = Interpolate(t.Y[bottomIndex],t.Y[topIndex],t.colors[bottomIndex][0],t.colors[topIndex][0],t.Y[middleIndex]);
        new_t->colors[rightFlatIndex][1] = Interpolate(t.Y[bottomIndex],t.Y[topIndex],t.colors[bottomIndex][1],t.colors[topIndex][1],t.Y[middleIndex]);
        new_t->colors[rightFlatIndex][2] = Interpolate(t.Y[bottomIndex],t.Y[topIndex],t.colors[bottomIndex][2],t.colors[topIndex][2],t.Y[middleIndex]);
    }
}
void RasterizeArbitraryTriangle(Triangle t, Screen screen){

    Triangle goingDown;
    goingDown.isGoingDown = true;
    goingDown.isGoingUp = false;
    FormNewTriangle(t,&goingDown);
    RasterizeTriangle(goingDown,screen);
 
    Triangle goingUp;
    goingUp.isGoingUp = true;
    goingUp.isGoingDown = false;
    FormNewTriangle(t,&goingUp);
    RasterizeTriangle(goingUp,screen);
    return;
}

bool isTriangle(Triangle t){
   if( (fabs(t.X[0] - t.X[1]) < 0.000000001 && fabs(t.Y[0] - t.Y[1]) < 0.0000000001) ||  (fabs(t.X[0] - t.X[2]) < 0.00000000001 && fabs(t.Y[0] - t.Y[2]) < 0.0000000001) || (fabs(t.X[1] - t.X[2]) < 0.00000001 && fabs(t.Y[1] - t.Y[2]) < 0.00000001)){
        return false;
    }
    return true;
}


bool isDoubleEqual(double x, double y){
    if(fabs(x - y) < 0.00000001){
        return true;
    }
    return false;
}

void IdentifyTriangle(Triangle* t, Screen screen){ 
   t->isGoingDown = false;  
   t->isGoingUp = false;
   
   if(!isTriangle(*t)){
      return;
   }
   //if any of following if statements are true, we have up or down triangle
   int v1,v2,v3;
   if(isDoubleEqual(t->Y[0],t->Y[1])){ //Y[0] and Y[1] are flat vertices!
     
      t->leftFlatIndex = getArrayIndexOfMinValue(t->X[0],t->X[1],0,1); //min X value holds the left flat vertex index
      t->rightFlatIndex = getArrayIndexOfMaxValue(t->X[0],t->X[1],0,1);
      t->downOrUpIndex = 2;    

      int indexOfMinYvalue = getArrayIndexOfMinValue(t->Y[0],t->Y[2],0,2); 
      if(indexOfMinYvalue == 2){ //Y[2] is smallest vertex... Y[2] is the bottom vertex! (going down)
         t->isGoingDown = true;
      }
      else{                     //Y[2] is highest vertex .. Y[2] is top vertex! (going up!)
         t->isGoingUp = true;
      }
   }
   else if(isDoubleEqual(t->Y[0],t->Y[2])){ //Y[0] and Y[2] are flat vertices!
      t->leftFlatIndex = getArrayIndexOfMinValue(t->X[0],t->X[2],0,2); //min X value holds the left flat vertex index
      t->rightFlatIndex = getArrayIndexOfMaxValue(t->X[0],t->X[2],0,2);
      t->downOrUpIndex = 1;
      int indexOfMinYvalue = getArrayIndexOfMinValue(t->Y[0],t->Y[1],0,1);
      if(indexOfMinYvalue == 1){  //Y[1] is the smallest vertex... Y[1] is the bottom vertex! (going down)
         t->isGoingDown = true;
      }
      else{
         t->isGoingUp = true;
      }
   }

   else if(isDoubleEqual(t->Y[1],t->Y[2])){ //Y[1] and Y[2] are flat vertices!
      t->leftFlatIndex = getArrayIndexOfMinValue(t->X[1],t->X[2],1,2); //min X value holds the left flat vertex index
      t->rightFlatIndex = getArrayIndexOfMaxValue(t->X[1],t->X[2],1,2);
      t->downOrUpIndex = 0;
      int indexOfMinYvalue = getArrayIndexOfMinValue(t->Y[0],t->Y[2],0,2); 
      if(indexOfMinYvalue == 0){  //Y[0] is the smallest vertex... Y[0] is the bottom vertex! (going down)
         t->isGoingDown = true;
      }
      else{
         t->isGoingUp = true;
      }
   }
   else{
      //is arbitrary
      return;
   }
}

void AllocateScreen(Screen* screen){
   screen->width = 1000;
   screen->height = 1000;
   screen->image = NewImage(screen->width, screen->height);
   screen->buffer = (unsigned char *) (screen->image)->GetScalarPointer(0,0,0);
   int npixels = screen->width * screen->height;
   screen->zbuffer = new double[npixels];
}

void InitializeScreen(Screen* screen){
   int npixels = screen->width * screen->height;
   for (int i = 0 ; i < npixels*3 ; i++){
       screen->buffer[i] = 0; //initialize depth buffer, set all depth values to -1
   }   
   for(int i = 0; i < npixels; i++){
       screen->zbuffer[i] = -1.0;
   }
}

void TransformTrianglesToDeviceSpace(std::vector<Triangle>* triangles, Camera c){
    
    // From Camera Frame, calculate Camera Transform
    // Calculate View Transform
    // Calculate Device Transform
    // Compose Camera Transform, View Transform, and Device Transform into 1 matrix... (M)
    // For each Triangle t, 
    //        apply M to each vertex of t, 
    //                     then apply rasterization / zbuffer 
   
    Matrix cameraSpace = c.CameraTransform();
    Matrix imageSpace = c.ViewTransform();
    Matrix deviceSpace = c.DeviceTransform();
    
    Matrix cameraSpaceToImageSpace = Matrix::ComposeMatrices(cameraSpace,imageSpace);
    Matrix M = Matrix::ComposeMatrices(cameraSpaceToImageSpace,deviceSpace);
   
    for(int i = 0; i < triangles->size(); i++){
        Triangle &t = (*triangles)[i]; //three points (t.X[0],t.Y[0]), (t.X[1],t.Y[1]), (t.X[2],t.Y[1])
        double  vertex1 [] = {t.X[0],t.Y[0],t.Z[0],1};
        double  vertex2 [] = {t.X[1],t.Y[1],t.Z[1],1};
        double  vertex3 [] = {t.X[2],t.Y[2],t.Z[2],1};
        double newV1[4];
        double newV2[4];
        double newV3[4];
        M.TransformPoint(vertex1,newV1);
        M.TransformPoint(vertex2,newV2);
        M.TransformPoint(vertex3,newV3);
        t.X[0] = newV1[0]/newV1[3];
        t.Y[0] = newV1[1]/newV1[3];
        t.Z[0] = newV1[2]/newV1[3];

        t.X[1] = newV2[0]/newV2[3];
        t.Y[1] = newV2[1]/newV2[3];
        t.Z[1] = newV2[2]/newV2[3];

        t.X[2] = newV3[0]/newV3[3];
        t.Y[2] = newV3[1]/newV3[3];
        t.Z[2] = newV3[2]/newV3[3];
    }

}
void RenderTriangles(std::vector<Triangle> triangles,Screen screen){
    for(int i = 0; i < triangles.size(); i++){
       Triangle &t = triangles[i]; 
       IdentifyTriangle(&t,screen);
       RasterizeTriangle(t,screen);
    }
}

void SaveImage(vtkImageData* image, int frameNumber){
   char filename[50];
   if(frameNumber == 0){
       sprintf(filename,"frame00%d",frameNumber);
   }
   else{
       sprintf(filename,"frame%d",frameNumber);
   }
   WriteImage(image,(const char*)filename); 
}

double CalculatePhongShading(LightingParameters &, double *viewDirection, double *normal){
    return 0.5;
}
int main(){
   Screen screen;
   AllocateScreen(&screen);

   for(int i = 0; i < 4; i++){
      std::vector<Triangle> triangles = GetTriangles();
      int f = 250*i;
      InitializeScreen(&screen); 
      Camera c = GetCamera(f, 1000);
      TransformTrianglesToDeviceSpace(&triangles,c); 
      RenderTriangles(triangles,screen);
      SaveImage(screen.image,f);
   }
}
