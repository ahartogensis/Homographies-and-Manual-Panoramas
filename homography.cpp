#include "homography.h"
#include "matrix.h"

using namespace std;

void applyHomography(const Image &source, const Matrix &H, Image &out,
                     bool bilinear) {
  // --------- HANDOUT  PS06 ------------------------------
  // Transform image source using the homography H, and composite in onto out.
  // if bilinear == true, using bilinear interpolation. Use nearest neighbor
  // otherwise.

  Matrix hInverse = H.inverse();

  for(int x = 0; x < out.width(); x++){
    for(int y = 0; y < out.height(); y++){
      for(int z = 0; z < out.channels(); z++){
          Matrix homography = Matrix(3,1);
          homography(0,0) = x;
          homography(1,0) = y; 
          homography(2,0) = 1; 

          Matrix result = hInverse * homography; 

          Matrix projective_coordinates = Matrix(2, 1); 
          projective_coordinates(0,0) = result(0,0) / result(2,0);
          projective_coordinates(1,0) = result(1,0) / result(2,0);

          float value = 0; 
          if(projective_coordinates(0,0) >= 0 && projective_coordinates(0,0) < source.width() - 1 
              && projective_coordinates(1,0) >= 0 && projective_coordinates(1,0) < source.height() - 1){
            if(bilinear){
              value = interpolateLin(source, projective_coordinates(0,0), projective_coordinates(1,0), z, false);
            }else{
              value = source(round(projective_coordinates(0,0)), round(projective_coordinates(1,0)), z);
            }
          }else{
            value = out(x,y,z);
          }
         
          out(x,y,z) = value; 

      }
    }
  }

}

Matrix computeHomography(const CorrespondencePair correspondences[4]) {
  // --------- HANDOUT  PS06 ------------------------------
  // Compute a homography from 4 point correspondences.
  Matrix A = Matrix::Zero(9,9); 
  Matrix B = Matrix::Zero(9,1); 
  Matrix homography = Matrix(3,3); 

  for(int i = 0; i < 4; i++){
    float x = correspondences[i].point1.x(); 
    float y = correspondences[i].point1.y();

    float n_x = correspondences[i].point2.x(); 
    float n_y = correspondences[i].point2.y();

    //assume i = 1
    A(2*i, 0) = x;
    A(2*i, 1) = y;
    A(2*i, 2) = 1;
    A(2*i, 3) = 0;
    A(2*i, 4) = 0;
    A(2*i, 5) = 0;
    A(2*i, 6) = -n_x*x;
    A(2*i, 7) = -n_x*y;
    A(2*i, 8) = -n_x;

    A(2*i + 1, 0) = 0;
    A(2*i + 1, 1) = 0;
    A(2*i + 1, 2) = 0;
    A(2*i + 1, 3) = x;
    A(2*i + 1, 4) = y;
    A(2*i + 1, 5) = 1;
    A(2*i + 1, 6) = -n_y*x;
    A(2*i + 1, 7) = -n_y*y;
    A(2*i + 1, 8) = -n_y;
  }
  A(8,8) = 1;
  B(8,0) = 1; 

  Matrix T = A.inverse() * B;

  homography(0,0) = T(0,0);
  homography(0,1) = T(1,0);
  homography(0,2) = T(2,0);
  homography(1,0) = T(3,0);
  homography(1,1) = T(4,0);
  homography(1,2) = T(5,0);
  homography(2,0) = T(6,0);
  homography(2,1) = T(7,0);
  homography(2,2) = T(8,0);

  return homography;
}

BoundingBox computeTransformedBBox(int imwidth, int imheight, Matrix H) {
  // --------- HANDOUT  PS06 ------------------------------
  // Predict the bounding boxes that encompasses all the transformed
  // coordinates for pixels frow and Image with size (imwidth, imheight)

  //corners of the picture 
  Vec3f upper_left(0,0,1);
  Vec3f lower_left(0,imheight,1);
  Vec3f upper_right(imwidth,0,1);
  Vec3f lower_right(imwidth,imheight,1);
  
  //new coordinaters
  Vec3f new_upperLeft = H * upper_left; 
  Vec3f new_lowerLeft = H * lower_left; 
  Vec3f new_upperRight = H * upper_right; 
  Vec3f new_lowerRight = H * lower_right; 

  //projection coordinates 
  Vec2f proj_upperleft;
  Vec2f proj_lowerleft; 
  Vec2f proj_upperright; 
  Vec2f proj_lowerright; 

  proj_upperleft.x() = new_upperLeft.x() / new_upperLeft.z();
  proj_upperleft.y() = new_upperLeft.y() / new_upperLeft.z();

  proj_lowerleft.x() = new_lowerLeft.x() / new_lowerLeft.z();
  proj_lowerleft.y() = new_lowerLeft.y() / new_lowerLeft.z(); 

  proj_upperright.x() = new_upperRight.x() / new_upperRight.z();
  proj_upperright.y() = new_upperRight.y() / new_upperRight.z();

  proj_lowerright.x() = new_lowerRight.x() / new_lowerRight.z();
  proj_lowerright.y() = new_lowerRight.y() / new_lowerRight.z();

  return BoundingBox(min(proj_lowerleft.x(), proj_upperleft.x()), max(proj_lowerright.x(), proj_upperright.x()),
                    min(proj_upperright.y(), proj_upperleft.y()), max(proj_lowerright.y(), proj_lowerleft.y()));
}

BoundingBox bboxUnion(BoundingBox B1, BoundingBox B2) {
  // --------- HANDOUT  PS06 ------------------------------
  // Compute the bounding box that tightly bounds the union of B1 an B2.
  return BoundingBox(min(B1.x1, B2.x1), max(B1.x2, B2.x2), min(B1.y1, B2.y1), max(B1.y2, B2.y2));
}

Matrix makeTranslation(BoundingBox B) {
  // --------- HANDOUT  PS06 ------------------------------
  // Compute a translation matrix (as a homography matrix) that translates the
  // top-left corner of B to (0,0).
  Matrix output = Matrix::Zero(3,3);
  output << 1, 0, -B.x1, 0, 1, -B.y1, 0, 0, 1;  
  return output;
}

Image stitch(const Image &im1, const Image &im2,
             const CorrespondencePair correspondences[4]) {
  // --------- HANDOUT  PS06 ------------------------------
  // Transform im1 to align with im2 according to the set of correspondences.
  // make sure the union of the bounding boxes for im2 and transformed_im1 is
  // translated properly (use makeTranslation)

  Matrix H = computeHomography(correspondences); 
  BoundingBox B1 = computeTransformedBBox(im1.width(), im1.height(), H);
  Matrix translation1 = makeTranslation(B1); 
  BoundingBox B2 (0, im2.width() - 1, 0, im2.height() - 1);
  BoundingBox B = bboxUnion(B1, B2);
  Matrix translation = makeTranslation(B); 
  Image BlackBox(B.x2 - B.x1, B.y2 - B.y1, im1.channels());
  BlackBox.set_color(0.0); 

  applyHomography(im1, translation * H, BlackBox, true);
  applyHomography(im2, translation, BlackBox, true); 

  return BlackBox;
}

// debug-useful
Image drawBoundingBox(const Image &im, BoundingBox bbox) {
  // --------- HANDOUT  PS06 ------------------------------
  /*
    ________________________________________
   / Draw me a bounding box!                \
   |                                        |
   | "I jumped to my                        |
   | feet, completely thunderstruck. I      |
   | blinked my eyes hard. I looked         |
   | carefully all around me. And I saw a   |
   | most extraordinary small person, who   |
   | stood there examining me with great    |
   | seriousness."                          |
   \              Antoine de Saint-Exupery  /
    ----------------------------------------
           \   ^__^
            \  (oo)\_______
               (__)\       )\/\
                   ||----w |
                   ||     ||
  */
  Image output = im; 
  for (int i = bbox.x1; i < bbox.x2; i++){
        output(i,bbox.y1,0) = 0;
        output(i,bbox.y1,1) = 1;
        output(i,bbox.y1,2) = 0;  

        output(i,bbox.y2,0) = 0;
        output(i,bbox.y2,1) = 1;
        output(i,bbox.y2,2) = 0;  
    }
  for (int j = bbox.y1; j < bbox.y2; j++){
        output(bbox.x1,j,0) = 0;
        output(bbox.x1,j,1) = 1;
        output(bbox.x1,j,2) = 0;  

        output(bbox.x2,j,0) = 0;
        output(bbox.x2,j,1) = 1;
        output(bbox.x2,j,2) = 0;   
    }
  return output;
}

void applyHomographyFast(const Image &source, const Matrix &H, Image &out,
                         bool bilinear) {
  // --------- HANDOUT  PS06 ------------------------------
  // Same as apply but change only the pixels of out that are within the
  // predicted bounding box (when H maps source to its new position).
  BoundingBox B = computeTransformedBBox(source.width(), source.height(), H); 
  Matrix hInverse = H.inverse();

  for(int x = B.x1; x < B.x2; x++){
    for(int y = B.y1; y < B.y2; y++){
      for(int z = 0; z < out.channels(); z++){
          Matrix homography = Matrix(3,1);
          homography(0,0) = x;
          homography(1,0) = y; 
          homography(2,0) = 1; 

          Matrix result = hInverse * homography; 

          Matrix projective_coordinates = Matrix(2, 1); 
          projective_coordinates(0,0) = result(0,0) / result(2,0);
          projective_coordinates(1,0) = result(1,0) / result(2,0);

          float value = 0; 
          if(projective_coordinates(0,0) >= 0 && projective_coordinates(0,0) < source.width() 
              && projective_coordinates(1,0) >= 0 && projective_coordinates(1,0) < source.height()){
            if(bilinear){
              value = interpolateLin(source, projective_coordinates(0,0), projective_coordinates(1,0), z, false);
            }else{
              value = source(round(projective_coordinates(0,0)), round(projective_coordinates(1,0)), z);
            }
          }else{
            value = out(x,y,z);
          }
         
          out(x,y,z) = value; 

      }
    }
  }
}
