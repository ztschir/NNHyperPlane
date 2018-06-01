#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <cmath>
#include <algorithm>
#include <cuda_profiler_api.h>

//#define USE_GPU

#define Real float
#define Point Real

struct Node{
  bool isLeaf;
  Node *left;
  Node *right;
  Point *leafPoints;
  int numberOfPoints;
};


void generatePoints(Point * points, int numberOfPoints, int numberOfDimensions){
  //srand (time(NULL));
  srand(0);
  for(int i = 0; i < numberOfPoints; i++)
    for(int j = 0; j < numberOfDimensions; j++)
      points[i*numberOfDimensions + j] = ((double) rand() / (RAND_MAX));
}

/*__global__
void findCentroidDevice(Point * points, Point * center, int numberOfPoints, int numberOfDimensions){
  int i = blockIdx.x;
  if(i < numberOfDimensions){
    Real sum = 0.0;
    //for(int j = 0; j < numberOfPoints; j++){
    //sum += points[j].values[i];
      //}
    sum += points[0].values[i];
    center->values[i] = sum / numberOfPoints;
  }
}
*/
Point * findCentroid(Point * points, int numberOfPoints, int numberOfDimensions){
  Point * center = new Point[numberOfDimensions];
  
  #ifdef USE_GPU
  Point * devicePoints;
  Point * deviceCenter;
  cudaMalloc((void **)&devicePoints, numberOfPoints*sizeof(Point));
  cudaMalloc((void **)&deviceCenter, 1*sizeof(Point));

  cudaMemcpy(devicePoints, points, numberOfPoints*sizeof(Point), cudaMemcpyHostToDevice);
  findCentroidDevice<<<numberOfDimensions, 1>>>(devicePoints, deviceCenter, numberOfPoints, numberOfDimensions);
  cudaMemcpy(center, deviceCenter, 1*sizeof(Point), cudaMemcpyDeviceToHost);

  cudaFree(devicePoints);
  cudaFree(deviceCenter);  
  
  #else
  
  for(int i = 0; i < numberOfDimensions; i++){
    Real sum = 0.0;
    for(int j = 0; j < numberOfPoints; j++)
      sum += points[j*numberOfDimensions + i];
    center[i] = sum / numberOfPoints;
  }
  #endif

  return center;
}

Real findSingleDistance(Point * firstPoint, Point * secondPoint, int numberOfDimensions){
  Real sum = 0.0;
  for(int i = 0; i < numberOfDimensions; i++){
    sum += (firstPoint[i] - secondPoint[i]) * (firstPoint[i] - secondPoint[i]);
  }
  return sqrt(sum);
}

Real * findDistance(Point * firstPoints, Point * diffPoint, int numberOfPoints, int numberOfDimensions){
  Real * distances = new Real[numberOfPoints];
  for(int i = 0; i < numberOfPoints; i++)
    distances[i] = findSingleDistance(&firstPoints[numberOfDimensions*i], diffPoint, numberOfDimensions);
  return distances;
}


int maxPointIndex(Real * distances, int numberOfPoints){
  Real max = 0.0;
  int ind = 0;
  for(int i = 0; i < numberOfPoints; i++){
    if(distances[i] > max){
      max = distances[i];
      ind = i;
    }
  }
  return ind;
}

Point * pointDiff(Point * point1, Point * point2, int numberOfDimensions){
  Point * diff = new Point[numberOfDimensions];
  
  for(int i = 0; i < numberOfDimensions; i++)
    diff[i] = point1[i] - point2[i];

  return diff;
}

Real * projectPoints(Point * points, Point * point2, Point * point1, int numberOfPoints, int numberOfDimensions){
  Real * distance1 = findDistance(points, point1, numberOfPoints, numberOfDimensions);
  Real * distance2 = findDistance(points, point2, numberOfPoints, numberOfDimensions);
  Real * distance = findDistance(point1, point2, 1, numberOfDimensions);

  Real * t = new Real[numberOfPoints];
  for(int i = 0; i < numberOfPoints; i++){
    t[i] = ((distance1[i]*distance1[i] - distance2[i]*distance2[i]) / distance[0]) + distance[0];
    t[i] /= 2.0;
    t[i] /= distance[0];
  }
  return t;
}

template <typename Iterator>
Real median(Iterator begin, Iterator end) {
  Iterator middle = begin + (end - begin) / 2;
   
  std::nth_element(begin, middle, end);
  
  if ((end - begin) % 2 != 0) {
    return *middle;
  } else {
    Iterator lower_middle = std::max_element(begin, middle);
    return (*middle + *lower_middle) / 2.0;
  }
}

Node * hyperplaneSplit(Point * points, int numberOfPoints, int numberOfDimensions, int minNumberOfPoints){
  
  if(numberOfPoints/2 < minNumberOfPoints){
    Node * leaf = new Node();
    leaf->isLeaf = true;
    leaf->leafPoints = points;
    leaf->numberOfPoints = numberOfPoints;
    return leaf;
  }
  
  Node * node = new Node();
  node->isLeaf = false;
  
  Point * center = findCentroid(points, numberOfPoints, numberOfDimensions);
  Real * distance = findDistance(points, center, numberOfPoints, numberOfDimensions);
  Point * Point1 = &points[maxPointIndex(distance, numberOfPoints)*numberOfDimensions];

  distance = findDistance(points, Point1, numberOfPoints, numberOfDimensions);
  Point * Point2 = &points[maxPointIndex(distance, numberOfPoints)*numberOfDimensions];
  Point * diffPoint = pointDiff(Point1, Point2, numberOfDimensions);
  Real * alpha = projectPoints(points, diffPoint, center, numberOfPoints, numberOfDimensions);

  Real * tempAlpha = new Real[numberOfPoints];
  std::copy(alpha, alpha + numberOfPoints, tempAlpha);
  Real mid = median(tempAlpha, tempAlpha + numberOfPoints);

  int numberOfPointsLeft = (int)std::floor((float)numberOfPoints/2.0);
  int numberOfPointsRight = numberOfPoints - numberOfPointsLeft;
  int j = 0;
  Real * temp = new Real[numberOfDimensions];
  for(int i = 0; i < numberOfPointsLeft; i++){
    if(alpha[i] > mid){
      for(j; j < numberOfPointsRight; j++){
        if(alpha[numberOfPointsLeft+j] < mid){
          memcpy(temp, &points[i*numberOfDimensions], numberOfDimensions*sizeof(Real));
          memcpy(&points[i*numberOfDimensions], &points[(numberOfPointsLeft+j)*numberOfDimensions],
                 numberOfDimensions*sizeof(Real));
          memcpy(&points[(numberOfPointsLeft+j)*numberOfDimensions], temp, numberOfDimensions*sizeof(Real));
          j++;
          break;
        }
      }
    }
  }

  node->left = hyperplaneSplit(points, numberOfPointsLeft, numberOfDimensions, minNumberOfPoints);
  node->right = hyperplaneSplit(&points[numberOfPointsLeft*numberOfDimensions], numberOfPointsRight,
                                numberOfDimensions, minNumberOfPoints);
  
  delete temp;
  delete center;
  delete tempAlpha;
  return node;
}


void traverseTree(Node * node, int numberOfDimensions){
  if(!node->isLeaf){
    traverseTree(node->left, numberOfDimensions);    
    traverseTree(node->right, numberOfDimensions);
  }
  else{
    printf("\n");
    for(int i = 0; i < node->numberOfPoints * numberOfDimensions; i += 2){
      printf("%f %f\n", i, node->leafPoints[i], node->leafPoints[i+1]);
    }
  }
}

int main(int argc, char **argv){

  int numberOfPoints = 16;
  int numberOfDimensions = 2;
  
  //int numberOfNeighbors = 4;
  int minNumberOfPoints = 4;
  
  Point points[numberOfPoints*numberOfDimensions];
  generatePoints(points, numberOfPoints, numberOfDimensions);

  Node * topLevelNode = hyperplaneSplit(points, numberOfPoints, numberOfDimensions, minNumberOfPoints);

  traverseTree(topLevelNode, numberOfDimensions);


  cudaProfilerStop();
  
  return 0;
}
