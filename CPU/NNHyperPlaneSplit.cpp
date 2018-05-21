#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <cmath>
#include <algorithm>

#define Real float

struct Point{
  Real *values;
  Point *neighbors;
};

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
  for(int i = 0; i < numberOfPoints; i++){
    points[i].values = new Real[numberOfDimensions];
    for(int j = 0; j < numberOfDimensions; j++)
      points[i].values[j] = ((double) rand() / (RAND_MAX));
  }
}

Point * findCentroid(Point * points, int numberOfPoints, int numberOfDimensions){
  Point * center = new Point();
  center->values = new Real[numberOfDimensions]();
  for(int i = 0; i < numberOfDimensions; i++){
    for(int j = 0; j < numberOfPoints; j++)
      center->values[i] += points[j].values[i];
    center->values[i] /= numberOfPoints;
  }
  return center;
}

Real findSingleDistance(Point * firstPoint, Point * secondPoint, int numberOfDimensions){
  Real sum = 0.0;
  for(int i = 0; i < numberOfDimensions; i++){
    sum += (firstPoint->values[i] - secondPoint->values[i]) * (firstPoint->values[i] - secondPoint->values[i]);
  }
  return sqrt(sum);
}

Real * findDistance(Point * firstPoints, Point * secondPoint, int numberOfFirstPoints, int numberOfDimensions){
  Real * distances = new Real[numberOfFirstPoints];
  for(int i = 0; i < numberOfFirstPoints; i++)
    distances[i] = findSingleDistance(firstPoints + i, secondPoint, numberOfDimensions);
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
  Point * diff = new Point();
  diff->values = new Real[numberOfDimensions]();
  
  for(int i = 0; i < numberOfDimensions; i++)
    diff->values[i] = point1->values[i] - point2->values[i];

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

/*int partitions(Real * alpha, int low,int high){
  int p=low;
  int r=high;
  Real x=alpha[r];
  int i=p-1;
  for(int j=p; j<=r-1; j++){
    if (alpha[j] <= x)
      std::swap(alpha[++i],alpha[j]);
  }
  std::swap(alpha[i+1],alpha[r]);
  return i+1;
}
Real median(Real * alpha, int left,int right,int kth){
  for(;;){
    int pivotIndex=partitions(alpha,left,right);
    int len=pivotIndex-left+1;
    if(kth == len)
      return alpha[pivotIndex];
    else if(kth < len)
      right=pivotIndex-1;
    else{
      kth=kth-len;
      left=pivotIndex+1;
    }
  }
  }*/
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
    leaf->leafPoints = new Point[numberOfPoints];
    for(int i = 0; i < numberOfPoints; i++)
      leaf->leafPoints[i] = points[i];
    leaf->numberOfPoints = numberOfPoints;
    return leaf;
  }

  Node * node = new Node();
  node->isLeaf = false;
  
  Point * center = findCentroid(points, numberOfPoints, numberOfDimensions);
  Real * distance = findDistance(points, center, numberOfPoints, numberOfDimensions);
  
  int maxPoint = maxPointIndex(distance, numberOfPoints);
  Point * Point1 = &points[maxPointIndex(distance, numberOfPoints)];

  distance = findDistance(points, Point1, numberOfPoints, numberOfDimensions);
  Point * Point2 = &points[maxPointIndex(distance, numberOfPoints)];
  Point * diffPoint = pointDiff(Point1, Point2, numberOfDimensions);
  Real * alpha = projectPoints(points, diffPoint, center, numberOfPoints, numberOfDimensions);

  Real * tempAlpha = new Real[numberOfPoints];
  std::copy(alpha, alpha + numberOfPoints, tempAlpha);
  Real mid = median(tempAlpha + 0, tempAlpha + numberOfPoints);
  delete tempAlpha;

  int numberOfPointsLeft = (int)std::floor((float)numberOfPoints/2.0);
  int numberOfPointsRight = (int)std::ceil((float)numberOfPoints/2.0);  
  
  Point * pointsLeft = new Point[numberOfPointsLeft];
  Point * pointsRight = new Point[numberOfPointsRight];

  for(int i = 0, j = 0, k = 0; i < numberOfPoints; i++){
    if(alpha[i] < mid)
      pointsLeft[j++] = points[i];
    else
      pointsRight[k++] = points[i];
  }

  node->left = hyperplaneSplit(pointsLeft, numberOfPointsLeft, numberOfDimensions, minNumberOfPoints);
  node->right = hyperplaneSplit(pointsRight, numberOfPointsRight, numberOfDimensions, minNumberOfPoints);

  delete pointsLeft;
  delete pointsRight;
  
  return node;
}


void traverseTree(Node * node){
  if(!node->isLeaf){
    traverseTree(node->left);    
    traverseTree(node->right);
  }
  else{
    printf("\n");
    for(int i = 0; i < node->numberOfPoints; i++){
      printf("%f %f\n", i, node->leafPoints[i].values[0], node->leafPoints[i].values[1]);
    }
  }
}

int main(int argc, char **argv){

  int numberOfPoints = 16;
  int numberOfDimensions = 2;
  int numberOfNeighbors = 4;
  int minNumberOfPoints = 4;
  
  Point points[numberOfPoints];
  generatePoints(points, numberOfPoints, numberOfDimensions);

  Node * topLevelNode = hyperplaneSplit(points, numberOfPoints, numberOfDimensions, minNumberOfPoints);

  traverseTree(topLevelNode);
  
  return 0;
}
