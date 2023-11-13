#include "component.h"
#include "cimage.h"
#include "bmplib.h"
#include <deque>
#include <iomanip>
#include <iostream>
#include <cmath>


using namespace std;

CImage::CImage(const char* bmp_filename)
{

    img_ = readRGBBMP(bmp_filename, h_, w_);


    if(img_ == NULL) {
        throw std::logic_error("Could not read input file");
    }

    // Set the background RGB color using the upper-left pixel
    for(int i=0; i < 3; i++) {
        bgColor_[i] = img_[0][0][i];
    }


    bfsBgrdThresh_ = 60;


    for (int i = 0; i < h_; i++) {
      std::vector<int> vecOfVectors;
      for (int j = 0; j < w_; j++) {
        vecOfVectors.push_back(-1);
      }
      labels_.push_back(vecOfVectors);
    }

}


CImage::~CImage()
{
    // Add code here if necessary
  if (img_ != nullptr) {
    deallocateImage(img_);
  }

}

bool CImage::isCloseToBground(uint8_t p1[3], double within) {
    // Computes "RGB" (3D Cartesian distance)
    double dist = sqrt(pow(p1[0]-bgColor_[0],2) +
                        pow(p1[1]-bgColor_[1],2) +
                        pow(p1[2]-bgColor_[2],2) );
    return dist <= within;
}

size_t CImage::findComponents()
{
  int compIndex = -1;

  for (int i = 0; i < h_; i++) {
    for (int j = 0; j < w_; j++) {
      if (labels_[i][j] == -1 && !isCloseToBground(img_[i][j], bfsBgrdThresh_)) {
        compIndex++;
        Component newComponent = bfsComponent(i, j, compIndex);
        components_.push_back(newComponent);
      }
    }
  }

  return components_.size();
}


void CImage::printComponents() const
{
    cout << "Height and width of image: " << h_ << "," << w_ << endl;
    cout << setw(4) << "Ord" << setw(4) << "Lbl" << setw(6) << "ULRow" << setw(6) << "ULCol" << setw(4) << "Ht." << setw(4) << "Wi." << endl;
    for(size_t i = 0; i < components_.size(); i++) {
        const Component& c = components_[i];
        cout << setw(4) << i << setw(4) << c.label << setw(6) << c.ulNew.row << setw(6) << c.ulNew.col
             << setw(4) << c.height << setw(4) << c.width << endl;
    }

}


int CImage::getComponentIndex(int mylabel) const
{
  for (unsigned i = 0; i < components_.size(); i++) {
    if (components_[i].label == mylabel) {
      return i;
    }
  }
  return -1;
}


void CImage::translate(int mylabel, int nr, int nc)
{
    // Get the index of specified component
    int cid = getComponentIndex(mylabel);
    if(cid < 0) {
        return;
    }
    int h = components_[cid].height;
    int w = components_[cid].width;

    int totalNC = nc + w;
    int totalNR = nr + h;

    if (nr < 0 || nc < 0 || totalNR > h_ || totalNC > w_ || totalNC > 0 || totalNR > 0 || nr > h_ || nr > w_) {
      return;
    }



    Location nl(nr, nc);
    components_[cid].ulNew = nl;
}


void CImage::forward(int mylabel, int delta)
{
    int cid = getComponentIndex(mylabel);
    if(cid < 0 || delta <= 0) {
        return;
    }


    Component compToMove = components_[cid];
    components_.erase(components_.begin() + cid);

    int newPos = delta + cid;
    if (newPos <= components_.size()) {
      components_.insert(components_.begin() + newPos, compToMove);
    }
    else {
      components_.push_back(compToMove);
    }

}


void CImage::backward(int mylabel, int delta)
{
    int cid = getComponentIndex(mylabel);
    if(cid < 0 || delta <= 0) {
        return;
    }

  Component compToMove = components_[cid];
  components_.erase(components_.begin() + cid);

  int newPos = cid - delta;
  if (newPos >= 0) {
    components_.insert(components_.begin()+ newPos, compToMove);
  }
  else {
    components_.insert(components_.begin(), compToMove);
  }

}

void CImage::save(const char* filename)
{
    // Create another image filled in with the background color
    uint8_t*** out = newImage(bgColor_);

  
    for (unsigned int j = 0; j < components_.size(); j++) {
      int deltaRow = components_[j].ulNew.row - components_[j].ulOrig.row;
      int deltaCol = components_[j].ulNew.col - components_[j].ulOrig.col;
      for (int k = 0; k < h_; k++) {
        int newRow = 0;
        int newCol = 0;
        for (int i = 0; i < w_; i++) {
          if (labels_[k][i] == components_[j].label) {
            newRow = k + deltaRow;
            newCol = i + deltaCol;
            out[newRow][newCol] = img_[j][k];
          }
        }
      }
    }
  
    writeRGBBMP(filename, out, h_, w_);

}


uint8_t*** CImage::newImage(uint8_t bground[3]) const
{
    uint8_t*** img = new uint8_t**[h_];
    for(int r=0; r < h_; r++) {
        img[r] = new uint8_t*[w_];
        for(int c = 0; c < w_; c++) {
            img[r][c] = new uint8_t[3];
            img[r][c][0] = bground[0];
            img[r][c][1] = bground[1];
            img[r][c][2] = bground[2];
        }
    }
    return img;
}

void CImage::deallocateImage(uint8_t*** img) const
{

  if (img != nullptr) {
    for (int i = 0; i < h_; i++) {
      for (int j = 0; j < w_; j++) {
        delete [] img[i][j];
      }
    delete [] img[i];
    }
  delete [] img;
  } 
}


Component CImage::bfsComponent(int pr, int pc, int mylabel)
{
    // Arrays to help produce neighbors easily in a loop
    // by encoding the **change** to the current location.
    // Goes in order N, NW, W, SW, S, SE, E, NE
    int neighbor_row[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
    int neighbor_col[8] = {0, -1, -1, -1, 0, 1, 1, 1};

    int right = pc;
    int bot = pr;
    int left = pc;
    int top = pr;

    deque<Location> locQueue = deque<Location>();
    locQueue.push_back(Location(pr, pc));
    labels_[pr][pc] = mylabel;
    
    while (!locQueue.empty()) {
      Location currentLoc = locQueue.front();
      locQueue.pop_front();

      right = max(right, currentLoc.col);
      bot = max(bot, currentLoc.row);
      left = min(left, currentLoc.col);
      top = min(top, currentLoc.row);

      for (int i = 0; i < 8; i++) {
        Location neighLoc = Location(currentLoc.row + neighbor_row[i], currentLoc.col + neighbor_col[i]);
        int newRow = neighLoc.row;
        int newCol = neighLoc.col;

        if (newRow < 0 || newCol >= h_ || newCol < 0 || newCol >= w_) {
          continue;
        }
         
        bool unLabeled = (labels_[newRow][newCol] == -1);
        bool closeToBackground = !isCloseToBground(img_[newRow][newCol], bfsBgrdThresh_);

        if (unLabeled && closeToBackground) {
          labels_[newRow][newCol] = mylabel;
          locQueue.push_back(Location(newRow, newCol));
        }
      }
    }
  return Component(Location(top, left), bot - top + 1, right - left + 1, mylabel);
}

void CImage::labelToRGB (const char* filename)
{
    //multiple ways to do this -- this is one way
    vector<uint8_t[3]> colors(components_.size());
    for(unsigned int i=0; i<components_.size(); i++) {
        colors[i][0] = rand() % 256;
        colors[i][1] = rand() % 256;
        colors[i][2] = rand() % 256;
    }

    for(int i=0; i<h_; i++) {
        for(int j=0; j<w_; j++) {
            int mylabel = labels_[i][j];
            if(mylabel >= 0) {
                img_[i][j][0] =  colors[mylabel][0];
                img_[i][j][1] =  colors[mylabel][1];
                img_[i][j][2] =  colors[mylabel][2];
            } else {
                img_[i][j][0] = 0;
                img_[i][j][1] = 0;
                img_[i][j][2] = 0;
            }
        }
    }
    writeRGBBMP(filename, img_, h_, w_);
}


const Component& CImage::getComponent(size_t i) const
{
    if(i >= components_.size()) {
        throw std::out_of_range("Index to getComponent is out of range");
    }
    return components_[i];
}


size_t CImage::numComponents() const
{
    return components_.size();
}


void CImage::drawBoundingBoxesAndSave(const char* filename)
{
    for(size_t i=0; i < components_.size(); i++){
        Location ul = components_[i].ulOrig;
        int h = components_[i].height;
        int w = components_[i].width;
        for(int i = ul.row; i < ul.row + h; i++){
            for(int k = 0; k < 3; k++){
                img_[i][ul.col][k] = 255-bgColor_[k];
                img_[i][ul.col+w-1][k] = 255-bgColor_[k];

            }
            // cout << "bb " << i << " " << ul.col << " and " << i << " " << ul.col+w-1 << endl; 
        }
        for(int j = ul.col; j < ul.col + w ; j++){
            for(int k = 0; k < 3; k++){
                img_[ul.row][j][k] = 255-bgColor_[k];
                img_[ul.row+h-1][j][k] = 255-bgColor_[k];

            }
            // cout << "bb2 " << ul.row << " " << j << " and " << ul.row+h-1 << " " << j << endl; 
        }
    }
    writeRGBBMP(filename, img_, h_, w_);
}


