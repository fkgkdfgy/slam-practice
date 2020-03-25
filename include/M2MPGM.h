#ifndef _M2MPGM_H_
#define _M2MPGM_H_
#include "common.h"
#include "PGM.h"


// Probability map of different pixelSize 
// whereas the pixelSize in the PGM is the smallest one
class M2MPGM
{
public:

    M2MPGM() {}

    M2MPGM(shared_ptr<PGMBuilder> pgmbuilder_,int height = 90);

    M2MPGM(const PGMBuilder& pgmbuilder_, int height = 90);
    // TEST PASS
    vector<float> ManyPixelPGMCreate(int height,const vector<uint16_t> & Origin_table_ ,
                                     const vector<float> & Corres_table_);
    
    // SlideWindow Function
    void AddValue(float temp);
    void RemoveValue(float temp);
    float GetMaxProb();

    // Get Part Function
    bool & GetFinished() { return CreatedFinished; }
    const vector<float>& GetMap(int height) { assert(height < (PGMset_.size())); return PGMset_[height]; }
    int  GetHeight() { return int(PGMset_.size() - 1); }
    const MapLimits& GetMapLimits()const { return maplimits_; }
    const AlignedBox2i& GetScale() const{ return scale_; }
    double GetPixelSize() const { return pixelSize_; }
private:
    double pixelSize_;
    vector<vector<float>> PGMset_;
    MapLimits maplimits_;
    int height_;
    bool CreatedFinished;
    deque<float> SlideWindow;
    AlignedBox2i scale_;
    static int max_height;
};

#endif
