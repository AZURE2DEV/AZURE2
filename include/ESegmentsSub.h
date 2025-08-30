#ifndef ESEGMENTSSUB_H
#define ESEGMENTSSUB_H

class CNuc;
class Config;
class EData;

class ESegmentsSub {
public:
    ESegmentsSub(int entranceKey, int exitKey);
    ~ESegmentsSub();
    
    int GetEntranceKey() const { return entranceKey_; }
    int GetExitKey() const { return exitKey_; }
    
    double CalculateTheoretical(int pointIndex, CNuc* cnuc, const Config& configure, EData* edata) const;
    
private:
    int entranceKey_;
    int exitKey_;
};

#endif