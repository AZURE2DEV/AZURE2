#ifndef AZUREOUTPUT_H
#define AZUREOUTPUT_H

#include "AZUREFBuffer.h"
#include <vector>

/// A class to assist in writing AZURE output files

/*!
 * The EData::WriteOutputFiles function simply loops over all ESegment and EPoint objects when
 * writing the output of an AZURE calculation. To ensure that all output for a given entrance and exit
 * pair combination is written to a single file, the AZUREOutput class is used.  The AZUREOutput object
 * is a container for a vector of AZUREFBuffer objects.
 */

class AZUREOutput {
 public:
  /// Write output files into this directory.
  AZUREOutput(std::string);
  ~AZUREOutput();
  /// Are these extrapolation outputs (AZUREOut_*.extrap) rather than data ones?
  bool IsExtrap() const;
  /// Buffer for one reaction's output file, opening it on first use.
  std::filebuf *operator()(int entranceKey, int exitKey, bool isAngDist = false);
  /// Number of open output buffers.
  int NumAZUREFBuffers() const;
  /// 1-based position of the buffer for this reaction, or 0 if not open.
  int IsAZUREFBuffer(int, int, bool);
  /// Directory the files are written to.
  std::string GetOutputDir() const;
  /// Append a buffer.
  void AddAZUREFBuffer(AZUREFBuffer *);
  /// Mark the output as extrapolation, changing the file suffix.
  void SetExtrap();
  /// Buffer \p i, 1-based.
  AZUREFBuffer *GetAZUREFBuffer(int);

 private:
  bool is_extrap_;
  std::string outputdir_;
  std::vector<AZUREFBuffer *> azurefbuffers_;
};

#endif
