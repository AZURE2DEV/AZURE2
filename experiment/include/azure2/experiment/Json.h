#ifndef AZURE2_EXPERIMENT_JSON_H
#define AZURE2_EXPERIMENT_JSON_H

#include <string>
#include <vector>

namespace azure2 {
namespace experiment {

/*!
 * A minimal JSON writer for the schema (Schema::ToJson, docs/upgrade/project_format3.md section 6.1) and the
 * diagnostics (section 6.2): objects keep the order in which keys are written, so the output is deterministic and the
 * committed copy tests/schema/schema_v3.json can be compared byte for byte.  Writing only; pyazr parses with Python's
 * json module.  Numbers are written with the fewest digits (up to 17) that read back to the same double.
 */
class JsonWriter {
 public:
  /// indent > 0: one member per line, nested by `indent` spaces; 0: everything on one line.
  explicit JsonWriter(int indent = 2) : indent_(indent) {}

  JsonWriter &BeginObject();
  JsonWriter &EndObject();
  JsonWriter &BeginArray();
  JsonWriter &EndArray();
  JsonWriter &Key(const std::string &key);
  JsonWriter &String(const std::string &value);
  JsonWriter &Number(double value);
  JsonWriter &Integer(long value);
  JsonWriter &Bool(bool value);
  JsonWriter &Null();
  /// An array of strings, on one line.
  JsonWriter &StringArray(const std::vector<std::string> &values);

  const std::string &Str() const { return out_; }

 private:
  void BeforeValue();
  void NewLine();
  struct Level {
    bool object;
    int count;
  };
  int indent_;
  std::vector<Level> stack_;
  bool afterKey_ = false;
  std::string out_;
};

std::string JsonEscape(const std::string &s);
std::string JsonNumber(double value);

}  // namespace experiment
}  // namespace azure2

#endif
