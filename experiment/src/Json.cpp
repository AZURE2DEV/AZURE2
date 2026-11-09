// A minimal deterministic JSON writer (Json.h).
#include "azure2/experiment/Json.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace azure2 {
namespace experiment {

std::string JsonEscape(const std::string &s) {
  std::string out = "\"";
  for (unsigned char c : s) {
    switch (c) {
      case '"': out += "\\\""; break;
      case '\\': out += "\\\\"; break;
      case '\n': out += "\\n"; break;
      case '\t': out += "\\t"; break;
      case '\r': out += "\\r"; break;
      default:
        if (c < 0x20) {
          char buf[8];
          std::snprintf(buf, sizeof buf, "\\u%04x", c);
          out += buf;
        } else {
          out += static_cast<char>(c);  // UTF-8 bytes pass through
        }
    }
  }
  return out + "\"";
}

std::string JsonNumber(double value) {
  if (!std::isfinite(value)) return "null";
  char buf[40];
  for (int precision = 1; precision <= 17; precision++) {
    std::snprintf(buf, sizeof buf, "%.*g", precision, value);
    if (std::strtod(buf, nullptr) == value) break;
  }
  return buf;
}

void JsonWriter::NewLine() {
  if (indent_ <= 0) return;
  out_ += "\n";
  out_ += std::string(stack_.size() * static_cast<std::size_t>(indent_), ' ');
}

void JsonWriter::BeforeValue() {
  if (afterKey_) {
    afterKey_ = false;
    return;
  }
  if (!stack_.empty()) {
    if (stack_.back().count++ > 0) out_ += indent_ <= 0 ? ", " : ",";
    NewLine();
  }
}

JsonWriter &JsonWriter::BeginObject() {
  BeforeValue();
  out_ += "{";
  stack_.push_back({true, 0});
  return *this;
}

JsonWriter &JsonWriter::EndObject() {
  const bool empty = stack_.back().count == 0;
  stack_.pop_back();
  if (!empty) NewLine();
  out_ += "}";
  return *this;
}

JsonWriter &JsonWriter::BeginArray() {
  BeforeValue();
  out_ += "[";
  stack_.push_back({false, 0});
  return *this;
}

JsonWriter &JsonWriter::EndArray() {
  const bool empty = stack_.back().count == 0;
  stack_.pop_back();
  if (!empty) NewLine();
  out_ += "]";
  return *this;
}

JsonWriter &JsonWriter::Key(const std::string &key) {
  BeforeValue();
  out_ += JsonEscape(key) + ": ";
  afterKey_ = true;
  return *this;
}

JsonWriter &JsonWriter::String(const std::string &value) {
  BeforeValue();
  out_ += JsonEscape(value);
  return *this;
}

JsonWriter &JsonWriter::Number(double value) {
  BeforeValue();
  out_ += JsonNumber(value);
  return *this;
}

JsonWriter &JsonWriter::Integer(long value) {
  BeforeValue();
  out_ += std::to_string(value);
  return *this;
}

JsonWriter &JsonWriter::Bool(bool value) {
  BeforeValue();
  out_ += value ? "true" : "false";
  return *this;
}

JsonWriter &JsonWriter::Null() {
  BeforeValue();
  out_ += "null";
  return *this;
}

JsonWriter &JsonWriter::StringArray(const std::vector<std::string> &values) {
  BeforeValue();
  out_ += "[";
  for (std::size_t i = 0; i < values.size(); i++) out_ += (i ? ", " : "") + JsonEscape(values[i]);
  out_ += "]";
  return *this;
}

}  // namespace experiment
}  // namespace azure2
