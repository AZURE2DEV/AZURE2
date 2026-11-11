#ifndef SHAREDTABLE_H
#define SHAREDTABLE_H

#include <atomic>
#include <cassert>
#include <memory>
#include <utility>

/*!
 * A read-only table shared between copies of a point (WorkingCopy clones, mapped points, re-placed sub-points), written in
 * place only by a sole holder: the copy-on-write handle of Unit's energy-dependent tables and Legendre polynomials.
 *
 * The lifetime is that of a std::shared_ptr.  The ownership test is not shared_ptr::use_count(): that is a relaxed
 * load, and an acquire fence after it, which makes it correct under the C++ memory model, is not modelled by
 * ThreadSanitizer (docs/upgrade/sanitizer_findings.md §4.2).  Each table instead carries its own count of handles,
 * `holders`:
 *
 *  - copying a handle increments it (relaxed, as shared_ptr does: the source handle is held, so the count is ≥ 1 and
 *    nothing is published by the increment);
 *  - dropping a handle (destructor, reset, assignment) decrements it with release order, after every access the
 *    holder made to the table;
 *  - SoleOwner() loads it with acquire order.  A load that reads 1 reads the last read-modify-write in the count's
 *    modification order; every earlier release decrement heads a release sequence that continues through all later
 *    read-modify-writes ([intro.races]), so the load synchronizes with each of them, and every other holder's reads
 *    of the table happen before the sole holder's write.
 *
 * As with use_count(), the test is exact only while no other thread copies a handle to the same table; within a
 * parallel step a point copies only from a point its own thread has finished (Unit::CopyEDependentValuesFrom).  A
 * Weak reference (NodeTableCache) can create a holder at any moment; a table published through one is never written
 * again (NodeTableCache.h).
 */
template <class T>
class SharedTable {
  struct Box {
    std::atomic<long> holders{1};
    T value;
    template <class... A>
    explicit Box(A &&...a) : value(std::forward<A>(a)...) {}
  };

 public:
  /// A non-owning reference; lock() gives a holding handle, or an empty one if the table is gone.
  class Weak {
   public:
    Weak() = default;
    explicit Weak(const SharedTable &t) : w_(t.p_) {}
    SharedTable lock() const {
      std::shared_ptr<Box> p = w_.lock();
      if (p) p->holders.fetch_add(1, std::memory_order_relaxed);
      return SharedTable(std::move(p));
    }
    bool expired() const { return w_.expired(); }

   private:
    std::weak_ptr<Box> w_;
  };

  SharedTable() = default;
  template <class... A>
  static SharedTable Make(A &&...a) { return SharedTable(std::make_shared<Box>(std::forward<A>(a)...)); }
  SharedTable(const SharedTable &o) : p_(o.p_) {
    if (p_) p_->holders.fetch_add(1, std::memory_order_relaxed);
  }
  SharedTable(SharedTable &&o) noexcept : p_(std::move(o.p_)) {}
  /// Copy and move assignment (the argument is a copy or the moved-from handle); the old table is dropped.
  SharedTable &operator=(SharedTable o) noexcept {
    std::swap(p_, o.p_);
    return *this;
  }
  ~SharedTable() { Drop(); }

  void reset() { Drop(); }
  explicit operator bool() const { return (bool)p_; }
  const T &operator*() const { return p_->value; }
  const T *operator->() const { return &p_->value; }

  /// True if this handle is the only holder of its table, which it may then write (Mutable).
  bool SoleOwner() const { return p_ && p_->holders.load(std::memory_order_acquire) == 1; }
  /// The table, writable; only after SoleOwner() returned true on this thread.
  T &Mutable() {
    assert(SoleOwner());
    return p_->value;
  }

 private:
  explicit SharedTable(std::shared_ptr<Box> p) : p_(std::move(p)) {}  // one holder, already counted
  void Drop() {
    if (!p_) return;
    p_->holders.fetch_sub(1, std::memory_order_release);  // after this holder's last access to the table
    p_.reset();
  }
  std::shared_ptr<Box> p_;
};

#endif
