#ifndef MTQueue_h
#define MTQueue_h

#include <atomic>
#include <cassert>
#include <condition_variable>
#include <mutex>
#include <type_traits>
#include <utility>

class semaphore
{
public:
  explicit semaphore(unsigned int count = 0) noexcept : m_count(count) {}

  void post() noexcept
  {
    {
      std::unique_lock<std::mutex> lock(m_mutex);
      ++m_count;
    }
    m_cv.notify_one();
  }

  void post(unsigned int count) noexcept
  {
    {
      std::unique_lock<std::mutex> lock(m_mutex);
      m_count += count;
    }
    m_cv.notify_all();
  }

  void wait() noexcept
  {
    std::unique_lock<std::mutex> lock(m_mutex);
    m_cv.wait(lock, [&]() { return m_count != 0; });
    --m_count;
  }

  template <typename T>
  bool wait_for(T&& t) noexcept
  {
    std::unique_lock<std::mutex> lock(m_mutex);
    if(!m_cv.wait_for(lock, t, [&]() { return m_count != 0; }))
      return false;
    --m_count;
    return true;
  }

  template <typename T>
  bool wait_until(T&& t) noexcept
  {
    std::unique_lock<std::mutex> lock(m_mutex);
    if(!m_cv.wait_until(lock, t, [&]() { return m_count != 0; }))
      return false;
    --m_count;
    return true;
  }

private:
  unsigned int m_count;
  std::mutex m_mutex;
  std::condition_variable m_cv;
};

class fast_semaphore
{
public:
  explicit fast_semaphore(int count = 0) noexcept : m_count(count), m_semaphore(0) { assert(count > -1); }

  void post() noexcept
  {
    int count = m_count.fetch_add(1, std::memory_order_release);
    if(count < 0)
      m_semaphore.post();
  }

  void wait() noexcept
  {
    int count = m_count.fetch_sub(1, std::memory_order_acquire);
    if(count < 1)
      m_semaphore.wait();
  }

private:
  std::atomic_int m_count;
  semaphore m_semaphore;
};

template <typename T>
class blocking_queue
{
public:
  blocking_queue(unsigned int size)
      : m_size(size), m_pushIndex(0), m_popIndex(0), m_count(0), m_data((T*)operator new(size * sizeof(T))),
        m_openSlots(size), m_fullSlots(0)
  {
  }

  blocking_queue(const blocking_queue&) = delete;
  blocking_queue(blocking_queue&&)      = delete;
  blocking_queue& operator=(const blocking_queue&) = delete;
  blocking_queue& operator=(blocking_queue&&) = delete;

  ~blocking_queue()
  {
    while(m_count--)
      {
        m_data[m_popIndex].~T();
        m_popIndex = ++m_popIndex % m_size;
      }
    operator delete(m_data);
  }

  void push(const T& item)
  {
    m_openSlots.wait();
    {
      std::lock_guard<std::mutex> lock(m_cs);
      new(m_data + m_pushIndex) T(item);
      m_pushIndex = ++m_pushIndex % m_size;
      ++m_count;
    }
    m_fullSlots.post();
  }

  void pop(T& item)
  {
    m_fullSlots.wait();
    {
      std::lock_guard<std::mutex> lock(m_cs);
      item = m_data[m_popIndex];
      m_data[m_popIndex].~T();
      m_popIndex = ++m_popIndex % m_size;
      --m_count;
    }
    m_openSlots.post();
  }

  bool empty()
  {
    std::lock_guard<std::mutex> lock(m_cs);
    return m_count == 0;
  }

private:
  unsigned int m_size;
  unsigned int m_pushIndex;
  unsigned int m_popIndex;
  unsigned int m_count;
  T* m_data;

  semaphore m_openSlots;
  semaphore m_fullSlots;
  std::mutex m_cs;
};

template <typename T>
class fast_blocking_queue
{
public:
  explicit fast_blocking_queue(unsigned int size)
    : m_size(size), m_pushIndex(0), m_popIndex(0), m_count(0),
      m_data((T*)operator new(size * sizeof(T))), m_openSlots(size), m_fullSlots(0)
  {
    assert(size != 0);
  }

  ~fast_blocking_queue() noexcept
  {
    while(m_count--)
      {
        m_data[m_popIndex].~T();
        m_popIndex = ++m_popIndex % m_size;
      }
    operator delete(m_data);
  }

  template <typename Q = T>
  typename std::enable_if<std::is_copy_constructible<Q>::value and std::is_nothrow_copy_constructible<Q>::value,
                          void>::type
  push(const T& item) noexcept
  {
    m_openSlots.wait();

    auto pushIndex = m_pushIndex.fetch_add(1);
    new(m_data + (pushIndex % m_size)) T(item);
    ++m_count;

    auto expected = m_pushIndex.load();
    while(!m_pushIndex.compare_exchange_strong(expected, m_pushIndex % m_size))
      expected = m_pushIndex.load();

    m_fullSlots.post();
  }

  template <typename Q = T>
  typename std::enable_if<std::is_move_constructible<Q>::value and std::is_nothrow_move_constructible<Q>::value,
                          void>::type
  push(T&& item) noexcept
  {
    m_openSlots.wait();

    auto pushIndex = m_pushIndex.fetch_add(1);
    new(m_data + (pushIndex % m_size)) T(std::move(item));
    ++m_count;

    auto expected = m_pushIndex.load();
    while(!m_pushIndex.compare_exchange_strong(expected, m_pushIndex % m_size))
      expected = m_pushIndex.load();

    m_fullSlots.post();
  }

  template <typename Q = T>
  typename std::enable_if<not std::is_move_assignable<Q>::value and std::is_nothrow_copy_assignable<Q>::value,
                          void>::type
  pop(T& item) noexcept
  {
    m_fullSlots.wait();

    auto popIndex = m_popIndex.fetch_add(1);
    item          = m_data[popIndex % m_size];
    m_data[popIndex % m_size].~T();
    --m_count;

    auto expected = m_popIndex.load();
    while(!m_popIndex.compare_exchange_strong(expected, m_popIndex % m_size))
      expected = m_popIndex.load();

    m_openSlots.post();
  }

  template <typename Q = T>
  typename std::enable_if<std::is_move_assignable<Q>::value and std::is_nothrow_move_assignable<Q>::value, void>::type
  pop(T& item) noexcept
  {
    m_fullSlots.wait();

    auto popIndex = m_popIndex.fetch_add(1);
    item          = std::move(m_data[popIndex % m_size]);
    m_data[popIndex % m_size].~T();
    --m_count;

    auto expected = m_popIndex.load();
    while(!m_popIndex.compare_exchange_strong(expected, m_popIndex % m_size))
      expected = m_popIndex.load();

    m_openSlots.post();
  }

  T pop() noexcept(std::__is_nothrow_invocable<void, decltype(&fast_blocking_queue<T>::pop<T>), T&>::value)
  {
    T item;
    pop(item);
    return item;
  }

  bool empty() { return m_count == 0; }

private:
  const unsigned int m_size;
  std::atomic_uint m_pushIndex;
  std::atomic_uint m_popIndex;
  std::atomic_uint m_count;
  T* m_data;

  fast_semaphore m_openSlots;
  fast_semaphore m_fullSlots;
};

#endif
